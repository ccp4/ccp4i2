/**
 * Capability detection + transport for integrated natural-language scene
 * generation (design: MOORHEN_SCENES_NL_UI.md).
 *
 * `useSceneNlCapability()` probes the deployment's `nlp/status` endpoint and
 * reports whether an in-app "Generate" affordance should appear. It is a
 * PROVIDER ABSTRACTION, not a Materia hard-dependency: today the only provider
 * is Materia's compounds `nlp_scene`; future desktop-byo / on-device providers
 * plug in behind the same shape. Desktop Electron (no compounds app) 404s the
 * probe → `available: false` → the panel falls back to copy-paste. That absence
 * IS the signal; nothing needs building into desktop.
 *
 * `generateScene()` is the transport for the POST. It uses a raw fetch (not the
 * apiFetch wrapper) so it can read the structured `{status:"error", kind}` body
 * on a non-2xx response — the wrapper throws on non-ok and discards the body.
 */
import useSWR from "swr";
import { getAccessToken } from "@ccp4/ccp4i2-api";

// Same-origin proxy paths. NO trailing slash — the Next proxy adds the one
// Django wants; a pre-added slash 308s the body away (see the compounds proxy
// convention).
const STATUS_URL = "/api/proxy/compounds/nlp/status";
const SCENE_URL = "/api/proxy/compounds/nlp/scene";

export type SceneNlProvider = "none" | "materia";

export interface SceneNlCapability {
  /** True when an in-app Generate affordance should be shown. */
  available: boolean;
  provider: SceneNlProvider;
  /** What generate_scene will actually use, so the UI can anticipate JSON vs
   *  YAML. Ingest strips nulls regardless, so this is informational. */
  responseFormat?: "strict" | "raw";
  model?: string;
  quota?: { limit?: number; remaining?: number };
}

interface NlStatusResponse {
  scene?: {
    enabled?: boolean;
    responseFormat?: "strict" | "raw";
    model?: string;
    dailyLimit?: number;
    dailyRemaining?: number;
  };
}

/** Materia's error kinds (mirrors nlp_query), plus a client-only "network". */
export type SceneGenerateErrorKind =
  | "disabled"
  | "rate_limited"
  | "bad_request"
  | "llm_error"
  | "network";

export class SceneGenerateError extends Error {
  constructor(
    public readonly kind: SceneGenerateErrorKind,
    message: string,
  ) {
    super(message);
    this.name = "SceneGenerateError";
  }
}

async function authHeaders(): Promise<Record<string, string>> {
  const headers: Record<string, string> = { "Content-Type": "application/json" };
  const token = await getAccessToken();
  if (token) headers["Authorization"] = `Bearer ${token}`;
  return headers;
}

async function statusFetcher(): Promise<NlStatusResponse | null> {
  try {
    const resp = await fetch(STATUS_URL, { headers: await authHeaders() });
    if (!resp.ok) return null; // 404 (no compounds app, e.g. desktop) → unavailable
    return (await resp.json()) as NlStatusResponse;
  } catch {
    return null; // network / offline → unavailable (safe: falls back to copy-paste)
  }
}

/**
 * Probe capability. SWR-cached; a 404/unreachable result is terminal for the
 * session (don't hammer a missing endpoint), so retry-on-error is off.
 */
export function useSceneNlCapability(): SceneNlCapability {
  const { data } = useSWR<NlStatusResponse | null>(STATUS_URL, statusFetcher, {
    revalidateOnFocus: false,
    shouldRetryOnError: false,
    dedupingInterval: 5 * 60_000,
  });
  const scene = data?.scene;
  const available = !!scene?.enabled;
  return {
    available,
    provider: available ? "materia" : "none",
    responseFormat: scene?.responseFormat,
    model: scene?.model,
    quota: scene
      ? { limit: scene.dailyLimit, remaining: scene.dailyRemaining }
      : undefined,
  };
}

const STATUS_TO_KIND: Record<number, SceneGenerateErrorKind> = {
  400: "bad_request",
  404: "disabled",
  429: "rate_limited",
  502: "llm_error",
};

/**
 * POST a request (+ optional grounding) to the scene endpoint and return the
 * model's RAW scene text (YAML or JSON) for the caller to normalise/parse.
 * Throws SceneGenerateError with the server's `kind` (or one mapped from the
 * HTTP status) on failure.
 */
export async function generateScene(
  request: string,
  grounding?: string,
): Promise<string> {
  let resp: Response;
  try {
    resp = await fetch(SCENE_URL, {
      method: "POST",
      headers: await authHeaders(),
      body: JSON.stringify({ request, grounding }),
    });
  } catch (e) {
    throw new SceneGenerateError(
      "network",
      `Could not reach the scene service: ${e instanceof Error ? e.message : "network error"}`,
    );
  }

  let body: { status?: string; scene?: string; kind?: SceneGenerateErrorKind; message?: string } = {};
  try {
    body = await resp.json();
  } catch {
    /* non-JSON body; handled below via status */
  }

  if (resp.ok && body.status === "scene" && typeof body.scene === "string") {
    return body.scene;
  }

  const kind: SceneGenerateErrorKind =
    body.kind ?? STATUS_TO_KIND[resp.status] ?? "llm_error";
  throw new SceneGenerateError(
    kind,
    body.message || `Scene generation failed (HTTP ${resp.status}).`,
  );
}
