/**
 * Client half of staged import for web deployments.
 *
 * When the server advertises `import_staging` (a served deployment with a staging
 * directory) and a file is larger than the threshold, deliver it into the
 * staging directory in chunks that stay under every body cap, then import it by
 * the returned owner-bound handle. Desktop (local_path) and small files are
 * untouched.
 *
 * Only the transport lives here; the call sites (`uploadFileParam`,
 * `import-project-content`) decide when to use it.
 */

import { getAccessToken } from "@ccp4/ccp4i2-api";

const PROXY_BASE = "/api/proxy/ccp4i2/";

/**
 * The bearer header every call needs, the way api-fetch.ts attaches it.
 * The proxy route answers 401 to any non-public path without one, and a
 * served deployment has no other way to authenticate a bare fetch.
 */
async function authHeaders(extra: Record<string, string> = {}): Promise<Record<string, string>> {
  const token = await getAccessToken();
  return token ? { ...extra, Authorization: `Bearer ${token}` } : extra;
}

export interface StagingCapability {
  chunk_bytes: number;
  max_bytes: number;
  threshold_bytes: number;
}

let capabilityPromise: Promise<StagingCapability | null> | null = null;

/** The deployment's staging capability, or null. Cached for the session. */
export function stagingCapability(): Promise<StagingCapability | null> {
  if (!capabilityPromise) {
    capabilityPromise = authHeaders()
      .then((headers) => fetch(`${PROXY_BASE}version/`, { headers }))
      .then((r) => (r.ok ? r.json() : null))
      .then((j) => (j && j.import_staging) || null)
      .catch(() => null);
  }
  return capabilityPromise;
}

export interface StageOptions {
  /** 0..1 as chunks complete. */
  onProgress?: (fraction: number) => void;
  signal?: AbortSignal;
  /** Chunks in flight at once. */
  concurrency?: number;
}

/** A human message for a staging status code, falling back to the server's. */
async function stagingError(resp: Response): Promise<Error> {
  let message = "";
  try {
    const body = await resp.json();
    message = body?.error || "";
  } catch {
    /* no JSON body */
  }
  const byStatus: Record<number, string> = {
    413: "The file is larger than this server allows.",
    429: "Too many uploads are in progress; please try again shortly.",
    410: "The upload session expired; please start the import again.",
    409: "Some chunks did not arrive; please try the import again.",
    422: "The uploaded file did not verify; please try the import again.",
  };
  return new Error(message || byStatus[resp.status] || `Upload failed (${resp.status}).`);
}

/**
 * Stage a file into the server's staging directory in chunks and return the
 * handle (`upload_id`) to import it by. Throws a human-readable error on any
 * refusal.
 */
export async function stageFile(
  file: File,
  cap: StagingCapability,
  opts: StageOptions = {},
): Promise<string> {
  const { onProgress, signal, concurrency = 3 } = opts;

  const beginResp = await fetch(`${PROXY_BASE}staged-uploads/`, {
    method: "POST",
    headers: await authHeaders({ "Content-Type": "application/json" }),
    body: JSON.stringify({ filename: file.name, size_bytes: file.size }),
    signal,
  });
  if (!beginResp.ok) throw await stagingError(beginResp);
  const { upload_id, chunk_bytes } = await beginResp.json();

  const nChunks = Math.max(1, Math.ceil(file.size / chunk_bytes));
  let done = 0;

  const putChunk = async (index: number): Promise<void> => {
    const start = index * chunk_bytes;
    const blob = file.slice(start, Math.min(start + chunk_bytes, file.size));
    for (let attempt = 1; ; attempt++) {
      try {
        const r = await fetch(
          `${PROXY_BASE}staged-uploads/${upload_id}/chunks/${index}/`,
          {
            method: "PUT",
            body: blob,
            headers: await authHeaders({ "Content-Type": "application/octet-stream" }),
            signal,
          },
        );
        if (!r.ok) throw await stagingError(r);
        done += 1;
        onProgress?.(done / nChunks);
        return;
      } catch (err) {
        // A 4xx from the server won't get better on retry; a transient network
        // drop might. Retry a few times, then give up.
        if (signal?.aborted || attempt >= 3 || err instanceof Error && /\b(413|409|422|410|404)\b/.test(err.message)) {
          throw err;
        }
        await new Promise((res) => setTimeout(res, 400 * attempt));
      }
    }
  };

  // Bounded-concurrency pool over the chunk indexes.
  const queue = Array.from({ length: nChunks }, (_, i) => i);
  const worker = async () => {
    let index: number | undefined;
    while ((index = queue.shift()) !== undefined) {
      await putChunk(index);
    }
  };
  await Promise.all(
    Array.from({ length: Math.min(concurrency, nChunks) }, worker),
  );

  const finResp = await fetch(`${PROXY_BASE}staged-uploads/${upload_id}/finish/`, {
    method: "POST",
    headers: await authHeaders(),
    signal,
  });
  if (!finResp.ok) throw await stagingError(finResp);
  return upload_id as string;
}

/**
 * The transport a call site should use for one file, or null to send bytes.
 * Returns "staged_upload" (with the handle) when a served deployment advertises
 * staging and the file is over the threshold; otherwise null (upload as before).
 * Desktop `local_path` is decided by the caller before this.
 */
export async function maybeStage(
  file: File,
  opts: StageOptions = {},
): Promise<{ field: "staged_upload"; value: string } | null> {
  const cap = await stagingCapability();
  if (!cap || file.size <= cap.threshold_bytes) return null;
  const handle = await stageFile(file, cap, opts);
  return { field: "staged_upload", value: handle };
}
