/**
 * A scene a job builds on request: one of its plugin methods, called through
 * object_method, returns ``{success: true, scene}`` (or ``{success: false,
 * error}``). Report links name one (report MoorhenLink); the project Moorhen
 * page takes it as ``?sceneJob=&sceneMethod=&sceneArgs=`` and applies what it
 * returns, so reloading the page builds the scene again from what the job
 * holds now.
 *
 * Only methods named ``...Scene`` are called: a link is a URL anyone can be
 * sent, and those methods build a scene and change nothing.
 */

export interface SceneSource {
  jobId: number;
  method: string;
  kwargs: Record<string, unknown>;
}

export const isSceneMethod = (method: string) => /^[A-Za-z_][A-Za-z0-9_]*Scene$/.test(method);

export function sceneSourceUrl(projectId: number, source: SceneSource): string | null {
  if (!Number.isFinite(projectId) || !isSceneMethod(source.method)) return null;
  const query = new URLSearchParams({
    sceneJob: String(source.jobId),
    sceneMethod: source.method,
    sceneArgs: JSON.stringify(source.kwargs ?? {}),
  });
  return `/ccp4i2/moorhen-page/project/${projectId}?${query.toString()}`;
}

/** The scene source a page URL names, or null when it names none (or a
 *  method that is not a scene method). */
export function sceneSourceFromParams(params: {
  get(name: string): string | null;
}): SceneSource | null {
  const jobId = parseInt(params.get("sceneJob") || "", 10);
  const method = params.get("sceneMethod") || "";
  if (!Number.isFinite(jobId) || !isSceneMethod(method)) return null;
  let kwargs: Record<string, unknown> = {};
  try {
    const parsed = JSON.parse(params.get("sceneArgs") || "{}");
    if (parsed && typeof parsed === "object" && !Array.isArray(parsed)) kwargs = parsed;
  } catch {
    return null;
  }
  return { jobId, method, kwargs };
}
