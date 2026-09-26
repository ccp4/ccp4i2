/**
 * The run targets this deployment registered, from GET version/.
 *
 * The server is the authority on where a program can run
 * (docs/run-target-dispatch.md): a client offers to dispatch only when some
 * registered target runs programs, and names the choice from this list.
 * Cached for the session, like the staging capability.
 */
import { PROXY_BASE, authHeaders } from "./staged-upload";

export interface RunTarget {
  name: string;
  runs_jobs: boolean;
  runs_programs: boolean;
  error?: string;
}

let promise: Promise<RunTarget[]> | null = null;

export function runTargets(): Promise<RunTarget[]> {
  if (!promise) {
    promise = authHeaders()
      .then((headers) => fetch(`${PROXY_BASE}version/`, { headers }))
      .then((r) => (r.ok ? r.json() : null))
      .then((j) => (j && Array.isArray(j.run_targets) ? (j.run_targets as RunTarget[]) : []))
      .catch(() => []);
  }
  return promise;
}

/** The registered targets that can run a program (axis B). */
export function programTargets(): Promise<RunTarget[]> {
  return runTargets().then((t) => t.filter((x) => x.runs_programs && !x.error));
}
