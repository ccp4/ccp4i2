/**
 * The recorded Moorhen session, from the window's side.
 *
 * A session window represents one job. This hook talks to the four
 * `interactive_*` endpoints: it reads the session state and load plan,
 * heartbeats while open (informational: no timer acts on it), saves models
 * into the job's drop directory, finishes the session, and on close tells
 * the server a window detached (which finishes an empty session and
 * otherwise leaves it open for reconnect). See docs/moorhen-task-design.md.
 */
import { useCallback, useEffect, useMemo, useRef, useState } from "react";
import { moorhen } from "moorhen/types/moorhen";
import { apiGet, apiPost, apiUpload } from "../api-fetch";
import { withAccessToken } from "../api";
import { serialiseMolecule } from "../lib/moorhen-serialise";
import { JobStatus } from "../types/models";

export interface SessionLoadPlanEntry {
  kind: "dictionary" | "coordinates" | "map_2fofc" | "map_fofc" | "map_anom" | "map";
  param: string;
  file_id: number | null;
  file_uuid: string | null;
  name: string;
  type: string | null;
  sub_type: number | null;
  label: string;
}

export interface SessionOutput {
  number: number;
  name: string;
  kind: "model" | "dictionary";
  annotation: string;
}

export interface SessionRow {
  requested_at: string;
  last_heartbeat: string | null;
  attached: boolean;
  dispatched: boolean;
  finished: boolean;
  finished_at: string | null;
}

export interface SessionState {
  job_id: number;
  task_name: string;
  status: number;
  session: SessionRow | null;
  load_plan: SessionLoadPlanEntry[];
  outputs: SessionOutput[];
}

export type FinishDisposition =
  | "dispatched"
  | "harvesting"
  | "deleted"
  | "kept_open"
  | "already_finished";

export interface SessionJobInfo {
  id: number;
  number: string;
  title: string;
  project: number;
  task_name: string;
}

export interface MoorhenSessionApi {
  jobId: number;
  job: SessionJobInfo | null;
  state: SessionState | null;
  error: string | null;
  /** The session exists, is not finished, and the job is still RUNNING. */
  isOpen: boolean;
  refresh: () => Promise<void>;
  saveMolecule: (mol: moorhen.Molecule, annotation: string) => Promise<SessionOutput>;
  finish: () => Promise<FinishDisposition>;
}

const HEARTBEAT_MS = 10_000;

/** The interactive endpoints answer in the {success, data} envelope. */
function unwrap<T>(payload: any): T {
  if (payload && typeof payload === "object" && "success" in payload) {
    if (payload.success === false) {
      throw new Error(payload.error || "request failed");
    }
    return payload.data as T;
  }
  return payload as T;
}

export function useMoorhenSession(
  jobId: number | null | undefined,
  enabled: boolean,
): MoorhenSessionApi | null {
  const [state, setState] = useState<SessionState | null>(null);
  const [job, setJob] = useState<SessionJobInfo | null>(null);
  const [error, setError] = useState<string | null>(null);
  // The close beacon cannot await a token, so the URL it needs is kept
  // resolved while the window is open (refreshed on every heartbeat).
  const beaconUrlRef = useRef<string | null>(null);
  const finishedRef = useRef(false);
  const active = enabled && !!jobId;

  const refresh = useCallback(async () => {
    if (!active) return;
    try {
      const next = unwrap<SessionState>(
        await apiGet(`jobs/${jobId}/interactive_session/`),
      );
      setState(next);
      setError(null);
      finishedRef.current = !!next.session?.finished;
    } catch (err) {
      setError(err instanceof Error ? err.message : String(err));
    }
  }, [active, jobId]);

  // Job identity for the panel header, once.
  useEffect(() => {
    if (!active) return;
    let cancelled = false;
    apiGet(`jobs/${jobId}`)
      .then((info: any) => {
        if (!cancelled && info?.id) {
          setJob({
            id: info.id,
            number: info.number,
            title: info.title,
            project: info.project,
            task_name: info.task_name,
          });
        }
      })
      .catch(() => {});
    return () => {
      cancelled = true;
    };
  }, [active, jobId]);

  useEffect(() => {
    if (active) void refresh();
  }, [active, refresh]);

  const isOpen =
    !!state?.session && !state.session.finished && state.status === JobStatus.RUNNING;

  // Heartbeat while open. Informational only, by design: it lets the job
  // menu say whether a window is attached. It also keeps the beacon URL
  // fresh.
  useEffect(() => {
    if (!active || !isOpen) return;
    let stopped = false;
    const beat = async () => {
      if (stopped) return;
      try {
        beaconUrlRef.current = await withAccessToken(
          `/api/proxy/ccp4i2/jobs/${jobId}/interactive_finish/`,
        );
        const next = unwrap<SessionState>(
          await apiPost(`jobs/${jobId}/interactive_heartbeat/`, {}),
        );
        if (!stopped) setState(next);
      } catch (err) {
        // A heartbeat that fails because the session was finished or
        // cancelled elsewhere is the cue to re-read the state.
        if (!stopped) void refresh();
      }
    };
    void beat();
    const timer = window.setInterval(beat, HEARTBEAT_MS);
    return () => {
      stopped = true;
      window.clearInterval(timer);
    };
  }, [active, isOpen, jobId, refresh]);

  // Window closed without Finish: say so. The server finishes an empty
  // session (nothing saved: mark the job for deletion) and otherwise keeps
  // it open for reconnect. Best-effort; a lost beacon leaves the session
  // open, which is the safe outcome.
  useEffect(() => {
    if (!active) return;
    const onPageHide = () => {
      if (finishedRef.current) return;
      const url = beaconUrlRef.current;
      if (!url || typeof navigator.sendBeacon !== "function") return;
      navigator.sendBeacon(
        url,
        new Blob([JSON.stringify({ finished: false })], { type: "application/json" }),
      );
    };
    window.addEventListener("pagehide", onPageHide);
    return () => window.removeEventListener("pagehide", onPageHide);
  }, [active]);

  const saveMolecule = useCallback(
    async (mol: moorhen.Molecule, annotation: string): Promise<SessionOutput> => {
      if (!active) throw new Error("No session");
      const serialised = await serialiseMolecule(mol);
      if (!serialised) throw new Error("Moorhen returned no coordinates for this molecule");
      const formData = new FormData();
      formData.append("kind", "model");
      formData.append("annotation", annotation || mol.name);
      formData.append("file", serialised.blob, serialised.filename);
      const dropped = unwrap<SessionOutput>(
        await apiUpload(`jobs/${jobId}/interactive_drop/`, formData),
      );
      await refresh();
      return dropped;
    },
    [active, jobId, refresh],
  );

  const finish = useCallback(async (): Promise<FinishDisposition> => {
    if (!active) throw new Error("No session");
    const result = unwrap<SessionState & { disposition: FinishDisposition }>(
      await apiPost(`jobs/${jobId}/interactive_finish/`, { finished: true }),
    );
    finishedRef.current = true;
    setState(result);
    return result.disposition;
  }, [active, jobId]);

  // One stable object per change of its parts: the wrapper hands it to
  // Moorhen's side panels through useMemo, and a fresh object every render
  // would rebuild those panels on every render.
  return useMemo(
    () =>
      !active || !jobId
        ? null
        : { jobId, job, state, error, isOpen, refresh, saveMolecule, finish },
    [active, jobId, job, state, error, isOpen, refresh, saveMolecule, finish],
  );
}
