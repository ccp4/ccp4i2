"use client";
import { useCallback, useEffect, useRef, useState } from "react";

/**
 * Poll the backend /health endpoint until it answers, so the UI can wait for the
 * Django/uvicorn server to actually be up rather than racing it.
 *
 * The health endpoint (api/ccp4i2/health/) is public (no auth) and returns 200
 * when the DB connection is live, 503 otherwise. It is reachable in both Electron
 * and web via the Next proxy, so we don't need the raw uvicorn port.
 *
 * Phases:
 *   "checking" — a probe is in flight or scheduled
 *   "ready"    — the server answered healthy
 *   "timeout"  — gave up after `timeoutMs` without a healthy response
 */
export type ServerReadyStatus = "checking" | "ready" | "timeout";

const HEALTH_URL = "/api/proxy/ccp4i2/health/";

export function useServerReady(options?: {
  intervalMs?: number;
  timeoutMs?: number;
  enabled?: boolean;
}) {
  const intervalMs = options?.intervalMs ?? 1000;
  const timeoutMs = options?.timeoutMs ?? 60_000;
  const enabled = options?.enabled ?? true;

  const [status, setStatus] = useState<ServerReadyStatus>("checking");
  const [attempts, setAttempts] = useState(0);
  // Bumped by retry() to restart the polling loop.
  const [nonce, setNonce] = useState(0);

  const retry = useCallback(() => {
    setAttempts(0);
    setStatus("checking");
    setNonce((n) => n + 1);
  }, []);

  useEffect(() => {
    if (!enabled) return;

    let cancelled = false;
    let timer: ReturnType<typeof setTimeout> | null = null;
    const startedAt = Date.now(); // browser Date is fine here

    const probe = async (): Promise<boolean> => {
      try {
        const res = await fetch(HEALTH_URL, {
          method: "GET",
          cache: "no-store",
          headers: { Accept: "application/json" },
        });
        return res.ok;
      } catch {
        return false; // network error while the server is still coming up
      }
    };

    const tick = async () => {
      if (cancelled) return;
      const ok = await probe();
      if (cancelled) return;

      if (ok) {
        setStatus("ready");
        return;
      }
      setAttempts((n) => n + 1);

      if (Date.now() - startedAt >= timeoutMs) {
        setStatus("timeout");
        return;
      }
      timer = setTimeout(tick, intervalMs);
    };

    tick();

    return () => {
      cancelled = true;
      if (timer) clearTimeout(timer);
    };
  }, [enabled, intervalMs, timeoutMs, nonce]);

  return { status, attempts, retry };
}
