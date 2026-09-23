import { useCallback, useEffect, useRef, useState } from "react";

import { Job } from "../../../../types/models";
import { useJob } from "../../../../utils";
import { ChainInfo } from "./dm-spec";

/**
 * What the server can tell us about the user's own model and the assembly
 * they have described, without running anything.
 *
 * `dm_multidomain.ncs_preview()` is pure gemmi, so this costs a superposition
 * or two and no binary. The numbers that matter are per body per copy: how
 * many CA atoms matched, and what RMSD the superposition achieved. They are
 * the difference between "I typed some residue ranges" and "I can see whether
 * these residues really move as one unit" — the question the old interface
 * could only answer by running the job.
 */

export interface PreviewCopy {
  label: string;
  nCA?: number;
  rmsd?: number;
  skipped?: boolean;
  error?: string;
}

export interface PreviewBody {
  index: number;
  spec: string;
  mode: string;
  segments: Array<{ role: string; lo: number; hi: number }>;
  copies: PreviewCopy[];
  nReferenceCA?: number;
  error: string | null;
}

export interface PreviewMessage {
  code: number;
  severity: "error" | "warning";
  text: string;
  path?: string;
}

export interface NcsPreview {
  ok: boolean;
  error?: string;
  model?: {
    chains: ChainInfo[];
    entities: string[][];
    nCopiesDetected: number;
  };
  suggestion?: { assembly: string[]; segments: string };
  assembly?: {
    source: "parameter" | "detected";
    rows: string[];
    instances: Array<{ label: string; roles: Record<string, string> }>;
  };
  bodies?: PreviewBody[];
  messages?: PreviewMessage[];
}

/**
 * Fetch the preview, debounced on `key` — a digest of everything the preview
 * depends on (the model, the assembly rows, the body specs). Editing a residue
 * range should not fire a superposition per keystroke, but it should settle
 * quickly enough that the numbers feel attached to the field.
 */
export const useNcsPreview = (
  job: Job,
  key: string,
  enabled: boolean,
  delayMs = 400
): { preview: NcsPreview | null; isLoading: boolean; refresh: () => void } => {
  const { callPluginMethod } = useJob(job.id);
  const [preview, setPreview] = useState<NcsPreview | null>(null);
  const [isLoading, setIsLoading] = useState(false);
  const [nonce, setNonce] = useState(0);
  // Only the newest request may write to state: a slow superposition must not
  // land on top of the answer for what the user has since typed.
  const latest = useRef(0);

  const refresh = useCallback(() => setNonce((n) => n + 1), []);

  useEffect(() => {
    if (!enabled) {
      setPreview(null);
      return undefined;
    }
    const mine = latest.current + 1;
    latest.current = mine;
    setIsLoading(true);
    const timer = setTimeout(async () => {
      const result = await callPluginMethod("ncs_preview");
      if (latest.current !== mine) return;
      setPreview(
        result && typeof result === "object"
          ? (result as NcsPreview)
          : { ok: false, error: "The server did not answer" }
      );
      setIsLoading(false);
    }, delayMs);
    return () => clearTimeout(timer);
  }, [key, enabled, nonce, delayMs, callPluginMethod]);

  return { preview, isLoading, refresh };
};

/** Colour for an RMSD, read as "does this body hold together". */
export const rmsdTone = (rmsd: number): "success" | "warning" | "error" => {
  if (rmsd <= 1.5) return "success";
  if (rmsd <= 3.0) return "warning";
  return "error";
};
