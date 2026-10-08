/**
 * The evidence for one PanDDA event, as the campaign viewer loads it.
 *
 * A filled box in the campaign's site matrix is an event PanDDA found at that
 * site in that dataset. Clicking it opens the dataset's current model at the
 * site; what makes the event judgeable is its own evidence on top of that:
 * the event map (contoured where the receipt says, in absolute map units --
 * a default sigma contour is wrong for a BDC-corrected map) and the autobuilt
 * pose, drawn with the dictionary PanDDA built it from.
 *
 * The box names the event as `<receipt job id>:<event idx>` in the page URL's
 * `event` parameter. The receipt (a `pandda_events` job in the member
 * project) answers for the rest through its `eventEvidence` plugin method,
 * reached through the generic object_method endpoint.
 *
 * Everything but `fetchEventEvidence` is pure, to be tested without a viewer.
 */

import { apiPost } from "../api-fetch";

/** One of the receipt's `File` rows, as the method returns it. */
export interface EvidenceFile {
  id: number;
  uuid: string;
  name: string;
  type: string | null;
  annotation: string;
}

/** What `pandda_events.eventEvidence` returns in `data`. */
export interface EventEvidence {
  receipt_job_id: number;
  dtag: string;
  event_idx: number;
  position: number;
  site_idx: number | null;
  centroid: [number, number, number] | null;
  ligand_id: string | null;
  display_contour: number | null;
  optimal_contour: number | null;
  /** Where to contour the event map: display, else optimal, else null. */
  contour: number | null;
  event_map: EvidenceFile | null;
  pose: EvidenceFile | null;
  dictionary: EvidenceFile | null;
  has_map: boolean;
  has_pose: boolean;
  /** The event map's colour, as the receipt's own scenes draw it. */
  colour: string;
  radius: number;
}

/** An event, as the `event` URL parameter names it. */
export interface EventRef {
  receiptJobId: number;
  eventIdx: number;
}

/** `<receipt job id>:<event idx>`, or null when either is missing. */
export function eventParam(
  receiptJobId: number | null | undefined,
  eventIdx: number | null | undefined
): string | null {
  if (receiptJobId == null || eventIdx == null) return null;
  if (!Number.isInteger(receiptJobId) || !Number.isInteger(eventIdx)) return null;
  return `${receiptJobId}:${eventIdx}`;
}

/** The inverse of `eventParam`; null for anything malformed. */
export function parseEventParam(value: string | null | undefined): EventRef | null {
  if (!value) return null;
  const m = /^(\d+):(-?\d+)$/.exec(value.trim());
  if (!m) return null;
  return { receiptJobId: parseInt(m[1], 10), eventIdx: parseInt(m[2], 10) };
}

/** Where the viewer downloads a project file from, by its database id. */
export function fileDownloadUrl(fileId: number): string {
  return `/api/proxy/ccp4i2/files/${fileId}/download/`;
}

/** "#rrggbb" (or "#rrggbbaa") as the {r, g, b} 0-1 triple Moorhen wants. */
export function hexToRgb01(hex: string): { r: number; g: number; b: number } | null {
  const m = /^#([0-9a-fA-F]{2})([0-9a-fA-F]{2})([0-9a-fA-F]{2})(?:[0-9a-fA-F]{2})?$/.exec(hex);
  if (!m) return null;
  return {
    r: parseInt(m[1], 16) / 255,
    g: parseInt(m[2], 16) / 255,
    b: parseInt(m[3], 16) / 255,
  };
}

/** What the viewer should load for an event, and under what names. */
export interface EvidenceLoadPlan {
  map: {
    fileId: number;
    name: string;
    /** Absolute map units; null leaves the viewer's own default. */
    contourLevel: number | null;
    colour: string;
    radius: number;
  } | null;
  pose: {
    fileId: number;
    name: string;
    /** The dictionary to attach to the pose molecule alone, if any. */
    dictionaryFileId: number | null;
  } | null;
}

export function evidenceLoadPlan(evidence: EventEvidence): EvidenceLoadPlan {
  const n = evidence.event_idx;
  return {
    map: evidence.event_map
      ? {
          fileId: evidence.event_map.id,
          name: `Event ${n} map`,
          contourLevel: evidence.contour ?? null,
          colour: evidence.colour,
          radius: evidence.radius,
        }
      : null,
    pose: evidence.pose
      ? {
          fileId: evidence.pose.id,
          name: `Event ${n} autobuild`,
          dictionaryFileId: evidence.dictionary?.id ?? null,
        }
      : null,
  };
}

/**
 * Ask the receipt for one event's evidence. Null when it has no such event
 * or cannot be asked; the caller then loads nothing more and says nothing
 * alarming, since the dataset's own model is already in view.
 */
export async function fetchEventEvidence(ref: EventRef): Promise<EventEvidence | null> {
  try {
    const response: any = await apiPost(`jobs/${ref.receiptJobId}/object_method/`, {
      object_path: "pandda_events",
      method_name: "eventEvidence",
      args: [ref.eventIdx],
      kwargs: {},
    });
    const result = response?.data?.result;
    if (!result || result.success === false) {
      console.warn("[event evidence]", result?.error ?? "no result");
      return null;
    }
    return (result.data ?? null) as EventEvidence | null;
  } catch (err) {
    console.warn("[event evidence] request failed:", err);
    return null;
  }
}
