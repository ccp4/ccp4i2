/**
 * A job status as a colour, keyed on the enum and in one place.
 *
 * There were four maps: hex colours in job-avatar.tsx and the campaigns
 * tables, MUI severity names in the hierarchy browser and the file browser.
 * None of them covered every status, and each stopped at a different point,
 * so a job could be one colour in the job list and another in a campaign.
 *
 * What that cost: RUNNING_REMOTELY (7) and UNSATISFACTORY (10) both fell
 * through to the default, which is the same grey as UNKNOWN (0) -- a job
 * running on an Azure Batch node and a job that finished badly were drawn
 * identically, and identically to a job whose state nobody knows. Neither
 * colour was ever chosen; both were the absence of one.
 *
 * Same treatment as job-status-label: one function, keyed on JobStatus, so
 * adding a status to the enum shows up here as a missing case rather than as
 * silence in the interface.
 */
import { JobStatus } from "../types/models";

/**
 * Background colours, in the palette job-avatar established. The five that
 * existed keep their values, so nothing familiar changes.
 */
export function getStatusColour(status: number): string {
  switch (status) {
    case JobStatus.PENDING:
      return "#FFF"; // not started: empty
    case JobStatus.QUEUED:
      return "#FFA"; // waiting its turn: amber
    case JobStatus.RUNNING:
      return "#AAF"; // working here: blue
    case JobStatus.RUNNING_REMOTELY:
      return "#AEF"; // working elsewhere: the same family, one step away
    case JobStatus.INTERRUPTED:
      return "#FDA"; // stopped on purpose: orange
    case JobStatus.FAILED:
      return "#FAA"; // red
    case JobStatus.FINISHED:
      return "#AFA"; // green
    case JobStatus.UNSATISFACTORY:
      return "#DDA"; // it finished, and it is not good: olive, not green
    case JobStatus.FILE_HOLDER:
      return "#EEE"; // not a run at all
    case JobStatus.TO_DELETE:
      return "#CCC";
    case JobStatus.UNKNOWN:
    default:
      return "#AAA";
  }
}

/** The MUI colour name, for the places that draw chips rather than avatars. */
export function getStatusSeverity(
  status: number
): "default" | "info" | "success" | "warning" | "error" {
  switch (status) {
    case JobStatus.QUEUED:
    case JobStatus.RUNNING:
    case JobStatus.RUNNING_REMOTELY:
      return "info";
    case JobStatus.FINISHED:
      return "success";
    case JobStatus.FAILED:
      return "error";
    case JobStatus.PENDING:
    case JobStatus.INTERRUPTED:
    case JobStatus.UNSATISFACTORY:
      return "warning";
    default:
      return "default";
  }
}

/**
 * Whether the job is doing something, so the interface can show that it is.
 *
 * RUNNING_REMOTELY belongs here: the job's program is running on a Batch
 * node, which is no less running for being somewhere else. The avatar's pulse
 * was gated on RUNNING alone, so a dispatched job sat looking inert for the
 * hours its run took.
 */
export function isWorking(status: number): boolean {
  return (
    status === JobStatus.RUNNING ||
    status === JobStatus.RUNNING_REMOTELY ||
    status === JobStatus.QUEUED
  );
}
