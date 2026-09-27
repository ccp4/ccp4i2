/**
 * A job status as a word, keyed on the enum and never on literals.
 *
 * The literal version in the task modal said 4=Failed, 5=Unsatisfactory
 * where Job.Status has 4=Interrupted, 5=Failed, 10=Unsatisfactory, so a
 * failed job read "Unsatisfactory" and an interrupted one "Failed"
 * (job-view.tsx and job-view-tabs.ts had already fixed the same off-by-one;
 * the modal was missed by that sweep). One function, one place.
 */
import { JobStatus } from "../types/models";

export function getStatusLabel(status: number): string {
  switch (status) {
    case JobStatus.UNKNOWN: return "Unknown";
    case JobStatus.PENDING: return "Pending";
    case JobStatus.QUEUED: return "Queued";
    case JobStatus.RUNNING: return "Running";
    case JobStatus.INTERRUPTED: return "Interrupted";
    case JobStatus.FAILED: return "Failed";
    case JobStatus.FINISHED: return "Finished";
    case JobStatus.RUNNING_REMOTELY: return "Running remotely";
    case JobStatus.FILE_HOLDER: return "File holder";
    case JobStatus.TO_DELETE: return "To delete";
    case JobStatus.UNSATISFACTORY: return "Unsatisfactory";
    default: return `Unknown (${status})`;
  }
}
