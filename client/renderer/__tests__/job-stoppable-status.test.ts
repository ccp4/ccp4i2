/**
 * "Stop job" in the job menu (#773) acts on queued and running jobs, local or
 * remote, and on nothing else: stopping a finished or pending job means
 * nothing, and the server's cancel would only mark it Interrupted.
 */
import { describe, expect, it } from "vitest";
import { JobStatus, isStoppableJobStatus } from "../types/models";

describe("isStoppableJobStatus", () => {
  it("is true for queued, running and running remotely", () => {
    for (const s of [JobStatus.QUEUED, JobStatus.RUNNING, JobStatus.RUNNING_REMOTELY]) {
      expect(isStoppableJobStatus(s)).toBe(true);
    }
  });

  it("is false for every other status", () => {
    for (const s of [
      JobStatus.UNKNOWN, JobStatus.PENDING, JobStatus.INTERRUPTED, JobStatus.FAILED,
      JobStatus.FINISHED, JobStatus.FILE_HOLDER, JobStatus.TO_DELETE, JobStatus.UNSATISFACTORY,
    ]) {
      expect(isStoppableJobStatus(s)).toBe(false);
    }
  });
});
