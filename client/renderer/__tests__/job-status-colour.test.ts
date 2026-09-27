/**
 * The colours that were never chosen.
 *
 * Four components carried their own status-to-colour map and none covered
 * every status, so RUNNING_REMOTELY and UNSATISFACTORY both fell through to
 * the default -- the same grey as UNKNOWN. A job running on a Batch node, a
 * job that finished badly, and a job in an unknown state were drawn
 * identically.
 */
import { describe, expect, it } from "vitest";

import { getStatusColour, getStatusSeverity, isWorking } from "../lib/job-status-colour";
import { JobStatus } from "../types/models";

describe("job status colours", () => {
  it("gives a remotely running job its own colour, not Unknown's grey", () => {
    expect(getStatusColour(JobStatus.RUNNING_REMOTELY)).not.toBe(
      getStatusColour(JobStatus.UNKNOWN)
    );
  });

  it("gives an unsatisfactory job its own colour, not Unknown's grey", () => {
    expect(getStatusColour(JobStatus.UNSATISFACTORY)).not.toBe(
      getStatusColour(JobStatus.UNKNOWN)
    );
  });

  it("does not draw unsatisfactory as finished", () => {
    expect(getStatusColour(JobStatus.UNSATISFACTORY)).not.toBe(
      getStatusColour(JobStatus.FINISHED)
    );
  });

  it("keeps the colours the interface already had", () => {
    expect(getStatusColour(JobStatus.RUNNING)).toBe("#AAF");
    expect(getStatusColour(JobStatus.FAILED)).toBe("#FAA");
    expect(getStatusColour(JobStatus.FINISHED)).toBe("#AFA");
    expect(getStatusColour(JobStatus.QUEUED)).toBe("#FFA");
    expect(getStatusColour(JobStatus.INTERRUPTED)).toBe("#FDA");
  });

  it("gives every status a colour of its own", () => {
    const statuses = [
      JobStatus.UNKNOWN, JobStatus.PENDING, JobStatus.QUEUED, JobStatus.RUNNING,
      JobStatus.INTERRUPTED, JobStatus.FAILED, JobStatus.FINISHED,
      JobStatus.RUNNING_REMOTELY, JobStatus.FILE_HOLDER, JobStatus.TO_DELETE,
      JobStatus.UNSATISFACTORY,
    ];
    const colours = statuses.map(getStatusColour);
    expect(new Set(colours).size).toBe(statuses.length);
  });

  it("treats a remotely running job as working, so the interface can say so", () => {
    expect(isWorking(JobStatus.RUNNING_REMOTELY)).toBe(true);
    expect(isWorking(JobStatus.FINISHED)).toBe(false);
  });

  it("maps severities for the chip-based views too", () => {
    expect(getStatusSeverity(JobStatus.RUNNING_REMOTELY)).toBe("info");
    expect(getStatusSeverity(JobStatus.UNSATISFACTORY)).toBe("warning");
    expect(getStatusSeverity(JobStatus.FAILED)).toBe("error");
    expect(getStatusSeverity(JobStatus.FINISHED)).toBe("success");
  });
});
