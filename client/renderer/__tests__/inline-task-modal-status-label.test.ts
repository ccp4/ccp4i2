// @vitest-environment node
import { describe, it, expect } from "vitest";
import { JobStatus } from "../types/models";
import { getStatusLabel } from "../lib/job-status-label";

/**
 * The task modal's status label was one off the enum: 4 read "Failed" and 5
 * "Unsatisfactory", so a job that genuinely failed said "Unsatisfactory" and
 * an interrupted one said "Failed" (Materia, 2026-09-27: an operator on
 * Batch reasonably resubmitted a failed job to learn nothing). Assert against
 * the enum, never against literals, for every member.
 */
describe("getStatusLabel", () => {
  it("labels the three statuses the sweep fixed elsewhere against the enum", () => {
    expect(getStatusLabel(JobStatus.INTERRUPTED)).toBe("Interrupted");
    expect(getStatusLabel(JobStatus.FAILED)).toBe("Failed");
    expect(getStatusLabel(JobStatus.UNSATISFACTORY)).toBe("Unsatisfactory");
    expect(getStatusLabel(JobStatus.RUNNING_REMOTELY)).toBe("Running remotely");
  });

  it("has a label for every enum member and never falls through for one", () => {
    const members = Object.values(JobStatus).filter((v): v is JobStatus => typeof v === "number");
    for (const status of members) {
      expect(getStatusLabel(status)).not.toMatch(/^Unknown \(/);
    }
    expect(getStatusLabel(99)).toBe("Unknown (99)");
  });
});
