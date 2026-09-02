// @vitest-environment node
import { describe, it, expect } from "vitest";
import { JobStatus } from "../types/models";
import { TAB, landingTab, visibleTabs } from "../components/job-view-tabs";

/**
 * Opening a job puts you on a tab. Which tab was decided by one expression of
 * bare numbers and which tabs exist by another, and the two disagreed: the
 * comment above the first read "4=Failed, 5=Unsatisfactory" where the enum has
 * 4=Interrupted, 5=Failed, 10=Unsatisfactory. So Interrupted and Unsatisfactory
 * jobs were sent to a Diagnostics tab shown only for Failed, and the clamp
 * bounced them to the task interface — a job that had run, opening on its
 * parameters, for no reason a user could see.
 *
 * The invariant is one line, and it is the one worth holding.
 */
const EVERY_STATUS = Object.values(JobStatus).filter(
  (v): v is JobStatus => typeof v === "number"
);

describe("landingTab", () => {
  it.each(EVERY_STATUS)("lands status %i on a tab that exists", (status) => {
    expect(visibleTabs(status, false).has(landingTab(status))).toBe(true);
  });

  it("sends a pending job to its parameters", () => {
    expect(landingTab(JobStatus.PENDING)).toBe(TAB.TASK_INTERFACE);
  });

  it("sends a failed job to its own record of why", () => {
    expect(landingTab(JobStatus.FAILED)).toBe(TAB.DIAGNOSTICS);
  });

  it.each([JobStatus.FINISHED, JobStatus.INTERRUPTED, JobStatus.UNSATISFACTORY])(
    "sends status %i to the results",
    (status) => {
      expect(landingTab(status)).toBe(TAB.REPORT);
    }
  );
});

describe("visibleTabs", () => {
  it("offers a failed job its report, which now carries the failure and the files", () => {
    expect(visibleTabs(JobStatus.FAILED, false).has(TAB.REPORT)).toBe(true);
  });

  it("keeps Diagnostics for the status that has one", () => {
    expect(visibleTabs(JobStatus.FAILED, false).has(TAB.DIAGNOSTICS)).toBe(true);
    expect(visibleTabs(JobStatus.FINISHED, false).has(TAB.DIAGNOSTICS)).toBe(false);
  });

  it("shows no report for a job that has not run", () => {
    expect(visibleTabs(JobStatus.PENDING, false).has(TAB.REPORT)).toBe(false);
    expect(visibleTabs(JobStatus.QUEUED, false).has(TAB.REPORT)).toBe(false);
  });

  it("shows everything in dev mode, whatever the status", () => {
    const visible = visibleTabs(JobStatus.PENDING, true);
    Object.values(TAB).forEach((tab) => expect(visible.has(tab)).toBe(true));
  });
});
