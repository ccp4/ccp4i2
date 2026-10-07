/**
 * The Judgement tab draws a job's verdict: the outcome in its tone, each
 * clause of the condition that decided it with a gauge, the numbers with what
 * they mean, and the next steps a user can turn into a job.
 */
import { describe, expect, it, vi, beforeEach } from "vitest";
import { fireEvent, render, screen, waitFor } from "@testing-library/react";
import { JudgementPanel, ThresholdGauge } from "@/components/judgement/judgement-panel";

const post = vi.fn();
const push = vi.fn();
const setMessage = vi.fn();
let verdictResponse: any;

vi.mock("next/navigation", () => ({ useRouter: () => ({ push }) }));
vi.mock("@/providers/popcorn-provider", () => ({ usePopcorn: () => ({ setMessage }) }));
vi.mock("@/api", () => ({
  useApi: () => ({
    get_endpoint: () => ({ data: verdictResponse, isLoading: false }),
    get: (url: string) => ({
      data: url.startsWith("agent/tasks/")
        ? { success: true, data: { judgement: { traps: ["Read the final table, not the first ones."] } } }
        : { modelcraft: { TASKTITLE: "Autobuild with ModelCraft" } },
    }),
    post,
  }),
}));

const job: any = { id: 7, number: "5", project: 3, task_name: "mrbump_basic", status: 6 };

beforeEach(() => {
  post.mockReset();
  push.mockReset();
  setMessage.mockReset();
  verdictResponse = {
    success: true,
    data: {
      task: "mrbump_basic",
      outcome: "placed",
      tone: "caution",
      because: "TFZ >= 8 and RFREE < 0.55",
      basis: "A correct placement of a distant model.",
      clauses: [
        { text: "TFZ >= 8", holds: true, values: { TFZ: 12.7 }, name: "TFZ", op: ">=", threshold: 8, value: 12.7 },
        { text: "RFREE < 0.55", holds: true, values: { RFREE: 0.52 }, name: "RFREE", op: "<", threshold: 0.55, value: 0.52 },
      ],
      results: { TFZ: 12.7, RFREE: 0.52, MODEL_2: null },
      meanings: { TFZ: "Phaser's refined TFZ.", RFREE: "R-free after refinement." },
      missing: [],
      next: [
        { when: 'outcome == "placed"', task: "modelcraft", inputs: { XYZIN: "[5].XYZOUT" }, advice: "Build it." },
        { when: "true", task: "mrbump_basic", rerun: true, inputs: { NCYC: "20" }, advice: "More cycles." },
        { when: "true", advice: "Look at the map." },
      ],
      note: "A draft no crystallographer has reviewed.",
      judgement_version: "abc123",
    },
  };
});

describe("JudgementPanel", () => {
  it("shows the outcome, why, the numbers and the traps", () => {
    render(<JudgementPanel job={job} />);
    expect(screen.getByText("placed")).toBeTruthy();
    expect(screen.getByText("TFZ >= 8")).toBeTruthy();
    expect(screen.getByText("RFREE < 0.55")).toBeTruthy();
    expect(screen.getByText("Phaser's refined TFZ.")).toBeTruthy();
    expect(screen.getByText("absent")).toBeTruthy(); // MODEL_2, optional and not there
    expect(screen.getByText(/Read the final table/)).toBeTruthy();
    expect(screen.getByText(/judgement version abc123/)).toBeTruthy();
  });

  it("offers a job for a step with a task, a rerun for a rerun, and only advice otherwise", () => {
    render(<JudgementPanel job={job} />);
    expect(screen.getByText("Autobuild with ModelCraft")).toBeTruthy();
    expect(screen.getByRole("button", { name: /Create job/ })).toBeTruthy();
    expect(screen.getByRole("button", { name: /Rerun with these changes/ })).toBeTruthy();
    expect(screen.getAllByRole("button")).toHaveLength(2); // the advice-only step has none
    expect(screen.getByText("XYZIN = [5].XYZOUT")).toBeTruthy();
  });

  it("makes the job and opens it, saying what was set", async () => {
    post.mockResolvedValue({
      success: true,
      data: { job: { id: 99, number: "6" }, rerun: false, inputs: [{ name: "XYZIN", ok: true }] },
    });
    render(<JudgementPanel job={job} />);
    fireEvent.click(screen.getByRole("button", { name: /Create job/ }));
    await waitFor(() => expect(push).toHaveBeenCalledWith("/ccp4i2/project/3/job/99"));
    expect(post).toHaveBeenCalledWith("jobs/7/apply_next/", { index: 0 });
    expect(setMessage.mock.calls[0][1]).toBe("success");
  });

  it("says which inputs could not be set", async () => {
    post.mockResolvedValue({
      success: true,
      data: { job: { id: 99, number: "6" }, rerun: true, inputs: [{ name: "NCYC", ok: false, error: "x" }] },
    });
    render(<JudgementPanel job={job} />);
    fireEvent.click(screen.getByRole("button", { name: /Rerun/ }));
    await waitFor(() => expect(setMessage).toHaveBeenCalled());
    expect(setMessage.mock.calls[0]).toEqual(["Job 6 made; not set: NCYC", "warning"]);
  });
});

describe("ThresholdGauge", () => {
  it("draws a value against its threshold, and nothing without numbers", () => {
    const { container, rerender } = render(
      <ThresholdGauge clause={{ text: "R < 0.5", holds: true, values: {}, name: "R", op: "<", threshold: 0.5, value: 0.3 }} />
    );
    expect(container.querySelector("circle")).toBeTruthy();
    rerender(<ThresholdGauge clause={{ text: "X == null", holds: true, values: {} }} />);
    expect(container.querySelector("svg")).toBeNull();
  });
});
