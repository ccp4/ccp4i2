// @vitest-environment node
import { describe, expect, it } from "vitest";
import { estimateRowHeight, ROW_BASE, ROW_WITH_SMILES } from "../lib/campaign-row-height";

const jobs = (n: number, sub = 0) => [
  ...Array.from({ length: n }, (_, i) => ({ number: String(i + 1) })),
  ...Array.from({ length: sub }, (_, i) => ({ number: `1.${i + 1}` })),
];
const opts = { jobsCollapsed: false, showSubJobs: false, hasSmiles: false };

describe("estimateRowHeight", () => {
  it("is one line for a few jobs", () => {
    expect(estimateRowHeight({ name: "x", jobs: jobs(3) }, opts)).toBe(ROW_BASE + 6);
  });

  it("grows with the lines fifteen job icons wrap onto", () => {
    expect(estimateRowHeight({ name: "x", jobs: jobs(15) }, opts)).toBe(12 + 3 * 50);
  });

  it("is one line when the jobs are folded, however many there are", () => {
    expect(estimateRowHeight({ name: "x", jobs: jobs(15) }, { ...opts, jobsCollapsed: true })).toBe(ROW_BASE);
  });

  it("counts subjobs only when they are shown", () => {
    const row = { name: "x", jobs: jobs(5, 10) };
    expect(estimateRowHeight(row, opts)).toBe(12 + 50);
    expect(estimateRowHeight(row, { ...opts, showSubJobs: true })).toBe(12 + 3 * 50);
  });

  it("leaves room for a ligand drawing", () => {
    expect(estimateRowHeight({ name: "x", jobs: jobs(1) }, { ...opts, hasSmiles: true })).toBe(ROW_WITH_SMILES);
  });
});
