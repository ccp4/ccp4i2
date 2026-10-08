/**
 * A guess at a campaign table row's height, before it is drawn and measured.
 *
 * The table is virtualised: rows off screen are not drawn, and the space they
 * will take is reserved from this guess until each is measured. A flat guess
 * (85 px) left gaps, or overlaps, that jumped as rows came into view: a row
 * with its jobs folded is ~56 px, one whose fifteen job icons wrap is three
 * lines tall, and one with a 2D drawing of its ligand is taller still.
 */
import { parseDatasetFilename } from "../types/campaigns";

/** A row with one line of everything. */
export const ROW_BASE = 56;
/** A row with a ligand drawing (75 px plus its padding). */
export const ROW_WITH_SMILES = 95;
/** One line of job icons, with their numbers under them. */
export const JOB_LINE = 50;
/** Job icons on one line at a typical width; a guess, measured afterwards. */
export const JOBS_PER_LINE = 5;

interface RowLike {
  name: string;
  jobs?: { number: string }[];
}

export function hasSmiles(project: RowLike | undefined, smilesMap: Record<number, string>): boolean {
  if (!project) return false;
  const parsed = parseDatasetFilename(project.name);
  const regId = parsed.nclId ? parseInt(parsed.nclId) : null;
  return Boolean(regId && smilesMap[regId]);
}

export function estimateRowHeight(
  project: RowLike | undefined,
  options: { jobsCollapsed: boolean; showSubJobs: boolean; hasSmiles: boolean }
): number {
  const base = options.hasSmiles ? ROW_WITH_SMILES : ROW_BASE;
  if (!project || options.jobsCollapsed) return base;
  const jobs = (project.jobs ?? []).filter(
    (j) => options.showSubJobs || !j.number.includes(".")
  ).length;
  const lines = Math.max(1, Math.ceil(jobs / JOBS_PER_LINE));
  return Math.max(base, 12 + lines * JOB_LINE);
}
