/**
 * Which tabs a job panel shows, and which one it opens on.
 *
 * Kept apart from the component so both halves can be checked against each
 * other. They used to be two expressions of bare numbers in different parts of
 * one file, and they disagreed: the comment above the landing-tab expression
 * read "4=Failed, 5=Unsatisfactory" where `JobStatus` has 4=Interrupted,
 * 5=Failed, 10=Unsatisfactory. Interrupted and Unsatisfactory jobs were
 * therefore sent to a Diagnostics tab that is shown only for Failed, and the
 * clamp bounced them to the task interface.
 */
import { JobStatus } from "../types/models";

/**
 * The job panel's tabs, by the value each `<Tab>` carries.
 *
 * Named because the visibility rules test job *statuses*, which are numbers in
 * the same range: `[3, 4, 5, 6, 7, 9, 10]` reads as a list of tabs and is a
 * list of statuses.
 */
export const TAB = {
  TASK_INTERFACE: 0,
  PARAMS_XML: 1,
  REPORT_XML: 2,
  REPORT: 3,
  DIAGNOSTICS: 4,
  DEF_XML: 5,
  VALIDATION: 6,
  JOB_CONTAINER: 7,
  COMMENTS: 8,
  DIRECTORY: 9,
  LOGS: 10,
} as const;

const DEV_ONLY_TABS = [
  TAB.PARAMS_XML,
  TAB.REPORT_XML,
  TAB.DEF_XML,
  TAB.JOB_CONTAINER,
];

/**
 * Statuses for which there is a report worth offering.
 *
 * Failed is included: the report now carries the rendering failure and the
 * files the job did glean before it stopped, so it is a status that needs the
 * tab rather than one that has no use for it.
 */
export const REPORT_STATUSES: JobStatus[] = [
  JobStatus.RUNNING,
  JobStatus.INTERRUPTED,
  JobStatus.FAILED,
  JobStatus.FINISHED,
  JobStatus.RUNNING_REMOTELY,
  JobStatus.TO_DELETE,
  JobStatus.UNSATISFACTORY,
];

/** Which tabs exist for a job in this state. */
export const visibleTabs = (
  status: JobStatus | undefined,
  devMode: boolean
): Set<number> => {
  const visible = new Set<number>([
    TAB.TASK_INTERFACE,
    TAB.COMMENTS,
    TAB.DIRECTORY,
    TAB.LOGS,
  ]);
  if (devMode) DEV_ONLY_TABS.forEach((v) => visible.add(v));
  if (devMode || REPORT_STATUSES.includes(status as JobStatus))
    visible.add(TAB.REPORT);
  if (devMode || status === JobStatus.FAILED) visible.add(TAB.DIAGNOSTICS);
  if (devMode || status === JobStatus.PENDING) visible.add(TAB.VALIDATION);
  return visible;
};

/**
 * Where to land when a job is opened, given what it has to say.
 *
 * Only ever names a tab that `visibleTabs` agrees exists. The caller clamps as
 * well, but a function that returns a tab it knows is absent is relying on
 * someone else to be careful, and that reliance is what broke: Queued, Unknown
 * and File-holder jobs all asked for a Report tab they do not have.
 */
export const landingTab = (status: JobStatus): number => {
  if (status === JobStatus.PENDING) return TAB.TASK_INTERFACE;
  // The job's own record of why it stopped is the thing to read first.
  if (status === JobStatus.FAILED) return TAB.DIAGNOSTICS;
  // Anything that has run: the results. Anything that has not: its parameters.
  return REPORT_STATUSES.includes(status) ? TAB.REPORT : TAB.TASK_INTERFACE;
};
