/**
 * Where a scene's deployment-coupled file refs live, as URLs.
 *
 * A `job`+`param` ref names a job's output by role; a `projectPath` ref
 * names a file by its path in a project directory. Both need a project to
 * resolve against: the one the ref names (by uuid or name), else the page's.
 * Shared by every Moorhen viewer (components/moorhen/use-moorhen-loaders),
 * so the two pages agree on what a ref means.
 */
import { apiGet } from "../api-fetch";
import { projectPathUrl } from "./scene-project-path";

/** The project a page is scoped to, as far as it is known. */
export interface ProjectContext {
  pk: number | null;
  uuid: string | null;
  name: string | null;
}

export const EMPTY_PROJECT_CONTEXT: ProjectContext = { pk: null, uuid: null, name: null };

/** The mime type a job's authored Moorhen scene registers under. */
export const SCENE_FILE_TYPE = "application/moorhen-scene";

/** Map a project UUID to its primary key (the form the REST filters need). The
 *  scene carries the uuid; if it matches the current page's project we already
 *  have the pk, otherwise look it up in the projects list. Returns null (and
 *  logs) if it can't be resolved. */
export async function resolveProjectPk(
  want: { uuid?: string; name?: string },
  current: ProjectContext,
): Promise<number | null> {
  if (!want.uuid && !want.name) return current.pk; // no scene project ref → page context
  // The current window's project is authoritative when the scene names it (by
  // uuid OR name) — use its pk directly, no lookup.
  if (current.pk != null) {
    if (want.uuid && want.uuid === current.uuid) return current.pk;
    if (want.name && want.name === current.name) return current.pk;
  }
  // Otherwise resolve the *named* project exactly via the server-side uuid/name
  // filter (no list paging, no same-name ambiguity from the full list).
  const query = want.uuid
    ? `projects/?uuid=${encodeURIComponent(want.uuid)}`
    : `projects/?name=${encodeURIComponent(want.name as string)}`;
  try {
    const resp = await apiGet(query);
    const list = Array.isArray(resp) ? resp : (resp?.results ?? []);
    if (list.length === 1 && list[0]?.id != null) return list[0].id;
    if (list.length > 1) {
      console.warn(`[scene] project+job+param: "${want.uuid ?? want.name}" matched ${list.length} projects — cannot disambiguate`);
    } else {
      console.warn(`[scene] project+job+param: project "${want.uuid ?? want.name}" not found`);
    }
  } catch (e) {
    console.warn(`[scene] project+job+param: project lookup failed for "${want.uuid ?? want.name}"`, e);
  }
  // A project was explicitly named but couldn't be resolved. Do NOT fall back to
  // the window's project — that silently loads files from the wrong project.
  return null;
}

/**
 * Resolve a scene `job` (number, e.g. "1" or "1.1") + output `param` (e.g.
 * "XYZOUT") to a project file's proxy download URL, using the existing REST:
 * resolve the project pk (preferring the scene's projectId) → jobs in that
 * project → match number → that job's output files (directory=1) → match
 * job_param_name. Returns null (and logs the reason) if anything is missing.
 */
export async function resolveJobParamUrl(
  ref: { job: number | string; param: string; projectId?: string; projectName?: string },
  current: ProjectContext,
): Promise<string | null> {
  const projectPk = await resolveProjectPk({ uuid: ref.projectId, name: ref.projectName }, current);
  if (projectPk == null) {
    console.warn(`[scene] project+job+param ${ref.job}/${ref.param}: no project context (projectId=${ref.projectId ?? "none"}, projectName=${ref.projectName ?? "none"}, page pk null)`);
    return null;
  }
  const jobsResp = await apiGet(`jobs/?project=${projectPk}`);
  const jobs = Array.isArray(jobsResp) ? jobsResp : (jobsResp?.results ?? []);
  if (jobs.length === 0) {
    console.warn(`[scene] project+job+param: no jobs returned for project pk ${projectPk} (resp shape: ${Array.isArray(jobsResp) ? "array" : typeof jobsResp})`);
    return null;
  }
  const job = jobs.find((j: { number?: string }) => String(j.number) === String(ref.job));
  if (!job) {
    const nums = jobs.map((j: { number?: string }) => j.number).join(", ");
    console.warn(`[scene] project+job+param: job number "${ref.job}" not found in project pk ${projectPk}. Job numbers present: ${nums}`);
    return null;
  }
  const filesResp = await apiGet(`files/?job=${job.id}`);
  const files = Array.isArray(filesResp) ? filesResp : (filesResp?.results ?? []);
  if (files.length === 0) {
    console.warn(`[scene] project+job+param: no files returned for job ${ref.job} (pk ${job.id})`);
    return null;
  }
  // Match by output param name only — do NOT filter on directory. Job outputs
  // live in the job dir (directory 1), but imported coordinates live in the
  // import dir (directory 2); both are legitimate job+param targets.
  const file = files.find(
    (f: { job_param_name?: string }) => f.job_param_name === ref.param,
  );
  if (!file) {
    const avail = (files as { job_param_name?: string; directory?: number; type?: string }[])
      .filter((f) => f.job_param_name)
      .map((f) => `${f.job_param_name} (dir ${f.directory}, ${f.type})`)
      .join("; ");
    console.warn(`[scene] project+job+param: param "${ref.param}" not found in job ${ref.job}. Available: ${avail || "none"}`);
    return null;
  }
  return `/api/proxy/ccp4i2/files/${file.id}/download/`;
}

/** The URL of a scene's projectPath ref: the scene's project (by uuid or
 *  name) when it names one, else this page's. */
export async function resolveProjectPathUrl(
  ref: { projectPath?: string; projectId?: string; projectName?: string },
  current: ProjectContext,
): Promise<string | null> {
  if (!ref.projectPath) return null;
  const projectPk = await resolveProjectPk({ uuid: ref.projectId, name: ref.projectName }, current);
  if (projectPk == null) {
    console.warn(`[scene] projectPath ${ref.projectPath}: no project context`);
    return null;
  }
  return projectPathUrl(projectPk, ref.projectPath);
}

/**
 * The URL a scene file ref's bytes come from, for every ref kind that has
 * one: pdb (PDBe proxy), fileId (download by id), job+param, projectPath,
 * url. Bundle and cifText refs have none (their bytes are in hand).
 */
export async function sceneRefUrl(
  ref: {
    pdb?: string;
    fileId?: number;
    job?: number | string;
    param?: string;
    projectPath?: string;
    projectId?: string;
    projectName?: string;
    url?: string;
  },
  current: ProjectContext,
): Promise<string | null> {
  if (ref.pdb) return `/api/proxy/pdbe/entry-files/download/${ref.pdb.toLowerCase()}.cif`;
  // A ccp4i2 fileId is globally unique and the download URL keys on it alone:
  // no project qualifier needed. projectId/projectName are advisory context.
  if (ref.fileId !== undefined) return `/api/proxy/ccp4i2/files/${ref.fileId}/download/`;
  if (ref.job !== undefined && ref.param) {
    return resolveJobParamUrl(
      { job: ref.job, param: ref.param, projectId: ref.projectId, projectName: ref.projectName },
      current,
    );
  }
  if (ref.projectPath) return resolveProjectPathUrl(ref, current);
  if (ref.url) return ref.url;
  return null;
}
