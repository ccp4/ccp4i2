/**
 * The default projects directory setting, as the server reports and accepts
 * it. The Preferences panel and the New Project page both go through here, so
 * only one place knows the endpoint and the shape of its answer.
 */
import { apiGet, apiPatch } from "../api-fetch";

const ENDPOINT = "config/default-project-parent/";

export interface DefaultProjectsDir {
  /** Where a new project goes unless the user picks somewhere else. */
  directory: string;
  /** What a reset restores. */
  default: string;
  /** False in a deployment, which sets this with CCP4I2_PROJECTS_DIR. */
  editable: boolean;
}

const UNKNOWN: DefaultProjectsDir = {
  directory: "",
  default: "",
  editable: false,
};

export async function getDefaultProjectsDir(): Promise<DefaultProjectsDir> {
  const resp = await apiGet<any>(ENDPOINT);
  return { ...UNKNOWN, ...(resp?.data ?? resp) };
}

/** Store a new default, or reset to the built-in one by passing null. */
export async function setDefaultProjectsDir(
  directory: string | null
): Promise<void> {
  await apiPatch(ENDPOINT, { directory });
}
