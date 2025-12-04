import {
  File as FileInfo,
  Job as JobInfo,
  Project as ProjectInfo,
} from "../../types/models";
import { moorhen } from "moorhen/types/moorhen";
import { apiGet } from "../../api-fetch";

export interface ItemMetadata {
  fileId: number;
  projectName?: string;
  jobNumber?: string;
  fileAnnotation?: string;
  isLoading: boolean;
  error?: string;
}

export function extractFileId(uniqueId: string): number | null {
  const match = uniqueId.match(/\/api\/proxy\/files\/(\d+)\/download\//);
  return match ? parseInt(match[1], 10) : null;
}

export async function fetchItemMetadata(
  item: moorhen.Molecule | moorhen.Map
): Promise<ItemMetadata | null> {
  const fileId = extractFileId(item.uniqueId || "");
  if (!fileId) return null;

  try {
    // Step 1: Fetch file information
    const fileInfo: FileInfo = await apiGet(`files/${fileId}`);

    // Step 2: Fetch job information
    const jobInfo: JobInfo = await apiGet(`jobs/${fileInfo.job}`);

    // Step 3: Fetch project information
    const projectInfo: ProjectInfo = await apiGet(`projects/${jobInfo.project}`);

    return {
      fileId,
      projectName: projectInfo.name,
      jobNumber: jobInfo.number,
      fileAnnotation: fileInfo.annotation || fileInfo.job_param_name,
      isLoading: false,
    };
  } catch (error) {
    return {
      fileId,
      isLoading: false,
      error: error instanceof Error ? error.message : "Unknown error",
    };
  }
}
