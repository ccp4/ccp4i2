import { useCallback, useEffect } from "react";
import { Job, Project } from "../types/models";
import { FileSystemItem } from "../components/directory-browser";
import { useFilePreviewContext } from "../providers/file-preview-context";
import { fileKind } from "../lib/file-kinds";

/**
 * Opens a file from a job's directory in the preview dialog. While the job
 * runs, a text file is shown live: the dialog reads it again every few
 * seconds, so a log can be watched as it grows (#611). When the job stops,
 * the open preview is read once more and then left alone.
 *
 * @param isJobActive whether the job is running (useJobDirectory gives it)
 */
export function useJobFilePreview(
  job: Job,
  project: Project,
  isJobActive: boolean
) {
  const { contentSpecification, setContentSpecification } =
    useFilePreviewContext();

  const openFile = useCallback(
    (item: FileSystemItem) => {
      if (item.type === "directory") return;
      const relative = item.path.slice(project.directory.length + 1);
      const language = fileKind(item.name).preview ?? "text";
      setContentSpecification({
        url: `/api/proxy/ccp4i2/projects/${project.id}/project_file?path=${encodeURIComponent(relative)}`,
        title: item.name,
        language,
        // Only plain text: a JSON or XML file read half-written does not parse.
        liveJobId: isJobActive && language === "text" ? job.id : undefined,
      });
    },
    [project.id, project.directory, job.id, isJobActive, setContentSpecification]
  );

  // The job has stopped: end the live preview. The new specification makes
  // the dialog read the file once more, so it shows the finished log.
  useEffect(() => {
    if (!isJobActive && contentSpecification?.liveJobId === job.id) {
      setContentSpecification({ ...contentSpecification, liveJobId: undefined });
    }
  }, [isJobActive, job.id, contentSpecification, setContentSpecification]);

  return openFile;
}
