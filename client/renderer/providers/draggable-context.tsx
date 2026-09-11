import { DndContext, DragOverlay } from "@dnd-kit/core";
import { PropsWithChildren, useCallback } from "react";
import { useCCP4i2Window } from "../app-context";
import { Avatar, Box, Typography } from "@mui/material";
import { File, Job, Project } from "../types/models";
import { useJob, useProjectJobs } from "../utils";
import { useApi } from "../api";
import { useTaskInterface } from "./task-provider";

export const DraggableContext: React.FC<PropsWithChildren> = (props) => {
  const { jobId } = useCCP4i2Window();

  const { setInFlight } = useTaskInterface();

  const api = useApi();

  const { job, mutateContainer, fileItemToParameterArg, setParameter } =
    useJob(jobId);

  const { activeDragItem, setActiveDragItem } = useCCP4i2Window();

  const { jobs: project_jobs } = useProjectJobs(job?.project);

  const { data: projects } = api.get<Project[]>("projects");

  const setContextJob = useCallback(
    (job: Job, context_job: Job) => {
      const formData = {
        context_job_uuid: context_job ? context_job.uuid : null,
      };
      api.post<Job>(`jobs/${job.id}/set_context_job`, formData).then(() => {
        if (mutateContainer) {
          mutateContainer();
        }
      });
    },
    [mutateContainer, api]
  );

  const setFileByDrop = useCallback(
    async (job: Job, objectPath: string, file: File) => {
      if (!job) return;
      if (!objectPath) return;
      if (!file) return;
      if (!project_jobs) return;
      if (!projects) return;
      const setParameterArg = fileItemToParameterArg(
        file,
        objectPath,
        project_jobs,
        projects
      );
      setInFlight(true);
      //Not sure how to trigger onChange
      await setParameter(setParameterArg);
      await mutateContainer();
      setInFlight(false);
    },
    [project_jobs, projects, fileItemToParameterArg]
  );

  const isValidDrop = (file: File, item: any) => {
    if (!file) return false;
    if (!item) return false;
    return file.type === item._qualifiers?.mimeTypeName;
  };

  const handleDragEnd = async (event: any) => {
    if (event.active.data?.current?.job) {
      const context_job = event.active.data.current.job as Job;
      if (!event.over.data?.current?.job) return;
      if (!event.over.data?.current?.item) {
        const job = event.over.data.current.job as Job;
        setContextJob(job, context_job);
      }
    } else if (event.active.data?.current?.file) {
      const file = event.active.data.current.file as File;
      if (
        !event.over.data?.current?.job ||
        event.over.data?.current?.job?.status !== 1
      )
        return;
      if (event.over.data?.current?.item) {
        if (!isValidDrop(file, event.over.data.current.item)) return;
        const job = event.over.data.current.job as Job;
        setFileByDrop(job, event.over.data?.current?.item._objectPath, file);
      }
    }
  };

  // The dnd-kit draggable stores { job, file } as its data; files are
  // disabled above, so a live overlay is always a job.
  const dragJob = (activeDragItem as { job?: Job } | null)?.job;

  return (
    <DndContext
      onDragEnd={handleDragEnd}
      onDragStart={({ active }) => {
        setActiveDragItem(active.data.current as File | Job);

        // Write reference to system clipboard for cross-window paste
        const data = active.data.current;
        if (data?.file) {
          const file = data.file as File;
          const ref = {
            ccp4i2_file: true,
            uuid: file.uuid,
            id: file.id,
            name: file.name,
            type: file.type,
            sub_type: file.sub_type,
            content: file.content,
            annotation: file.annotation,
            job: file.job,
            job_param_name: file.job_param_name,
          };
          navigator.clipboard.writeText(JSON.stringify(ref)).catch(() => {});
        } else if (data?.job) {
          const dragJob = data.job as Job;
          const ref = {
            ccp4i2_job: true,
            id: dragJob.id,
            uuid: dragJob.uuid,
            task_name: dragJob.task_name,
            title: dragJob.title,
            number: dragJob.number,
          };
          navigator.clipboard.writeText(JSON.stringify(ref)).catch(() => {});
        }
      }}
    >
      {props.children}
      <DragOverlay>
        {dragJob ? (
          // Show what is actually being dragged — the job as a compact,
          // row-like chip — rather than a generic offset icon. dnd-kit is
          // disabled for files (they use native HTML5 drag), so the overlay
          // only ever previews a job. Keeping it small and anchored reads as
          // the row itself moving, which is what the drag does (issue #432).
          <Box
            sx={{
              display: "inline-flex",
              alignItems: "center",
              gap: 1,
              px: 1.5,
              py: 0.5,
              maxWidth: 320,
              bgcolor: "background.paper",
              border: 1,
              borderColor: "divider",
              borderRadius: 1,
              boxShadow: 4,
              opacity: 0.95,
              pointerEvents: "none",
              cursor: "grabbing",
            }}
          >
            <Avatar
              src="/svgicons/ccp4i2.svg"
              sx={{ width: 20, height: 20 }}
            />
            <Typography
              variant="body2"
              sx={{
                fontWeight: 600,
                whiteSpace: "nowrap",
                overflow: "hidden",
                textOverflow: "ellipsis",
              }}
            >
              {dragJob.number ? `${dragJob.number}. ` : ""}
              {dragJob.title || dragJob.task_name || "Job"}
            </Typography>
          </Box>
        ) : null}
      </DragOverlay>
    </DndContext>
  );
};
