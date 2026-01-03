import { DndContext, DragOverlay } from "@dnd-kit/core";
import { PropsWithChildren, useCallback } from "react";
import { useCCP4i2Window } from "../app-context";
import { Avatar } from "@mui/material";
import { File, Job, Project } from "../types/models";
import { useJob } from "../utils";
import { useApi } from "../api";
import { useTaskInterface } from "./task-provider";

export const DraggableContext: React.FC<PropsWithChildren> = (props) => {
  const { jobId } = useCCP4i2Window();

  const { setInFlight } = useTaskInterface();

  const api = useApi();

  const { job, mutateContainer, fileItemToParameterArg, setParameter } =
    useJob(jobId);

  const { activeDragItem, setActiveDragItem } = useCCP4i2Window();

  const { data: project_jobs } = api.get_endpoint<Job[]>({
    type: "projects",
    id: job?.project,
    endpoint: "jobs",
  });

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
    console.log(event);
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

  return (
    <DndContext
      onDragEnd={handleDragEnd}
      onDragStart={({ active }) => {
        console.log("Drag start", active);
        setActiveDragItem(active.data.current as File | Job);
      }}
    >
      {props.children}
      <DragOverlay>
        {activeDragItem ? (
          <Avatar
            src="/qticons/ccp4i2.png"
            sx={{
              width: 64,
              height: 64,
              boxShadow: 3,
              opacity: 0.8,
              pointerEvents: "none",
            }}
          />
        ) : null}
      </DragOverlay>
    </DndContext>
  );
};
