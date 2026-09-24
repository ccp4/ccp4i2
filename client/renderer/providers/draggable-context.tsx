import {
  DndContext,
  DragOverlay,
  KeyboardSensor,
  PointerSensor,
  useSensor,
  useSensors,
} from "@dnd-kit/core";
import { PropsWithChildren, useCallback } from "react";
import { useCCP4i2Window } from "../app-context";
import { File, Job, Project } from "../types/models";
import { useJob, useProjectJobs } from "../utils";
import { useApi } from "../api";
import { useTaskInterface } from "./task-provider";
import { JobDragChip } from "../components/job-drag-chip";

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

  // Rows in the job list are both draggable and clickable. Without a distance
  // threshold the pointer sensor claims the pointerdown and the click that
  // would open the job never arrives.
  const sensors = useSensors(
    useSensor(PointerSensor, { activationConstraint: { distance: 8 } }),
    useSensor(KeyboardSensor)
  );

  const isValidDrop = (file: File, item: any) => {
    if (!file) return false;
    if (!item) return false;
    return file.type === item._qualifiers?.mimeTypeName;
  };

  const handleDragEnd = async (event: any) => {
    setActiveDragItem(null);
    // A click, or a drag released away from any droppable, ends with no `over`.
    const over = event.over?.data?.current;
    if (!over) return;
    if (event.active.data?.current?.job) {
      const context_job = event.active.data.current.job as Job;
      if (!over.job) return;
      if (!over.item) {
        setContextJob(over.job as Job, context_job);
      }
    } else if (event.active.data?.current?.file) {
      const file = event.active.data.current.file as File;
      if (over.job?.status !== 1) return;
      if (over.item) {
        if (!isValidDrop(file, over.item)) return;
        setFileByDrop(over.job as Job, over.item._objectPath, file);
      }
    }
  };

  // The dnd-kit draggable stores { job, file } as its data; files are
  // disabled above, so a live overlay is always a job.
  const dragJob = (activeDragItem as { job?: Job } | null)?.job;

  return (
    <DndContext
      sensors={sensors}
      onDragEnd={handleDragEnd}
      onDragCancel={() => setActiveDragItem(null)}
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
        {/* dnd-kit is disabled for files (they use native HTML5 drag), so a
            live overlay only ever previews a job. */}
        {dragJob ? <JobDragChip job={dragJob} /> : null}
      </DragOverlay>
    </DndContext>
  );
};
