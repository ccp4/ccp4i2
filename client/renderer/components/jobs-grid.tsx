import { Grid2, GridSize } from "@mui/material";
import { Job } from "../types/models";
import { JobCard } from "./job-card";
import { useProjectJobs } from "../utils";

interface SizeProps {
  xs?: GridSize | null;
  sm?: GridSize | null;
  md?: GridSize | null;
  lg?: GridSize | null;
  xl?: GridSize | null;
}
interface JobsGridProps {
  projectId: number;
  size: SizeProps;
  parent?: number;
  withSubtitles?: boolean;
}
export const JobsGrid: React.FC<JobsGridProps> = ({
  projectId,
  size = 4,
  parent = null,
  withSubtitles = false,
}) => {
  const { jobs } = useProjectJobs(projectId);
  return (
    <Grid2 container>
      {jobs &&
        jobs
          .reverse()
          .filter((item) => item.parent === parent)
          .map((job: Job) => (
            <Grid2 key={job.number} size={size}>
              <JobCard key={job.id} job={job} withSubtitle={withSubtitles} />
            </Grid2>
          ))}
    </Grid2>
  );
};
