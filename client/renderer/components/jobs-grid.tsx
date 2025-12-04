import { useApi } from "../api";
import { Grid2, GridSize } from "@mui/material";
import { Job } from "../types/models";
import { JobCard } from "./job-card";

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
  const api = useApi();
  const { data: jobs } = api.get_endpoint<Job[]>({
    type: "projects",
    id: projectId,
    endpoint: "jobs",
  });
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
