"use client";
import { use } from "react";
import { Container, LinearProgress } from "@mui/material";
import { JobsGrid } from "../../../../components/jobs-grid";
import { useProject } from "../../../../utils";

export default function JobsPage({
  params,
}: {
  params: Promise<{ id: string }>;
}) {
  const { id } = use(params);
  const { project } = useProject(parseInt(id));
  if (!project) return <LinearProgress />;
  return (
    <Container>
      <JobsGrid projectId={parseInt(id)} size={{ xs: 12, md: 4, lg: 2 }} />
    </Container>
  );
}
