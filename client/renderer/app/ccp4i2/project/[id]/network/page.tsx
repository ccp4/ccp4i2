"use client";
import { useParams } from "next/navigation";
import { ProjectNetwork } from "../../../../../components/project-network";
import { useProject } from "../../../../../utils";
import ToolBar from "../../../../../components/tool-bar";
import Container from "@mui/material/Container/Container";

export default function Page() {
  const params = useParams();
  const { id } = params as { id: string };
  const { project } = useProject(parseInt(id));
  return (
    <div>
      <ToolBar />
      <Container>
        <h1>Project {project?.name} Network</h1>
        {id && <ProjectNetwork projectId={parseInt(id)} />}
      </Container>
    </div>
  );
}
