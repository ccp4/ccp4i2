"use client";
import { Badge, Container, Skeleton, Stack } from "@mui/material";
import ProjectsToolbar from "../../components/projects-toolbar";
import ProjectsTable from "../../components/projects-table";
import CCP4i2TopBar from "../../components/ccp4i2-topbar";
import { useApi } from "../../api";
import { Project } from "../../types/models";

export default function ProjectsPage() {
  const api = useApi();
  const { data: projects } = api.get<Project[]>("projects", 1000);
  const { data: task_tree } = api.get<any>(`task_tree/`);

  return (
    <>
      <CCP4i2TopBar title="CCP4i2 Projects" />
      <Container sx={{ my: 3 }}>
      <Stack spacing={2}>
        {task_tree?.task_tree?.tree && (
          <Badge badgeContent={task_tree.task_tree.tree.length} color="primary">
            {" "}
            Tasks available
          </Badge>
        )}
        <ProjectsToolbar />
        {projects ? (
          projects.length > 0 ? (
            <ProjectsTable />
          ) : (
            <Stack alignItems="center" spacing={2} sx={{ py: 8 }}>
              <Badge color="secondary" badgeContent={0}>
                <span style={{ fontSize: 18, color: "#888" }}>
                  No projects found
                </span>
              </Badge>
              <span style={{ color: "#666" }}>
                <svg width="64" height="64" viewBox="0 0 24 24" fill="none">
                  <rect
                    x="3"
                    y="7"
                    width="18"
                    height="13"
                    rx="2"
                    fill="#e0e0e0"
                  />
                  <path
                    d="M3 7V5a2 2 0 0 1 2-2h3.17a2 2 0 0 1 1.41.59l1.83 1.83H19a2 2 0 0 1 2 2v1"
                    stroke="#bdbdbd"
                    strokeWidth="2"
                    strokeLinecap="round"
                    strokeLinejoin="round"
                  />
                </svg>
              </span>
              <span
                style={{
                  color: "#666",
                  fontSize: 16,
                  textAlign: "center",
                }}
              >
                Use the <b>New</b> or <b>Import</b> buttons above to create or
                load a project.
              </span>
            </Stack>
          )
        ) : (
          <Skeleton variant="rectangular" width="100%" height={500} />
        )}
      </Stack>
    </Container>
    </>
  );
}
