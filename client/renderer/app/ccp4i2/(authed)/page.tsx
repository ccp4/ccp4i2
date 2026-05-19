"use client";
import { Box, Container, Skeleton, Stack, Typography } from "@mui/material";
import ProjectsToolbar from "@/components/projects-toolbar";
import ProjectsTable from "@/components/projects-table";
import CCP4i2TopBar from "@/components/ccp4i2-topbar";
import { useApi } from "@/api";
import { Project } from "@/types/models";

export default function ProjectsPage() {
  const api = useApi();
  // Note: second param is refreshInterval, not timeout. Set to 0 to disable auto-refresh.
  const { data: projects } = api.get<Project[]>("projects");

  return (
    <Stack
      sx={{
        height: "100vh",
        overflow: "hidden",
      }}
    >
      <CCP4i2TopBar title="CCP4i2 Projects" />
      <Container
        sx={{
          flex: 1,
          display: "flex",
          flexDirection: "column",
          overflow: "hidden",
          py: 2,
        }}
      >
        <Stack gap={2} sx={{ flex: 1, overflow: "hidden" }}>
          <ProjectsToolbar />
          {projects ? (
            projects.length > 0 ? (
              <Box sx={{ flex: 1, overflow: "hidden" }}>
                <ProjectsTable />
              </Box>
            ) : (
              <Stack alignItems="center" gap={2} sx={{ py: 8 }}>
                <Typography variant="h6">No projects found</Typography>
                <Typography variant="body1" color="textSecondary">
                  Use the <b>New</b> or <b>Import</b> buttons above
                  to create or load a project.
                </Typography>
              </Stack>
            )
          ) : (
            <Skeleton variant="rectangular" width="100%" height={500} />
          )}
        </Stack>
      </Container>
    </Stack>
  );
}
