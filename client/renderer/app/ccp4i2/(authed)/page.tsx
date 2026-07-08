"use client";
import { Box, Container, Skeleton, Stack } from "@mui/material";
import ProjectsToolbar from "@/components/projects-toolbar";
import ProjectsTable from "@/components/projects-table";
import CCP4i2TopBar from "@/components/ccp4i2-topbar";
import { LaunchGate } from "@/components/launch-gate";
import { WelcomeChooser } from "@/components/welcome-chooser";
import { useApi } from "@/api";
import { Project } from "@/types/models";

function ProjectsHome() {
  const api = useApi();
  // Note: second param is refreshInterval, not timeout. Set to 0 to disable auto-refresh.
  const { data: projects } = api.get<Project[]>("projects");

  // Empty database => first-project welcome. This is keyed on live state, not a
  // one-shot flag, so it also covers a failed migration or a rolled-back
  // projects directory that leaves the user back at zero projects.
  const isEmpty = projects && projects.length === 0;

  return (
    <Stack sx={{ height: "100vh", overflow: "hidden" }}>
      <CCP4i2TopBar title="CCP4i2 Projects" />
      {isEmpty ? (
        <Box sx={{ flex: 1, overflow: "auto" }}>
          <WelcomeChooser />
        </Box>
      ) : (
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
              <Box sx={{ flex: 1, overflow: "hidden" }}>
                <ProjectsTable />
              </Box>
            ) : (
              <Skeleton variant="rectangular" width="100%" height={500} />
            )}
          </Stack>
        </Container>
      )}
    </Stack>
  );
}

export default function ProjectsPage() {
  return (
    <LaunchGate>
      <ProjectsHome />
    </LaunchGate>
  );
}
