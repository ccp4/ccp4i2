"use client";
import { PropsWithChildren, use, useEffect, useState } from "react";
import {
  Paper,
  Stack,
  Tab,
  Tabs,
  useMediaQuery,
  useTheme,
  Box,
} from "@mui/material";
import { Panel, PanelGroup, PanelResizeHandle } from "react-resizable-panels";
import { useCCP4i2Window } from "../../../app-context";
import { useApi } from "../../../api";
import { CCP4i2DirectoryViewer } from "../../../components/directory-viewer";
import { useProject } from "../../../utils";
import { ClassicJobList } from "../../../components/classic-jobs-list";
import { DraggableContext } from "../../../providers/draggable-context";
import { FilePreviewProvider } from "../../../providers/file-preview-context";
import { JobMenuProvider } from "../../../providers/job-context-menu";
import { FileMenuProvider } from "../../../providers/file-context-menu";
import MenuBar from "../../../components/menu-bar";
import { NavigationShortcutsProvider } from "../../../providers/navigation-shortcuts-provider";
import { FileSystemFileBrowserProvider } from "../../../providers/file-system-file-browser-context";
import { JobTabProvider } from "../../../providers/job-tab-provider";
import { ParameterChangeIntentProvider } from "../../../providers/parameter-change-intent-provider";

export interface ProjectLayoutProps extends PropsWithChildren {
  params: Promise<{ id: string }>; // Removed jobid since it's not available at this route level
}

export default function ProjectLayout(props: ProjectLayoutProps) {
  const { setProjectId, setJobPanelSize } = useCCP4i2Window();
  const api = useApi();
  const theme = useTheme();
  const isMobile = useMediaQuery(theme.breakpoints.down("md")); // Mobile: < 900px

  // For desktop left panel tabs
  const [leftTabValue, setLeftTabValue] = useState(0);
  // For mobile main tabs (left panel + right panel)
  const [mobileTabValue, setMobileTabValue] = useState(0);

  const { id } = use(props.params);
  const { project } = useProject(parseInt(id));

  useEffect(() => {
    const asyncFunc = async () => {
      if (project && setProjectId) {
        setProjectId(project.id);
      }
    };
    asyncFunc();
  }, [project, setProjectId]);

  const handleLeftTabChange = (event: React.SyntheticEvent, value: number) => {
    setLeftTabValue(value);
  };

  const handleMobileTabChange = (
    event: React.SyntheticEvent,
    value: number
  ) => {
    setMobileTabValue(value);
  };

  return (
    <DraggableContext>
      <NavigationShortcutsProvider>
        <FileSystemFileBrowserProvider>
          <FilePreviewProvider>
            <JobMenuProvider>
              <JobTabProvider>
                <FileMenuProvider>
                  <ParameterChangeIntentProvider>
                    <Stack
                      spacing={2}
                      sx={{
                        height: "calc(100vh - 4rem)",
                        paddingTop: "1rem",
                        width: "100%",
                      }}
                    >
                      <MenuBar />

                      {isMobile ? (
                        // Mobile: Tabbed interface
                        <Box sx={{ height: "calc(100vh - 10rem)" }}>
                          <Paper sx={{ height: "100%" }}>
                            <Tabs
                              value={mobileTabValue}
                              onChange={handleMobileTabChange}
                              variant="fullWidth"
                            >
                              <Tab value={0} label="Jobs" />
                              <Tab value={1} label="Directory" />
                              <Tab value={2} label="Content" />
                            </Tabs>
                            <Box
                              sx={{
                                height: "calc(100% - 48px)",
                                overflow: "auto",
                                // Theme-aware scrollbar styling
                                scrollbarColor: `${theme.palette.action.disabled} transparent`,
                                scrollbarWidth: "thin",
                                "&::-webkit-scrollbar": {
                                  width: 8,
                                },
                                "&::-webkit-scrollbar-track": {
                                  background: "transparent",
                                },
                                "&::-webkit-scrollbar-thumb": {
                                  backgroundColor: theme.palette.action.disabled,
                                  borderRadius: 4,
                                },
                              }}
                            >
                              {mobileTabValue === 0 && project && (
                                <ClassicJobList projectId={project.id} />
                              )}
                              {mobileTabValue === 1 && project && (
                                <CCP4i2DirectoryViewer projectId={project.id} />
                              )}
                              {mobileTabValue === 2 && (
                                <Box sx={{ height: "100%" }}>
                                  {props.children}
                                </Box>
                              )}
                            </Box>
                          </Paper>
                        </Box>
                      ) : (
                        // Desktop: Side-by-side panels
                        <PanelGroup direction="horizontal">
                          <Panel defaultSize={30} minSize={20}>
                            <Paper
                              sx={{
                                overflowY: "auto",
                                height: "calc(100vh - 10rem)",
                                // Theme-aware scrollbar styling
                                scrollbarColor: `${theme.palette.action.disabled} transparent`,
                                scrollbarWidth: "thin",
                                "&::-webkit-scrollbar": {
                                  width: 8,
                                },
                                "&::-webkit-scrollbar-track": {
                                  background: "transparent",
                                },
                                "&::-webkit-scrollbar-thumb": {
                                  backgroundColor: theme.palette.action.disabled,
                                  borderRadius: 4,
                                },
                              }}
                            >
                              <Tabs
                                value={leftTabValue}
                                onChange={handleLeftTabChange}
                                variant="fullWidth"
                              >
                                <Tab value={0} label="Job list" />
                                <Tab value={1} label="Project directory" />
                              </Tabs>
                              {leftTabValue === 0 && project && (
                                <ClassicJobList projectId={project.id} />
                              )}
                              {leftTabValue === 1 && project && (
                                <CCP4i2DirectoryViewer projectId={project.id} />
                              )}
                            </Paper>
                          </Panel>
                          <PanelResizeHandle
                            style={{
                              width: 10,
                              backgroundColor: "transparent",
                              display: "flex",
                              alignItems: "center",
                              justifyContent: "center",
                              cursor: "col-resize",
                            }}
                          >
                            <div
                              style={{
                                width: 4,
                                height: "50%",
                                backgroundColor: "gray",
                                borderRadius: 2,
                              }}
                            />
                          </PanelResizeHandle>
                          <Panel
                            defaultSize={70}
                            minSize={20}
                            onResize={(size) =>
                              setJobPanelSize && setJobPanelSize(size)
                            }
                          >
                            {props.children}
                          </Panel>
                        </PanelGroup>
                      )}
                    </Stack>
                  </ParameterChangeIntentProvider>
                </FileMenuProvider>
              </JobTabProvider>
            </JobMenuProvider>
          </FilePreviewProvider>
        </FileSystemFileBrowserProvider>
      </NavigationShortcutsProvider>
    </DraggableContext>
  );
}
