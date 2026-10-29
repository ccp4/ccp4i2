"use client";
import { PropsWithChildren, use, useEffect, useState } from "react";
import { usePathname } from "next/navigation";
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
import { useCCP4i2Window } from "@/app-context";
import { CCP4i2DirectoryViewer } from "@/components/directory-viewer";
import { useProject } from "@/utils";
import { ClassicJobList } from "@/components/classic-jobs-list";
import { PanelErrorBoundary } from "@/components/panel-error-boundary";
import { DraggableContext } from "@/providers/draggable-context";
import { FilePreviewProvider } from "@/providers/file-preview-context";
import { JobMenuProvider } from "@/providers/job-context-menu";
import { FileMenuProvider } from "@/providers/file-context-menu";
import { FileSystemFileBrowserProvider } from "@/providers/file-system-file-browser-context";
import { JobTabProvider } from "@/providers/job-tab-provider";

export interface ProjectLayoutProps extends PropsWithChildren {
  params: Promise<{ id: string }>; // Removed jobid since it's not available at this route level
}

export default function ProjectLayout(props: ProjectLayoutProps) {
  const { setProjectId } = useCCP4i2Window();
  const theme = useTheme();
  const isMobile = useMediaQuery(theme.breakpoints.down("md")); // Mobile: < 900px

  // For desktop left panel tabs
  const [leftTabValue, setLeftTabValue] = useState(0);
  // For mobile main tabs (left panel + right panel)
  const [mobileTabValue, setMobileTabValue] = useState(0);

  const { id } = use(props.params);
  const { project } = useProject(parseInt(id));
  // A page that failed to draw gets a fresh start on navigation.
  const pathname = usePathname();

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
      <FileSystemFileBrowserProvider>
        <FilePreviewProvider>
          <JobMenuProvider>
            <JobTabProvider>
              <FileMenuProvider>
                <Stack
                    spacing={2}
                    sx={{
                      flex: 1,
                      minHeight: 0,
                      width: "100%",
                      pt: 2,
                    }}
                  >
                    {isMobile ? (
                      // Mobile: Tabbed interface
                      <Box sx={{ flex: 1, minHeight: 0 }}>
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
                              <PanelErrorBoundary name="job list" resetKeys={[project.id]}>
                                <ClassicJobList projectId={project.id} />
                              </PanelErrorBoundary>
                            )}
                            {mobileTabValue === 1 && project && (
                              <PanelErrorBoundary name="project directory" resetKeys={[project.id]}>
                                <CCP4i2DirectoryViewer projectId={project.id} />
                              </PanelErrorBoundary>
                            )}
                            {mobileTabValue === 2 && (
                              <Box sx={{ height: "100%" }}>
                                <PanelErrorBoundary name="page" resetKeys={[pathname]}>
                                  {props.children}
                                </PanelErrorBoundary>
                              </Box>
                            )}
                          </Box>
                      </Box>
                    ) : (
                      // Desktop: Side-by-side panels
                      // Explicit ids: the generated ones can differ between the
                      // server render and the client, and the panels then look
                      // up a group the DOM does not have ("No group found for
                      // id"), which crashes the whole page.
                      <PanelGroup id="project-panels" direction="horizontal" style={{flex: "auto"}}>
                        <Panel id="project-sidebar" order={1} defaultSize={30} minSize={20}>
                          <Stack sx={{ height: "100%" }}>
                            <Tabs
                              value={leftTabValue}
                              onChange={handleLeftTabChange}
                              variant="fullWidth"
                            >
                              <Tab value={0} label="Job list" />
                              <Tab value={1} label="Project directory" />
                            </Tabs>
                            {leftTabValue === 0 && project && (
                              <PanelErrorBoundary name="job list" resetKeys={[project.id]}>
                                <ClassicJobList projectId={project.id} />
                              </PanelErrorBoundary>
                            )}
                            {leftTabValue === 1 && project && (
                              <PanelErrorBoundary name="project directory" resetKeys={[project.id]}>
                                <CCP4i2DirectoryViewer projectId={project.id} />
                              </PanelErrorBoundary>
                            )}
                          </Stack>
                        </Panel>
                        <PanelResizeHandle
                          id="project-panels-handle"
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
                        <Panel id="project-content" order={2} defaultSize={70} minSize={20}>
                          <PanelErrorBoundary name="page" resetKeys={[pathname]}>
                            {props.children}
                          </PanelErrorBoundary>
                        </Panel>
                      </PanelGroup>
                    )}
                </Stack>
              </FileMenuProvider>
            </JobTabProvider>
          </JobMenuProvider>
        </FilePreviewProvider>
      </FileSystemFileBrowserProvider>
    </DraggableContext>
  );
}
