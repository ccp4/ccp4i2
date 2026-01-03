"use client";
import { PropsWithChildren } from "react";
import { CootProvider } from "../../providers/coot-provider";
import { RdkitProvider } from "../../providers/rdkit-provider";
import { RunningProcessesProvider } from "../../providers/running-processes";
import { DraggableContext } from "../../providers/draggable-context";
import { NavigationShortcutsProvider } from "../../providers/navigation-shortcuts-provider";
import { FileSystemFileBrowserProvider } from "../../providers/file-system-file-browser-context";
import Stack from "@mui/material/Stack/Stack";
import { FilePreviewProvider } from "../../providers/file-preview-context";
import { JobMenuProvider } from "../../providers/job-context-menu";
import { JobTabProvider } from "../../providers/job-tab-provider";
import { FileMenuProvider } from "../../providers/file-context-menu";
import MenuBar from "../../components/menu-bar";
import { Panel, PanelGroup } from "react-resizable-panels";

export default function JobLayout(props: PropsWithChildren) {
  return (
    <CootProvider>
      <RdkitProvider>
        <RunningProcessesProvider>
          <DraggableContext>
            <NavigationShortcutsProvider>
              <FileSystemFileBrowserProvider>
                <FilePreviewProvider>
                  <JobMenuProvider>
                    <JobTabProvider>
                      <FileMenuProvider>
                        <MenuBar />
                        {/* Children components will be rendered here */}
                        <Stack
                          spacing={2}
                          sx={{
                            height: "calc(100vh - 4rem)",
                            paddingTop: "1rem",
                            width: "100%",
                          }}
                        >
                          <PanelGroup direction="horizontal">
                            <Panel>{props.children}</Panel>
                          </PanelGroup>
                        </Stack>{" "}
                      </FileMenuProvider>
                    </JobTabProvider>
                  </JobMenuProvider>
                </FilePreviewProvider>
              </FileSystemFileBrowserProvider>
            </NavigationShortcutsProvider>
          </DraggableContext>
        </RunningProcessesProvider>
      </RdkitProvider>
    </CootProvider>
  );
}
