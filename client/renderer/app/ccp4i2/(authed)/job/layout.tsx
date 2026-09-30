"use client";
import { PropsWithChildren } from "react";
import { DraggableContext } from "@/providers/draggable-context";
import { FileSystemFileBrowserProvider } from "@/providers/file-system-file-browser-context";
import { FilePreviewProvider } from "@/providers/file-preview-context";
import { JobMenuProvider } from "@/providers/job-context-menu";
import { JobTabProvider } from "@/providers/job-tab-provider";
import { FileMenuProvider } from "@/providers/file-context-menu";
import { Panel, PanelGroup } from "react-resizable-panels";

export default function JobLayout(props: PropsWithChildren) {
  return (
    <DraggableContext>
      <FileSystemFileBrowserProvider>
        <FilePreviewProvider>
          <JobMenuProvider>
            <JobTabProvider>
              <FileMenuProvider>
                {/* Explicit ids, as in project/[id]/layout.tsx: generated
                    ones can differ between server and client and crash the page. */}
                <PanelGroup id="job-panels" direction="horizontal">
                  <Panel id="job-content">{props.children}</Panel>
                </PanelGroup>
              </FileMenuProvider>
            </JobTabProvider>
          </JobMenuProvider>
        </FilePreviewProvider>
      </FileSystemFileBrowserProvider>
    </DraggableContext>
  );
}
