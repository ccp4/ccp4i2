"use client";
import { PropsWithChildren } from "react";
import { CootProvider } from "../../providers/coot-provider";
import { RdkitProvider } from "../../providers/rdkit-provider";
import { RunningProcessesProvider } from "../../providers/running-processes";

export default function ProjectLayout(props: PropsWithChildren) {
  return (
    <CootProvider>
      <RdkitProvider>
        <RunningProcessesProvider>
          {/* Children components will be rendered here */}
          {props.children}
        </RunningProcessesProvider>
      </RdkitProvider>
    </CootProvider>
  );
}
