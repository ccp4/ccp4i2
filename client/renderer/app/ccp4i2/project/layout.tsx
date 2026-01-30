"use client";
import { PropsWithChildren } from "react";
import { CootProvider } from "../../../providers/coot-provider";
import { RunningProcessesProvider } from "../../../providers/running-processes";

export default function ProjectLayout(props: PropsWithChildren) {
  return (
    <CootProvider>
      <RunningProcessesProvider>
        {/* Children components will be rendered here */}
        {props.children}
      </RunningProcessesProvider>
    </CootProvider>
  );
}
