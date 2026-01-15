"use client";
import { PropsWithChildren } from "react";
import { CootProvider } from "../../../providers/coot-provider";
import { RDKitProvider } from "../../../providers/rdkit-provider";
import { RunningProcessesProvider } from "../../../providers/running-processes";

export default function ProjectLayout(props: PropsWithChildren) {
  return (
    <CootProvider>
      <RDKitProvider>
        <RunningProcessesProvider>
          {/* Children components will be rendered here */}
          {props.children}
        </RunningProcessesProvider>
      </RDKitProvider>
    </CootProvider>
  );
}
