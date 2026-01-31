"use client";
import { PropsWithChildren } from "react";
import { RunningProcessesProvider } from "../../../providers/running-processes";

export default function ProjectLayout(props: PropsWithChildren) {
  return (
    <RunningProcessesProvider>
      {/* Children components will be rendered here */}
      {props.children}
    </RunningProcessesProvider>
  );
}
