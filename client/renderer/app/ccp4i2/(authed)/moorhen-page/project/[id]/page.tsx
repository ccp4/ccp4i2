"use client";

import dynamic from "next/dynamic";

const ProjectMoorhenClient = dynamic(
  () => import("./project-moorhen-client"),
  { ssr: false }
);

export default function ProjectMoorhenPage() {
  return <ProjectMoorhenClient />;
}
