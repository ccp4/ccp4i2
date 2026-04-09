"use client";

import dynamic from "next/dynamic";

const JobByIdClient = dynamic(
  () => import("./job-by-id-client"),
  { ssr: false }
);

export default function JobByIdPage() {
  return <JobByIdClient />;
}
