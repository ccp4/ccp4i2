"use client";
import { useParams } from "next/dist/client/components/navigation";
import { JobView } from "../../../../../components/job-view";

export default function JobPage() {
  const params = useParams();
  const { jobid } = params as { jobid: string };
  return <JobView jobid={parseInt(jobid)} />;
}
