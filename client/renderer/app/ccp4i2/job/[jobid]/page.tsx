"use client";
import { JobView } from "../../../../components/job-view";
import { useParams } from "next/navigation";

export default function JobPage() {
  const params = useParams();
  const { jobid } = params as { jobid: string };
  return <JobView jobid={parseInt(jobid)} />;
}
