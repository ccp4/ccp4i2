"use client";
import { Suspense, useMemo, useRef } from "react";
import { useParams, useSearchParams } from "next/navigation";
import MoorhenLoader from "../../../../../components/moorhen/client-side-moorhen-loader";
import { ClientStoreProvider } from "../../../../../providers/client-store-provider";
import { useJob, useProject } from "../../../../../utils";

// Inner component that uses useSearchParams (requires Suspense boundary)
function JobByIdPageContent() {
  const { id } = useParams();
  const searchParams = useSearchParams();
  const viewParam = searchParams.get("view");

  const { job } = useJob(parseInt(id as string));
  const { files } = useProject(job?.project);

  // Use ref to store the fileIds once determined for a valid job
  const stableFileIds = useRef<number[]>([]);
  const lastValidJobId = useRef<number | null>(null);

  // Calculate fileIds and stabilize them once determined for a valid job
  const fileIds = useMemo(() => {
    if (!job?.id || !files) {
      return stableFileIds.current; // Return existing stable fileIds if job/files not ready
    }

    // If this is a new job, reset and calculate new fileIds
    if (lastValidJobId.current !== job.id) {
      const newFileIds = files
        .filter((file) => file.job === job.id)
        .map((file) => file.id);

      stableFileIds.current = newFileIds;
      lastValidJobId.current = job.id;
    }

    return stableFileIds.current;
  }, [job?.id, files]);

  return (
    <ClientStoreProvider>
      <MoorhenLoader fileIds={fileIds} viewParam={viewParam} />
    </ClientStoreProvider>
  );
}

// Outer component with Suspense boundary for useSearchParams
const JobByIdPage = () => {
  return (
    <Suspense fallback={<div>Loading...</div>}>
      <JobByIdPageContent />
    </Suspense>
  );
};

export default JobByIdPage;
