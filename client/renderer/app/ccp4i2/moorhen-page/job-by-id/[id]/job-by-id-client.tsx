/*
 * Copyright (C) 2026 Newcastle University
 *
 * This file is part of CCP4i2.
 *
 * CCP4i2 is free software: you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License version 3,
 * modified in accordance with the provisions of the license to address
 * the requirements of UK law.
 *
 * See https://www.ccp4.ac.uk/ccp4license.php for details.
 */
"use client";
import { Suspense, useMemo, useRef } from "react";
import { useParams, useSearchParams } from "next/navigation";
import MoorhenLoader from "../../../../../components/moorhen/client-side-moorhen-loader";
import { ClientStoreProvider } from "../../../../../providers/client-store-provider";
import { useJob, useProject } from "../../../../../utils";

// Inner component that uses useSearchParams (requires Suspense boundary)
function JobByIdPageContent() {
  const params = useParams();
  const id = params?.id as string | undefined;
  const searchParams = useSearchParams();
  const viewParam = searchParams?.get("view");

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
      <MoorhenLoader fileIds={fileIds} viewParam={viewParam} jobId={job?.id ?? null} />
    </ClientStoreProvider>
  );
}

// Outer component with Suspense boundary for useSearchParams
const JobByIdClient = () => {
  return (
    <Suspense fallback={<div>Loading...</div>}>
      <JobByIdPageContent />
    </Suspense>
  );
};

export default JobByIdClient;
