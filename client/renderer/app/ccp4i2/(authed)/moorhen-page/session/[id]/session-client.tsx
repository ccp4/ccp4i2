"use client";
/**
 * The recorded Moorhen session window. The route parameter is the job the
 * session belongs to: the window loads that job's inputs from the server's
 * load plan and saves back into that job. Sibling of job-by-id, which shows
 * a finished job's outputs; this one works on a running job's inputs.
 */
import { Suspense } from "react";
import { useParams, useSearchParams } from "next/navigation";
import MoorhenLoader from "@/components/moorhen/client-side-moorhen-loader";
import { ClientStoreProvider } from "@/providers/client-store-provider";

function SessionPageContent() {
  const params = useParams();
  const id = parseInt(params?.id as string);
  const searchParams = useSearchParams();
  const viewParam = searchParams?.get("view");
  if (!Number.isFinite(id)) return <div>No job given</div>;
  return (
    <ClientStoreProvider>
      <MoorhenLoader fileIds={[]} viewParam={viewParam} jobId={id} sessionJobId={id} />
    </ClientStoreProvider>
  );
}

const SessionClient = () => (
  <Suspense fallback={<div>Loading...</div>}>
    <SessionPageContent />
  </Suspense>
);

export default SessionClient;
