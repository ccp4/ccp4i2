"use client";
import { Suspense } from "react";
import { useParams, useSearchParams } from "next/navigation";
import MoorhenLoader from "../../../../../components/moorhen/client-side-moorhen-loader";
import { ClientStoreProvider } from "../../../../../providers/client-store-provider";
import { useApi } from "../../../../../api";

// Inner component that uses useSearchParams (requires Suspense boundary)
function FileByIdPageContent() {
  const { id } = useParams();
  const searchParams = useSearchParams();
  const viewParam = searchParams.get("view");

  const api = useApi();
  // Use centralized API hook for file fetching (data unused but kept for future use)
  const { data: file } = api.file(id as string);

  return (
    <ClientStoreProvider>
      <MoorhenLoader
        fileIds={id ? [parseInt(id as string)] : []}
        viewParam={viewParam}
      />
    </ClientStoreProvider>
  );
}

// Outer component with Suspense boundary for useSearchParams
const FileByIdPage = () => {
  return (
    <Suspense fallback={<div>Loading...</div>}>
      <FileByIdPageContent />
    </Suspense>
  );
};

export default FileByIdPage;
