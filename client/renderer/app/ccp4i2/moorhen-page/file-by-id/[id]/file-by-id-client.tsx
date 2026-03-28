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
import { Suspense } from "react";
import { useParams, useSearchParams } from "next/navigation";
import MoorhenLoader from "../../../../../components/moorhen/client-side-moorhen-loader";
import { ClientStoreProvider } from "../../../../../providers/client-store-provider";
import { useApi } from "../../../../../api";

// Inner component that uses useSearchParams (requires Suspense boundary)
function FileByIdPageContent() {
  const params = useParams();
  const id = params?.id as string | undefined;
  const searchParams = useSearchParams();
  const viewParam = searchParams?.get("view");

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
const FileByIdClient = () => {
  return (
    <Suspense fallback={<div>Loading...</div>}>
      <FileByIdPageContent />
    </Suspense>
  );
};

export default FileByIdClient;
