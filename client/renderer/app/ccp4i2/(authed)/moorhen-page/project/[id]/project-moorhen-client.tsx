"use client";

import { Suspense } from "react";
import { useParams, useSearchParams } from "next/navigation";
import MoorhenLoader from "@/components/moorhen/client-side-moorhen-loader";
import { ClientStoreProvider } from "@/providers/client-store-provider";

// Project-scoped Moorhen page: the project (pk from the route) provides context
// for scene authoring (the prompt manifest) and job+param file resolution, but
// no specific file or job is loaded — scenes do the loading. The route `[id]` is
// the project primary key, matching /ccp4i2/project/[id].
function ProjectMoorhenContent() {
  const params = useParams();
  const id = params?.id as string | undefined;
  const searchParams = useSearchParams();
  const viewParam = searchParams?.get("view");
  const projectId = id ? parseInt(id, 10) : null;

  return (
    <ClientStoreProvider>
      <MoorhenLoader fileIds={[]} viewParam={viewParam} projectId={projectId} />
    </ClientStoreProvider>
  );
}

const ProjectMoorhenClient = () => (
  <Suspense fallback={<div>Loading…</div>}>
    <ProjectMoorhenContent />
  </Suspense>
);

export default ProjectMoorhenClient;
