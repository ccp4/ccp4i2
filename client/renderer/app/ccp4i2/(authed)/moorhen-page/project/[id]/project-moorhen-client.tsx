"use client";

import { Suspense, useEffect, useMemo } from "react";
import { useParams, useRouter, useSearchParams } from "next/navigation";
import MoorhenLoader from "@/components/moorhen/client-side-moorhen-loader";
import { ClientStoreProvider } from "@/providers/client-store-provider";
import { sceneSourceFromParams } from "@/lib/scene-source";
import { listenAsViewer } from "@/lib/viewer-channel";

// Project-scoped Moorhen page: the project (pk from the route) provides context
// for scene authoring (the prompt manifest) and job+param file resolution. It
// loads no file or job unless asked: ``?path=<relative>`` (repeatable) opens
// files or folders by their path in the project, as the directory browser's
// "Open in Moorhen" does. The route `[id]` is the project primary key,
// matching /ccp4i2/project/[id].
function ProjectMoorhenContent() {
  const params = useParams();
  const id = params?.id as string | undefined;
  const searchParams = useSearchParams();
  const viewParam = searchParams?.get("view");
  const projectId = id ? parseInt(id, 10) : null;
  const pathsKey = (searchParams?.getAll("path") ?? []).join("\n");
  const filePaths = useMemo(() => (pathsKey ? pathsKey.split("\n") : []), [pathsKey]);
  // ?sceneJob=&sceneMethod=&sceneArgs=: a scene the job builds on request
  // (a report's "open in Moorhen" link; lib/scene-source.ts).
  const sceneKey = ["sceneJob", "sceneMethod", "sceneArgs"].map((k) => searchParams?.get(k) ?? "").join("\n");
  const sceneSource = useMemo(
    () => (searchParams ? sceneSourceFromParams(searchParams) : null),
    // eslint-disable-next-line react-hooks/exhaustive-deps
    [sceneKey],
  );

  // Take links other windows send this project's viewer (lib/viewer-channel):
  // replacing the URL changes the scene source, which the open viewer applies.
  const router = useRouter();
  useEffect(() => listenAsViewer((href) => router.replace(href)), [router]);

  return (
    <ClientStoreProvider>
      <MoorhenLoader fileIds={[]} viewParam={viewParam} projectId={projectId} filePaths={filePaths} sceneSource={sceneSource} />
    </ClientStoreProvider>
  );
}

const ProjectMoorhenClient = () => (
  <Suspense fallback={<div>Loading…</div>}>
    <ProjectMoorhenContent />
  </Suspense>
);

export default ProjectMoorhenClient;
