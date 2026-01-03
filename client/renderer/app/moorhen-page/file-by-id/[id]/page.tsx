"use client";
import { useParams } from "next/navigation";
import MoorhenLoader from "../../../../components/moorhen/client-side-moorhen-loader";
import { ClientStoreProvider } from "../../../../providers/client-store-provider";
import { useApi } from "../../../../api";

const FileByIdPage = () => {
  const { id } = useParams();
  const api = useApi();
  // Use centralized API hook for file fetching (data unused but kept for future use)
  const { data: file } = api.file(id as string);
  return (
    <ClientStoreProvider>
      <MoorhenLoader fileIds={id ? [parseInt(id as string)] : []} />
    </ClientStoreProvider>
  );
};
export default FileByIdPage;
