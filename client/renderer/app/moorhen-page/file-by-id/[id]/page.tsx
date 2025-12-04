"use client";
import { useParams } from "next/navigation";
import useSWR from "swr";
import MoorhenWrapper from "../../../../components/moorhen/moorhen-wrapper";
import MoorhenLoader from "../../../../components/moorhen/client-side-moorhen-loader";
import { ClientStoreProvider } from "../../../../providers/client-store-provider";
import { swrFetcher } from "../../../../api-fetch";

const FileByIdPage = () => {
  const { id } = useParams();
  const { data: file } = useSWR(`/api/proxy/files/${id}`, swrFetcher);
  return (
    <ClientStoreProvider>
      <MoorhenLoader fileIds={id ? [parseInt(id as string)] : []} />
    </ClientStoreProvider>
  );
};
export default FileByIdPage;
