"use client";

import MoorhenLoader from "../../../components/moorhen/client-side-moorhen-loader";
import MoorhenWrapper from "../../../components/moorhen/moorhen-wrapper";
import { ClientStoreProvider } from "../../../providers/client-store-provider";

export default function MoorhenPage() {
  return (
    <ClientStoreProvider>
      <MoorhenLoader fileIds={[]} />
    </ClientStoreProvider>
  );
}
