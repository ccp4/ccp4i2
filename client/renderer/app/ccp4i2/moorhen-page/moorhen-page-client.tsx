"use client";

import MoorhenLoader from "../../../components/moorhen/client-side-moorhen-loader";
import { ClientStoreProvider } from "../../../providers/client-store-provider";

export default function MoorhenPageClient() {
  return (
    <ClientStoreProvider>
      <MoorhenLoader fileIds={[]} />
    </ClientStoreProvider>
  );
}
