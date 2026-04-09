"use client";

import { ClientStoreProvider } from "../../providers/client-store-provider";
import MoorhenWrapper, { MoorhenWrapperProps } from "./moorhen-wrapper";

const ClientSideMoorhenComponent: React.FC<MoorhenWrapperProps> = (props) => {
  return (
    <ClientStoreProvider>
      <MoorhenWrapper {...props} />
    </ClientStoreProvider>
  );
};

export default ClientSideMoorhenComponent;
