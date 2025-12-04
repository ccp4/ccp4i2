"use client";

import { ClientStoreProvider } from "../../providers/client-store-provider";
import MoorhenWrapper, { MoorhenWrapperProps } from "./moorhen-wrapper";
import { useCCP4i2Window } from "../../app-context";
import { useStore } from "react-redux";

const ClientSideMoorhenComponent: React.FC<MoorhenWrapperProps> = (props) => {
  const { cootModule } = useCCP4i2Window();
  const store = useStore();
  return (
    <ClientStoreProvider>
      {cootModule && store && <MoorhenWrapper {...props} />}
    </ClientStoreProvider>
  );
};

export default ClientSideMoorhenComponent;
