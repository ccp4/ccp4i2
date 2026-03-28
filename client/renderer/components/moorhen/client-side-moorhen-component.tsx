/*
 * Copyright (C) 2025-2026 Newcastle University
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
