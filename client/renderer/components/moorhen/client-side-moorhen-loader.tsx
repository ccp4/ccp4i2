/*
 * Copyright (C) 2025 Newcastle University
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
import dynamic from "next/dynamic";
import { MoorhenWrapperProps } from "./moorhen-wrapper";

const MoorhenClient = dynamic(() => import("./client-side-moorhen-component"), {
  ssr: false,
});

export default function DynamicMoorhenClient(props: MoorhenWrapperProps) {
  return <MoorhenClient {...props} />;
}
