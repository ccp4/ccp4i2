/*
 * Copyright (C) 2026 Newcastle University
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

import dynamic from "next/dynamic";

const CampaignPageClient = dynamic(
  () => import("./campaign-page-client"),
  { ssr: false }
);

export default function CampaignPage() {
  return <CampaignPageClient />;
}
