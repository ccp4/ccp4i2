"use client";

import dynamic from "next/dynamic";

const CampaignPageClient = dynamic(
  () => import("./campaign-page-client"),
  { ssr: false }
);

export default function CampaignPage() {
  return <CampaignPageClient />;
}
