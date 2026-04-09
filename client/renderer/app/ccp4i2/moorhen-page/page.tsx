"use client";

import dynamic from "next/dynamic";

// Dynamically import the entire moorhen page client component with SSR disabled.
// Moorhen's UMD bundle references HTMLElement at module load time,
// which is not available during Next.js server-side prerendering.
const MoorhenPageClient = dynamic(
  () => import("./moorhen-page-client"),
  { ssr: false }
);

export default function MoorhenPage() {
  return <MoorhenPageClient />;
}
