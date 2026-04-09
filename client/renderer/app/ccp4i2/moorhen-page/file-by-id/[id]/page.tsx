"use client";

import dynamic from "next/dynamic";

const FileByIdClient = dynamic(
  () => import("./file-by-id-client"),
  { ssr: false }
);

export default function FileByIdPage() {
  return <FileByIdClient />;
}
