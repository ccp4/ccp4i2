"use client";

import dynamic from "next/dynamic";

const SessionClient = dynamic(() => import("./session-client"), { ssr: false });

export default function SessionPage() {
  return <SessionClient />;
}
