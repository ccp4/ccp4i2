"use client";
import { use } from "react";
import { EditProjectContent } from "@/components/edit-project-content";
import CCP4i2TopBar from "@/components/ccp4i2-topbar";

export default function EditProjectPage({
  params,
}: {
  params: Promise<{ id: string }>;
}) {
  const { id } = use(params);
  return (
    <>
      <CCP4i2TopBar title="Edit Project" showBackButton backPath="/ccp4i2" />
      <EditProjectContent projectId={parseInt(id)} />
    </>
  );
}
