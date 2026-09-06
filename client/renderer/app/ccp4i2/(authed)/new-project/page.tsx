"use client";
import { NewProjectContent } from "@/components/new-project-content";
import CCP4i2TopBar from "@/components/ccp4i2-topbar";

export default function NewProjectPage() {
  return (
    <>
      <CCP4i2TopBar title="New Project" showBackButton backPath="/ccp4i2" />
      <NewProjectContent />
    </>
  );
}
