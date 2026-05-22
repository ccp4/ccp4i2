"use client";
import { use } from "react";
import CCP4i2TopBar from "@/components/ccp4i2-topbar";
import { NavigationShortcutsProvider } from "@/providers/navigation-shortcuts-provider";
import { EditProjectContent } from "@/components/edit-project-content";

export default function EditProjectPage({
  params,
}: {
  params: Promise<{ id: string }>;
}) {
  const { id } = use(params);
  return (
    <NavigationShortcutsProvider>
      <CCP4i2TopBar title="Edit Project" showBackButton />
      <EditProjectContent id={id} />
    </NavigationShortcutsProvider>
  );
}
