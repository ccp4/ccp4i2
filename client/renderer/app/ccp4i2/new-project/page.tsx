"use client";
import { NewProjectContent } from "../../../components/new-project-content";
import { NavigationShortcutsProvider } from "../../../providers/navigation-shortcuts-provider";
import CCP4i2TopBar from "../../../components/ccp4i2-topbar";

export default function NewProjectPage() {
  return (
    <NavigationShortcutsProvider>
      <CCP4i2TopBar title="New Project" showBackButton backPath="/ccp4i2" />
      <NewProjectContent />
    </NavigationShortcutsProvider>
  );
}
