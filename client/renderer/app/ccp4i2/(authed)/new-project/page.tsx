"use client";
import { NewProjectContent } from "@/components/new-project-content";
import { useTopBar } from "@/providers/top-bar-context";

export default function NewProjectPage() {
  useTopBar({ title: "New Project" });
  return <NewProjectContent />;
}
