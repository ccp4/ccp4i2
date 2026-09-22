"use client";
import { use } from "react";
import { EditProjectContent } from "@/components/edit-project-content";
import { useTopBar } from "@/providers/top-bar-context";

export default function EditProjectPage({
  params,
}: {
  params: Promise<{ id: string }>;
}) {
  const { id } = use(params);
  useTopBar({ title: "Edit Project" });
  return <EditProjectContent projectId={parseInt(id)} />;
}
