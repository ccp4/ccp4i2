"use client";
import { useEffect, useState } from "react";
import { useRouter } from "next/navigation";
import {
  Button,
  CircularProgress,
  Container,
  Skeleton,
  Stack,
  Typography,
} from "@mui/material";
import { useApi } from "../api";
import { apiPost, apiDelete } from "../api-fetch";
import { Project, ProjectTag } from "../types/models";
import EditTags from "./edit-tags";
import { ProjectNameField, getProjectNameError } from "./project-name-field";
import { ProjectDirectoryField } from "./project-directory-field";
import { useProjectDirectory } from "../hooks/use-project-directory";

interface EditProjectContentProps {
  id: string;
}

export function EditProjectContent({ id }: EditProjectContentProps) {
  const api = useApi();
  const router = useRouter();
  const { data: project, isLoading } = api.get<Project>(`projects/${id}`);
  const { data: projects } = api.get<Project[]>("projects");
  const [name, setName] = useState("");
  const [tagIds, setTagIds] = useState<number[]>([]);
  const [initialTagIds, setInitialTagIds] = useState<number[]>([]);
  const [isSaving, setIsSaving] = useState(false);
  const [initialized, setInitialized] = useState(false);

  useEffect(() => {
    if (project && !initialized) {
      setName(project.name);
      const ids = (project.tags as ProjectTag[]).map((t) => t.id);
      setTagIds(ids);
      setInitialTagIds(ids);
      setInitialized(true);
    }
  }, [project, initialized]);

  const dir = useProjectDirectory(name, "move");

  const nameError = getProjectNameError(name, projects, project?.id);

  async function saveProject() {
    if (!project) return;
    setIsSaving(true);
    try {
      if (name !== project.name) {
        await api.patch(`projects/${project.id}/`, { name });
      }

      const added = tagIds.filter((tid) => !initialTagIds.includes(tid));
      const removed = initialTagIds.filter((tid) => !tagIds.includes(tid));
      for (const tagId of added) {
        await apiPost(`projects/${project.id}/tags/`, { tag_id: tagId });
      }
      for (const tagId of removed) {
        await apiDelete(`projects/${project.id}/tags/${tagId}/`);
      }

      if (dir.effectiveDirectory) {
        await apiPost(`projects/${project.id}/move_directory/`, {
          new_directory: dir.effectiveDirectory,
        });
      }

      router.push(`/ccp4i2/project/${project.id}`);
    } catch (err) {
      console.error("Error saving project:", err);
      alert("Error saving project: " + err);
      setIsSaving(false);
    }
  }

  if (isLoading || !initialized) {
    return (
      <Container maxWidth="sm" sx={{ my: 3 }}>
        <Stack spacing={2}>
          <Skeleton variant="text" height={48} />
          <Skeleton variant="rounded" height={56} />
          <Skeleton variant="rounded" height={56} />
        </Stack>
      </Container>
    );
  }

  return (
    <Container maxWidth="sm" sx={{ my: 3 }}>
      <Stack spacing={2}>
        <Typography variant="h4">Edit Project</Typography>
        <ProjectNameField
          value={name}
          onChange={setName}
          onSubmit={saveProject}
          error={nameError}
          autoFocus
        />
        <ProjectDirectoryField
          mode="move"
          currentDirectory={project!.directory}
          {...dir}
        />
        <EditTags tags={tagIds} onChange={setTagIds} />
        <Stack direction="row" spacing={2} justifyContent="flex-end">
          <Button
            variant="outlined"
            onClick={() => router.push(`/ccp4i2/project/${project!.id}`)}
            disabled={isSaving}
          >
            Cancel
          </Button>
          <Button
            variant="contained"
            onClick={saveProject}
            disabled={
              nameError.length > 0 || dir.directoryHasError || isSaving
            }
            startIcon={isSaving ? <CircularProgress size={16} /> : undefined}
          >
            {isSaving ? "Saving..." : "Save"}
          </Button>
        </Stack>
      </Stack>
    </Container>
  );
}

