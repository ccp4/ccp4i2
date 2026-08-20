"use client";
import { useEffect, useMemo, useState } from "react";
import {
  Button,
  CircularProgress,
  Container,
  Divider,
  LinearProgress,
  Stack,
  TextField,
  Typography,
} from "@mui/material";
import { DriveFileMove } from "@mui/icons-material";
import { useApi } from "../api";
import { Project } from "../types/models";
import { usePopcorn } from "../providers/popcorn-provider";
import { isElectron } from "../utils/platform";
import { MoveProjectDialog } from "./move-project-dialog";
import { TagsOfProject } from "./tags-of-project";

/**
 * View and change what a project is, as opposed to what is in it: its name,
 * its description, its tags and where it lives on disk.
 *
 * The name and the description are a form, saved together. Tags and the
 * directory are not: tags are applied as they are picked, and moving a project
 * rewrites files as well as rows, so it keeps its own dialog and its own
 * confirmation rather than being one unsaved edit among several.
 */
export const EditProjectContent: React.FC<{ projectId: number }> = ({
  projectId,
}) => {
  const api = useApi();
  const { setMessage } = usePopcorn();
  const { data: project, mutate: mutateProject } = api.get<Project>(
    `projects/${projectId}`
  );
  const { data: projects } = api.get<Project[]>("projects");

  const [name, setName] = useState("");
  const [description, setDescription] = useState("");
  const [isSaving, setIsSaving] = useState(false);
  const [moveOpen, setMoveOpen] = useState(false);
  // Only the desktop build has a filesystem to move a project around in.
  const [canMove, setCanMove] = useState(false);
  useEffect(() => setCanMove(isElectron()), []);

  // Seed the fields from the project once. Later revalidations -- a tag added,
  // a move finished -- must not overwrite what is being typed.
  const [seeded, setSeeded] = useState(false);
  useEffect(() => {
    if (project && !seeded) {
      setName(project.name);
      setDescription(project.description ?? "");
      setSeeded(true);
    }
  }, [project, seeded]);

  const nameError = useMemo(() => {
    if (name.length === 0) return "Name is required";
    if (!name.match("^[A-z0-9_-]+$"))
      return "Name can only contain letters, numbers, underscores, and hyphens";
    if (projects?.find((p) => p.name === name && p.id !== projectId))
      return "Name is already taken";
    return "";
  }, [name, projects, projectId]);

  const isDirty =
    project !== undefined &&
    (name !== project.name || description !== (project.description ?? ""));

  function discard() {
    if (!project) return;
    setName(project.name);
    setDescription(project.description ?? "");
  }

  async function save() {
    if (!project) return;
    setIsSaving(true);
    try {
      await api.patch(`projects/${project.id}`, { name, description });
      await mutateProject();
      setMessage("Project updated", "success");
    } catch (err) {
      setMessage(
        `Error updating project: ${err instanceof Error ? err.message : String(err)}`,
        "error"
      );
    } finally {
      setIsSaving(false);
    }
  }

  if (!project) return <LinearProgress />;

  return (
    <Container maxWidth="sm" sx={{ my: 3 }}>
      <Stack spacing={2}>
        <Typography variant="h4">Edit Project</Typography>

        <TextField
          label="Name"
          value={name}
          onChange={(event) => setName(event.target.value)}
          required
          error={nameError.length > 0}
          helperText={nameError}
          autoFocus
        />
        <TextField
          label="Description"
          value={description}
          onChange={(event) => setDescription(event.target.value)}
          multiline
          minRows={3}
          helperText="What the project is for; shown wherever it is listed."
        />
        <Stack direction="row" spacing={2} justifyContent="flex-end">
          <Button variant="outlined" onClick={discard} disabled={!isDirty || isSaving}>
            Discard changes
          </Button>
          <Button
            variant="contained"
            onClick={save}
            disabled={nameError.length > 0 || !isDirty || isSaving}
            startIcon={isSaving ? <CircularProgress size={16} /> : undefined}
          >
            {isSaving ? "Saving..." : "Save changes"}
          </Button>
        </Stack>

        <Divider />

        <Typography variant="h6">Tags</Typography>
        <TagsOfProject projectId={project.id} />

        <Divider />

        <Typography variant="h6">Location</Typography>
        <TextField
          label="Project directory"
          value={project.directory}
          disabled
          helperText={
            canMove
              ? "Renaming a project does not move it; use Move to relocate it on disk."
              : "Moving a project on disk is available in the desktop app."
          }
        />
        {canMove && (
          <Stack direction="row" justifyContent="flex-end">
            <Button
              variant="outlined"
              startIcon={<DriveFileMove />}
              onClick={() => setMoveOpen(true)}
            >
              Move...
            </Button>
          </Stack>
        )}

        <MoveProjectDialog
          project={moveOpen ? project : null}
          open={moveOpen}
          onClose={() => setMoveOpen(false)}
          onMoved={() => mutateProject()}
        />
      </Stack>
    </Container>
  );
};
