"use client";
import { useState } from "react";
import { useRouter } from "next/navigation";
import {
  Button,
  CircularProgress,
  Container,
  Stack,
  Typography,
} from "@mui/material";
import { useApi } from "../api";
import { apiPost } from "../api-fetch";
import { Project } from "../types/models";
import EditTags from "./edit-tags";
import { ProjectNameField, getProjectNameError } from "./project-name-field";
import { ProjectDirectoryField } from "./project-directory-field";
import { useProjectDirectory } from "../hooks/use-project-directory";
import {
  DroppedFile,
  FileDropZone,
  TASK_FOR_TYPE,
  PARAM_FOR_TYPE,
  AUTO_RUN_FOR_TYPE,
} from "./file-drop-zone";

export const NewProjectContent: React.FC = () => {
  const api = useApi();
  const router = useRouter();
  const [name, setName] = useState("");
  const [tags, setTags] = useState<number[]>([]);
  const [droppedFiles, setDroppedFiles] = useState<DroppedFile[]>([]);
  const [isCreating, setIsCreating] = useState(false);
  const { data: projects } = api.get<Project[]>("projects");

  const dir = useProjectDirectory(name, "create");

  async function createProject() {
    setIsCreating(true);
    try {
      const formData = new FormData();
      formData.append("name", name);
      formData.append("directory", dir.effectiveDirectory ?? "__default__");
      const project = await api.post<Project>("projects", formData);

      // Apply tags to the new project
      for (const tagId of tags) {
        try {
          await apiPost(`projects/${project.id}/tags/`, { tag_id: tagId });
        } catch (err) {
          console.error(`Error applying tag ${tagId}:`, err);
        }
      }

      // If files were dropped, create import jobs sequentially
      const importableFiles = droppedFiles.filter(
        (df) => TASK_FOR_TYPE[df.detectedType] !== null
      );

      for (const df of importableFiles) {
        try {
          await createImportJob(project.id, df);
          // Small delay between jobs to avoid DB contention (SQLite)
          if (importableFiles.length > 1) {
            await new Promise((r) => setTimeout(r, 500));
          }
        } catch (err) {
          console.error(`Error importing ${df.file.name}:`, err);
        }
      }

      router.push(`/ccp4i2/project/${project.id}`);
    } catch (err) {
      console.error("Error creating project:", err);
      alert("Error creating project: " + err);
      setIsCreating(false);
    }
  }

  async function createImportJob(projectId: number, df: DroppedFile) {
    const taskName = TASK_FOR_TYPE[df.detectedType];
    const paramPath = PARAM_FOR_TYPE[df.detectedType];
    if (!taskName || !paramPath) return;

    // 1. Create the job
    const result = await apiPost<any>(`projects/${projectId}/create_task/`, {
      task_name: taskName,
      title: df.detectedType === "unmerged"
        ? `Merge ${df.file.name}`
        : `Import ${df.file.name}`,
    });
    if (!result?.success || !result.data?.new_job) {
      throw new Error(`Failed to create ${taskName} job`);
    }
    const jobId = result.data.new_job.id;

    // 2. Upload the file
    if (df.detectedType === "sequence") {
      // Sequences: read text content and set as parameter
      const text = await df.file.text();
      await apiPost(`jobs/${jobId}/set_parameter/`, {
        object_path: `${taskName}.container.${paramPath}`,
        value: text,
      });
    } else if (df.detectedType === "unmerged") {
      // Unmerged data: add list item, then upload to the list slot
      await apiPost(`jobs/${jobId}/set_parameter/`, {
        object_path: `${taskName}.container.${paramPath}`,
        value: [{}],
      });
      const uploadForm = new FormData();
      uploadForm.append("file", df.file, df.file.name);
      uploadForm.append(
        "objectPath",
        `${taskName}.container.${paramPath}[0].file`
      );
      await apiPost(`jobs/${jobId}/upload_file_param/`, uploadForm);
    } else if (df.detectedType === "ligand") {
      // Ligands: upload file and set the mode selector
      const ext = df.file.name.toLowerCase().split(".").pop();
      const mode = ext === "mol2" ? "MOL2" : "MOL";
      const actualParam =
        ext === "mol2" ? "inputData.MOL2IN" : "inputData.MOLIN";
      await apiPost(`jobs/${jobId}/set_parameter/`, {
        object_path: `${taskName}.container.inputData.MOLSMILESORSKETCH`,
        value: mode,
      });
      const uploadForm = new FormData();
      uploadForm.append("file", df.file, df.file.name);
      uploadForm.append("object_path", `${taskName}.container.${actualParam}`);
      await apiPost(`jobs/${jobId}/upload_file_param/`, uploadForm);
    } else {
      // Standard file upload
      const uploadForm = new FormData();
      uploadForm.append("file", df.file, df.file.name);
      uploadForm.append("object_path", `${taskName}.container.${paramPath}`);
      await apiPost(`jobs/${jobId}/upload_file_param/`, uploadForm);
    }

    // 3. Only auto-run simple import jobs; leave data reduction jobs
    //    for the user to review parameters before running
    if (AUTO_RUN_FOR_TYPE[df.detectedType]) {
      await apiPost(`jobs/${jobId}/run/`, {});
    }
  }

  const nameError = getProjectNameError(name, projects);

  return (
    <Container maxWidth="sm" sx={{ my: 3 }}>
      <Stack spacing={2}>
        <Typography variant="h4">Create Project</Typography>
        <ProjectNameField
          value={name}
          onChange={setName}
          onSubmit={createProject}
          error={nameError}
          autoFocus
        />
        <ProjectDirectoryField mode="create" {...dir} />

        <EditTags tags={tags} onChange={setTags} />

        <FileDropZone files={droppedFiles} onChange={setDroppedFiles} />

        <Stack direction="row" spacing={2} justifyContent="flex-end">
          <Button
            variant="outlined"
            onClick={() => router.push("/ccp4i2")}
            disabled={isCreating}
          >
            Cancel
          </Button>
          <Button
            variant="contained"
            onClick={createProject}
            disabled={
              nameError.length > 0 || dir.directoryHasError || isCreating
            }
            startIcon={isCreating ? <CircularProgress size={16} /> : undefined}
          >
            {isCreating
              ? droppedFiles.length > 0
                ? "Creating & importing..."
                : "Creating..."
              : droppedFiles.length > 0
                ? `Create & import ${droppedFiles.length} file${droppedFiles.length > 1 ? "s" : ""}`
                : "Create"}
          </Button>
        </Stack>
      </Stack>
    </Container>
  );
};
