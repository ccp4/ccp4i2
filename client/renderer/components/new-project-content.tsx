"use client";
import { ChangeEvent, useEffect, useMemo, useState } from "react";
import path from "path";
import { useRouter } from "next/navigation";
import {
  Button,
  CircularProgress,
  Container,
  Stack,
  TextField,
  ToggleButton,
  ToggleButtonGroup,
  Typography,
} from "@mui/material";
import { Folder } from "@mui/icons-material";
import { useApi } from "../api";
import { apiPost } from "../api-fetch";
import { Project } from "../types/models";
import EditTags from "./edit-tags";
import {
  DroppedFile,
  FileDropZone,
  TASK_FOR_TYPE,
  PARAM_FOR_TYPE,
} from "./file-drop-zone";

export const NewProjectContent: React.FC = () => {
  const api = useApi();
  const router = useRouter();
  const [name, setName] = useState("");
  const [customDirectory, setCustomDirectory] = useState(false);
  const [ccp4i2ProjectDirectory, setCcp4i2ProjectDirectory] = useState<string>(
    "/home/user/CCP4X_PROJECTS"
  );
  const [directoryExists, setDirectoryExists] = useState(true);
  const [electronAPIAvailable, setElectronAPIAvailable] =
    useState<boolean>(false);
  const [tags, setTags] = useState<number[]>([]);
  const [droppedFiles, setDroppedFiles] = useState<DroppedFile[]>([]);
  const [isCreating, setIsCreating] = useState(false);
  const { data: projects } = api.get<Project[]>("projects");

  useEffect(() => {
    // Send a message to the main process to get the config
    if (window.electronAPI) {
      setElectronAPIAvailable(true);
      window.electronAPI.sendMessage("get-config");
      // Listen for messages from the main process
      window.electronAPI.onMessage(
        "message-from-main",
        (event: any, data: any) => {
          if (data.message === "get-config") {
            setCcp4i2ProjectDirectory(data.config.CCP4I2_PROJECTS_DIR);
          }
          if (data.message === "check-file-exists") {
            setDirectoryExists(data.exists);
          }
        }
      );
    } else {
      console.log("window.electron is not available");
      setDirectoryExists(false); // Assume directory does not exist in web mode
    }
  }, []);

  const directory = useMemo(() => {
    const result = path.join(
      ccp4i2ProjectDirectory || "",
      name.toLocaleLowerCase()
    );
    if (typeof window !== "undefined" && window.electronAPI) {
      window.electronAPI.sendMessage("check-file-exists", { path: result });
    }
    return result;
  }, [ccp4i2ProjectDirectory, name]);

  async function createProject() {
    setIsCreating(true);
    try {
      const formData = new FormData();
      formData.append("name", name);
      formData.append("directory", customDirectory ? directory : "__default__");
      const project = await api.post<Project>("projects", formData);

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
      title: `Import ${df.file.name}`,
    });
    if (!result?.success || !result.data?.new_job) {
      throw new Error(`Failed to create ${taskName} job`);
    }
    const jobId = result.data.new_job.id;

    // 2. Upload the file
    // For sequences, we read the text content and set it as a parameter
    if (df.detectedType === "sequence") {
      const text = await df.file.text();
      await apiPost(`jobs/${jobId}/set_parameter/`, {
        object_path: `${taskName}.container.${paramPath}`,
        value: text,
      });
    } else {
      const uploadForm = new FormData();
      uploadForm.append("file", df.file, df.file.name);
      uploadForm.append("objectPath", `${taskName}.container.${paramPath}`);
      await apiPost(`jobs/${jobId}/upload_file_param/`, uploadForm);
    }

    // 3. Run the job (fire-and-forget)
    await apiPost(`jobs/${jobId}/run/`, {});
  }

  function handleNameChange(event: ChangeEvent<HTMLInputElement>) {
    setName(event.target.value);
  }

  function handleNameKeyDown(event: React.KeyboardEvent<HTMLInputElement>) {
    if (event.key === "Enter" && nameError.length === 0 && directoryError.length === 0) {
      createProject();
    }
  }

  function handleCustomDirectoryChange(
    event: React.MouseEvent<HTMLElement>,
    value: any
  ) {
    if (value !== null) {
      setCustomDirectory(value);
    }
  }

  function handleDirectoryChange() {
    if (typeof window !== "undefined") {
      if (window.electronAPI) {
        window.electronAPI.sendMessage("locate-ccp4i2-project-directory");
      } else {
        console.error("Electron API is not available");
      }
    }
  }

  let nameError = "";
  if (name.length === 0) nameError = "Name is required";
  else if (!name.match("^[A-z0-9_-]+$"))
    nameError =
      "Name can only contain letters, numbers, underscores, and hyphens";
  else if (projects?.find((p) => p.name === name))
    nameError = "Name is already taken";

  const directoryError = useMemo(() => {
    if (customDirectory && directory.length === 0)
      return "Directory is required";
    else if (directory.length > 0 && directoryExists)
      return "Directory already exist";
    return "";
  }, [directoryExists, customDirectory, directory]);

  return (
    <Container maxWidth="sm" sx={{ my: 3 }}>
      <Stack spacing={2}>
        <Typography variant="h4">Create Project</Typography>
        <TextField
          label="Name"
          value={name}
          onChange={handleNameChange}
          onKeyDown={handleNameKeyDown}
          required
          error={nameError.length > 0}
          helperText={nameError}
          autoFocus
        />
        <ToggleButtonGroup
          exclusive
          value={customDirectory}
          onChange={handleCustomDirectoryChange}
          fullWidth
        >
          <ToggleButton value={false}>Default Directory</ToggleButton>
          <ToggleButton value={true}>Custom Directory</ToggleButton>
        </ToggleButtonGroup>
        {customDirectory && electronAPIAvailable && (
          <Stack direction="row" spacing={2} sx={{ alignItems: "center" }}>
            <TextField
              label="Parent directory where the project directory will be created"
              value={ccp4i2ProjectDirectory}
              disabled={true}
              sx={{ flexGrow: 1 }}
              required
            />
            <Button
              variant="outlined"
              startIcon={<Folder />}
              onClick={handleDirectoryChange}
            >
              Select
            </Button>
          </Stack>
        )}
        <Stack direction="row">
          <TextField
            label="Resulting name for project directory"
            value={directory}
            disabled={true}
            error={directoryError.length > 0}
            helperText={directoryError}
            sx={{ flexGrow: 1 }}
          />
        </Stack>

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
              nameError.length > 0 || directoryError.length > 0 || isCreating
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
