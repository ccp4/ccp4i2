"use client";
import { ChangeEvent, useEffect, useMemo, useState } from "react";
import path from "path";
import { useRouter } from "next/navigation";
import {
  Button,
  Container,
  Stack,
  TextField,
  ToggleButton,
  ToggleButtonGroup,
  Typography,
} from "@mui/material";
import { Folder } from "@mui/icons-material";
import { useApi } from "../api";
import { Project } from "../types/models";
import EditTags from "./edit-tags";

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

  function createProject() {
    const formData = new FormData();
    formData.append("name", name);
    formData.append("directory", customDirectory ? directory : "__default__");
    try {
      api.post<Project>("projects", formData).then((project) => {
        router.push(`/project/${project.id}`);
      });
    } catch (err) {
      console.error("Error creating project:", err);
      alert("Error creating project: " + err);
    }
  }

  function handleNameChange(event: ChangeEvent<HTMLInputElement>) {
    setName(event.target.value);
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
          required
          error={nameError.length > 0}
          helperText={nameError}
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
        <Stack direction="row" spacing={2} justifyContent="flex-end">
          <Button variant="outlined" onClick={() => router.push("/")}>
            Cancel
          </Button>
          <Button
            variant="contained"
            onClick={createProject}
            disabled={nameError.length > 0 || directoryError.length > 0}
          >
            Create
          </Button>
        </Stack>
      </Stack>
    </Container>
  );
};
