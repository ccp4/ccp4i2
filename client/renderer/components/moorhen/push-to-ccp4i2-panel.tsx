import { useApi } from "../../api";
import { Project as ProjectInfo, Job as JobInfo } from "../../types/models";
import { moorhen } from "moorhen/types/moorhen";
import {
  Autocomplete,
  TextField,
  Box,
  Typography,
  Button,
} from "@mui/material";
import React, { useCallback, useState, useEffect } from "react";
import { CreateTaskResponse } from "../../utils";
import { usePopcorn } from "../../providers/popcorn-provider";
import { ItemMetadata } from "./item-metadata-utils";

interface PushToCCP4i2Props {
  project?: ProjectInfo;
  molNo?: number;
  item?: moorhen.Molecule | moorhen.Map | null;
  itemMetadata?: ItemMetadata;
  onClose: () => void;
}

type ProjectInfoOrNull = ProjectInfo | null;
function detectCoordinateFormat(text: string): "pdb" | "mmcif" | "unknown" {
  const trimmedText = text.replace(/^\s+/, ""); // Remove leading whitespace and blank lines
  if (/^(HEADER|TITLE|ATOM  |HETATM)/m.test(trimmedText)) {
    return "pdb";
  }
  if (/^data_/m.test(trimmedText) && /_atom_site\./.test(trimmedText)) {
    return "mmcif";
  }
  return "unknown";
}

export const PushToCCP4i2Panel: React.FC<PushToCCP4i2Props> = ({
  project,
  molNo,
  item,
  itemMetadata,
  onClose,
}) => {
  const api = useApi();
  const { data: projects } = api.get<ProjectInfo[]>("projects") || [];
  const [selectedProject, setSelectedProject] = useState<
    ProjectInfo | undefined
  >(project);
  const { setMessage } = usePopcorn();

  // Set initial project based on itemMetadata.projectName if available
  useEffect(() => {
    if (
      itemMetadata &&
      itemMetadata.projectName &&
      Array.isArray(projects) &&
      projects.length > 0
    ) {
      const matchedProject = projects.find(
        (proj) => proj.name === itemMetadata.projectName
      );
      if (matchedProject && selectedProject?.id !== matchedProject.id) {
        setSelectedProject(matchedProject);
      }
    }
  }, [itemMetadata, projects]);

  const { data: jobs, mutate: mutateJobs } = api.get<JobInfo[]>(
    `projects/${selectedProject?.id}/jobs/`
  );

  const handlePushToCCP4i2 = useCallback(async () => {
    console.log({ molNo, item });
    // Implement your push logic here
    if (selectedProject && item) {
      try {
        setMessage("Pushing model coordinates to CCP4i2...", "info");
        console.log("Pushing to CCP4i2:", selectedProject);
        const result = await api.post<CreateTaskResponse>(
          `projects/${selectedProject.id}/create_task/`,
          {
            task_name: "coordinate_selector",
          }
        );

        if ((result as any)?.success === false) {
          setMessage(`Failed to create task: ${(result as any)?.error || 'Unknown error'}`, "error");
          return;
        }

        mutateJobs();
        const modelCoords =
          item.type === "molecule"
            ? await (item as moorhen.Molecule).getAtoms()
            : null;
        if (!modelCoords) return;

        const format = detectCoordinateFormat(modelCoords);
        setMessage(`Detected coordinate format: ${format}`, "info");

        const slugify = (name: string) =>
          name
            .replace(/[/\\?%*:|"<>]/g, "") // Remove illegal filename chars
            .replace(/\s+/g, "_") // Replace whitespace with underscores
            .replace(/[^a-zA-Z0-9._-]/g, "") // Remove other non-safe chars
            .replace(/^_+|_+$/g, ""); // Trim leading/trailing underscores

        const moleculeName =
          slugify((item as moorhen.Molecule).name) +
          (format === "mmcif" ? ".cif" : ".pdb");

        const formData = new FormData();
        formData.append("objectPath", "coordinate_selector.inputData.XYZIN");
        formData.append(
          "file",
          new Blob([modelCoords], {
            type: format === "mmcif" ? "chemical/x-cif" : "chemical/x-pdb",
          }),
          moleculeName
        );
        const newJobId = result.data?.new_job?.id;
        if (!newJobId) {
          setMessage(`Failed to create job: ${(result as any)?.error || 'Unknown error'}`, "error");
          return;
        }
        // Note: Direct API call is OK here because:
        // 1. This uploads to a NEWLY created job (not an existing one with SWR cache)
        // 2. No concurrent SWR fetching exists for this job yet
        // 3. The job runs immediately after upload, so no UI race conditions
        // For existing jobs, use uploadFileParam from useJob() for proper intent tracking
        const uploadResult = await api.post<any>(
          `jobs/${newJobId}/upload_file_param`,
          formData
        );

        if (uploadResult?.success === false) {
          setMessage(`Upload failed: ${uploadResult?.error || 'Unknown error'}`, "error");
          return;
        }

        const run_result = await api.post<CreateTaskResponse>(
          //Call run_local for more responsive task execution of this (which should be faster )
          `jobs/${newJobId}/run_local/`,
          {
            task_name: "coordinate_selector",
          }
        );

        if ((run_result as any)?.success === false) {
          setMessage(`Run failed: ${(run_result as any)?.error || 'Unknown error'}`, "error");
          return;
        }
        setMessage(`Model pushed to CCP4i2 successfully`, "success");
        if (onClose) onClose();
      } catch (error) {
        setMessage(`Error pushing to CCP4i2: ${error instanceof Error ? error.message : String(error)}`, "error");
      }
    }
  }, [selectedProject, molNo, item, api, mutateJobs, onClose, setMessage]);

  return (
    <Box sx={{ p: 2 }}>
      <Typography variant="h6" gutterBottom>
        Push to CCP4i2
      </Typography>

      {/* Show nicely formatted metadata if provided */}
      {itemMetadata && (
        <Box
          sx={{
            mb: 2,
            p: 2,
            border: "1px solid #eee",
            borderRadius: 2,
            background: "#fafafa",
          }}
        >
          <Typography variant="subtitle1" gutterBottom>
            File to push was fetched from CCP4i2 with the following metadata:
          </Typography>
          <Typography variant="body2">
            <strong>File ID:</strong> {itemMetadata.fileId}
          </Typography>
          {itemMetadata.projectName && (
            <Typography variant="body2">
              <strong>Project:</strong> {itemMetadata.projectName}
            </Typography>
          )}
          {itemMetadata.jobNumber && (
            <Typography variant="body2">
              <strong>Job Number:</strong> {itemMetadata.jobNumber}
            </Typography>
          )}
          {itemMetadata.fileAnnotation && (
            <Typography variant="body2">
              <strong>Annotation:</strong> {itemMetadata.fileAnnotation}
            </Typography>
          )}
          {itemMetadata.isLoading && (
            <Typography variant="body2" color="text.secondary">
              Loading metadata...
            </Typography>
          )}
          {itemMetadata.error && (
            <Typography variant="body2" color="error">
              Error: {itemMetadata.error}
            </Typography>
          )}
        </Box>
      )}

      <Autocomplete
        options={projects || []}
        getOptionLabel={(option) => option.name}
        //@ts-ignore
        value={selectedProject || null}
        onChange={(_, value) => setSelectedProject(value)}
        renderInput={(params) => (
          <TextField
            {...params}
            label="Select Project"
            variant="outlined"
            fullWidth
          />
        )}
        isOptionEqualToValue={(option, value) => option.id === value?.id}
        disableClearable
      />
      <Button
        sx={{ mt: 2 }}
        variant="contained"
        color="primary"
        disabled={!selectedProject}
        onClick={handlePushToCCP4i2}
      >
        Push to CCP4i2
      </Button>
      {/* Additional panel content can go here */}
    </Box>
  );
};
