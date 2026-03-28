/*
 * Copyright (C) 2026 Newcastle University
 *
 * This file is part of CCP4i2.
 *
 * CCP4i2 is free software: you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License version 3,
 * modified in accordance with the provisions of the license to address
 * the requirements of UK law.
 *
 * See https://www.ccp4.ac.uk/ccp4license.php for details.
 */
/**
 * CFreeRDataFileElement - FreeR data file selector with inline generation.
 *
 * Wraps CMiniMtzDataFileElement (for MTZ column selection and upload)
 * and adds a "Generate Free-R Flags" button that opens an inline modal
 * to run the freerflag task without navigating away.
 */
import {
  Box,
  Button,
  Stack,
  Tooltip,
  Typography,
} from "@mui/material";
import { Add as AddIcon } from "@mui/icons-material";
import { CMiniMtzDataFileElement } from "./cminimtzdatafile";
import { CCP4i2TaskElementProps } from "./task-element";
import { useCallback, useMemo, useState } from "react";
import { useApi } from "../../../api";
import { useJob, useProject, useProjectFiles } from "../../../utils";
import { File as DjangoFile, Project } from "../../../types/models";
import { InlineTaskModal } from "./inline-task-modal";

export const CFreeRDataFileElement: React.FC<CCP4i2TaskElementProps> = (
  props
) => {
  const { job, itemName, qualifiers } = props;
  const api = useApi();
  const { useTaskItem, setParameter, fileItemToParameterArg, mutateContainer } = useJob(job.id);
  const { item, value } = useTaskItem(itemName);
  const { files: projectFiles } = useProjectFiles(job.project);
  const { jobs: projectJobs } = useProject(job.project);
  const { data: projects } = api.get<Project[]>("projects");

  // State for inline task modal
  const [modalOpen, setModalOpen] = useState(false);

  // Check if a file is already selected
  const hasFile = Boolean(value?.dbFileId);

  // Check if there are any existing FreeR files in the project
  const existingFreeRFiles = useMemo(() => {
    if (!projectFiles) return [];
    return projectFiles.filter((file) => file.type === "application/CCP4-mtz-freerflag");
  }, [projectFiles]);

  /**
   * Handle the output file from the inline freerflag task.
   * Auto-selects it in this widget via setParameter.
   */
  const handleOutputFileReady = useCallback(
    async (outputFile: DjangoFile) => {
      if (!item?._objectPath || !projectJobs) return;

      try {
        const paramArg = fileItemToParameterArg(
          outputFile,
          item._objectPath,
          projectJobs,
          projects || []
        );

        const result = await setParameter(paramArg);
        if (result?.success) {
          props.onChange?.(result.data?.updated_item);
        }

        await mutateContainer();
      } catch (error) {
        console.error("Error auto-selecting FreeR output file:", error);
      }
    },
    [item?._objectPath, projectJobs, projects, fileItemToParameterArg, setParameter, props.onChange, mutateContainer]
  );

  return (
    <>
      <CMiniMtzDataFileElement
        {...props}
        qualifiers={{ ...qualifiers, downloadModes: ["ebiSFs"] }}
      >
        {/* Generate Free-R Flags button - shown when no file selected and job is editable */}
        {!hasFile && job.status === 1 && (
          <Box sx={{ mt: 1 }}>
            <Stack direction="row" spacing={1} alignItems="center">
              <Tooltip title="Generate Free-R flags by running the freerflag task">
                <Button
                  variant="outlined"
                  size="small"
                  startIcon={<AddIcon />}
                  onClick={() => setModalOpen(true)}
                  sx={{ textTransform: "none" }}
                >
                  Generate Free-R Flags
                </Button>
              </Tooltip>
              {existingFreeRFiles.length === 0 && (
                <Typography variant="caption" color="text.secondary">
                  No Free-R files in project
                </Typography>
              )}
            </Stack>
          </Box>
        )}
      </CMiniMtzDataFileElement>

      {/* Inline modal for running freerflag */}
      <InlineTaskModal
        open={modalOpen}
        onClose={() => setModalOpen(false)}
        taskName="freerflag"
        parentJob={job}
        onOutputFileReady={handleOutputFileReady}
        title="Generate Free-R Flags"
      />
    </>
  );
};
