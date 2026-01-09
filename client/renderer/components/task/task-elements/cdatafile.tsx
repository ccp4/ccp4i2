import React, {
  ChangeEvent,
  PropsWithChildren,
  ReactNode,
  SyntheticEvent,
  useCallback,
  useEffect,
  useMemo,
  useState,
} from "react";
import {
  Autocomplete,
  AutocompleteChangeReason,
  Avatar,
  Box,
  IconButton,
  LinearProgress,
  Stack,
  TextField,
  Tooltip,
} from "@mui/material";
import {
  Menu as MenuIcon,
  ChevronRight as ChevronRightIcon,
  Clear,
} from "@mui/icons-material";
import { useDndContext, useDroppable } from "@dnd-kit/core";

import { useApi } from "../../../api";
import { useJob, useProject, useProjectFiles } from "../../../utils";
import { CCP4i2TaskElementProps } from "./task-element";
import { File as CCP4i2File, nullFile, Project } from "../../../types/models";
import { useTaskInterface } from "../../../providers/task-provider";
import { useFileMenu } from "../../../providers/file-context-menu";
import { ErrorTrigger } from "./error-info";
import { InputFileFetch } from "./input-file-fetch";
import { InputFileUpload } from "./input-file-upload";
import { FIELD_SPACING } from "./field-sizes";
import { ExpandableSection } from "./expandable-section";

// Types
export interface CCP4i2DataFileElementProps
  extends CCP4i2TaskElementProps,
    PropsWithChildren {
  setFileContent?: (fileContent: ArrayBuffer | string | File | null) => void;
  setFiles?: (files: FileList | null) => void;
  infoContent?: ReactNode;
  onChange?: (updatedItem: any) => void;
  hasValidationError?: boolean;
  forceExpanded?: boolean;
}

/** Get a friendly file type label */
const getFileTypeLabel = (className: string | undefined): string => {
  if (!className) return "File";
  // Remove leading 'C' and trailing 'DataFile' or 'File'
  const clean = className
    .replace(/^C/, "")
    .replace(/DataFile$/, "")
    .replace(/File$/, "");
  return clean || "File";
};

// Main component
export const CDataFileElement: React.FC<CCP4i2DataFileElementProps> = ({
  job,
  sx,
  itemName,
  onChange,
  setFiles,
  children,
  visibility,
  disabled: disabledProp,
  qualifiers: propsQualifiers,
  hasValidationError: overrideValidationError,
  forceExpanded = false,
}) => {
  const api = useApi();
  const {
    useTaskItem,
    setParameter,
    getValidationColor,
    fileItemToParameterArg,
    mutateContainer,
    useFileDigest,
    useFileContent,
  } = useJob(job.id);

  const { item } = useTaskItem(itemName);
  const { inFlight, setInFlight } = useTaskInterface();
  const { setFileMenuAnchorEl, setFile } = useFileMenu();

  // Data and state
  // Use polling for files so task widgets see newly created output files
  const { files: projectFiles } = useProjectFiles(job.project, true);
  const { jobs: projectJobs } = useProject(job.project);
  const { data: projects } = api.get<Project[]>("projects");
  const { mutate: mutateDigest } = useFileDigest(`${item?._objectPath}`);
  const { mutate: mutateContent } = useFileContent(`${item?._objectPath}`);
  const [value, setValue] = useState<CCP4i2File>(nullFile);
  const [isManuallyExpanded, setIsManuallyExpanded] = useState(false);

  // Computed values
  const qualifiers = useMemo(
    () => ({
      ...item?._qualifiers,
      ...propsQualifiers,
    }),
    [item?._qualifiers, propsQualifiers]
  );

  const fileConfig = useMemo(() => {
    const allowedTypes = qualifiers?.mimeTypeName
      ? Array.isArray(qualifiers.mimeTypeName)
        ? qualifiers.mimeTypeName
        : [qualifiers.mimeTypeName]
      : null;
    const acceptedExtensions =
      qualifiers?.fileExtensions?.map((ext: string) => `.${ext}`).join(",") ||
      "";
    // requiredContentFlag filters files by their content flag (e.g., [1,2] for anomalous pairs)
    // This is critical for tasks like SAD/MAD phasing that need specific data types
    const requiredContentFlag = qualifiers?.requiredContentFlag
      ? Array.isArray(qualifiers.requiredContentFlag)
        ? qualifiers.requiredContentFlag
        : [qualifiers.requiredContentFlag]
      : null;
    return { allowedTypes, acceptedExtensions, requiredContentFlag };
  }, [qualifiers]);

  const fileOptions = useMemo(() => {
    if (!projectFiles || !fileConfig.allowedTypes) return [];
    return projectFiles
      .filter((file) => {
        const fileJob = projectJobs?.find((job) => job.id === file.job);
        const isValidType =
          fileConfig.allowedTypes!.includes(file.type) ||
          fileConfig.allowedTypes!.includes("Unknown");
        const isNotParentJob = fileJob ? !fileJob.parent : true;
        // Filter by requiredContentFlag if specified (and non-empty)
        // This prevents selecting incompatible files (e.g., IMEAN for SAD phasing)
        // null, undefined, or empty array means no filtering
        const hasValidContentFlag =
          !fileConfig.requiredContentFlag ||
          fileConfig.requiredContentFlag.length === 0 ||
          fileConfig.requiredContentFlag.includes(file.content);
        return isValidType && isNotParentJob && hasValidContentFlag;
      })
      .sort((a, b) => b.job - a.job);
  }, [projectFiles, projectJobs, fileConfig.allowedTypes, fileConfig.requiredContentFlag]);

  const borderColor = getValidationColor(item);
  const hasError = borderColor === "error.light";
  const computedValidationError = hasError;
  const hasValidationError = overrideValidationError ?? computedValidationError;
  const hasChildren = React.Children.count(children) > 0;
  const isExpanded = hasValidationError || isManuallyExpanded || forceExpanded;

  const guiLabel =
    qualifiers?.guiLabel || item?._objectPath?.split(".").at(-1) || "";
  const isDisabled =
    (typeof disabledProp === "function" ? disabledProp() : disabledProp) ||
    inFlight ||
    job.status !== 1;
  const isVisible =
    typeof visibility === "function" ? visibility() : visibility !== false;

  // Drag and drop
  const { isOver, setNodeRef } = useDroppable({
    id: `job_${job.uuid}_${itemName}`,
    data: { job, item },
  });
  const { active } = useDndContext();
  const isValidDrop =
    active?.data?.current?.file &&
    (fileConfig.allowedTypes?.includes(active.data.current.file.type) ||
      false) &&
    job.status === 1;

  const fileTypeLabel = getFileTypeLabel(item?._class);
  const hasFile = value && value !== nullFile;

  // Update value when item changes
  useEffect(() => {
    if (!item?._objectPath || !fileOptions) return;
    const dbFileId = item._value?.dbFileId?._value?.trim();
    if (!dbFileId) {
      setValue(nullFile);
      return;
    }
    const selectedFile = fileOptions.find(
      (file) => file.uuid.replace(/-/g, "") === dbFileId.replace(/-/g, "")
    );
    setValue(selectedFile || nullFile);
  }, [item, fileOptions]);

  // Reset expansion when children disappear
  useEffect(() => {
    if (!hasChildren) setIsManuallyExpanded(false);
  }, [hasChildren]);

  // Event handlers
  const handleFileSelect = useCallback(
    async (
      event: SyntheticEvent,
      selectedFile: CCP4i2File | null,
      reason: AutocompleteChangeReason
    ) => {
      const objectPath = item?._objectPath;
      if (!objectPath || !projects) return;

      const parameterArg =
        reason === "clear" || selectedFile === nullFile
          ? { value: null, object_path: objectPath }
          : fileItemToParameterArg(
              selectedFile!,
              objectPath,
              projectJobs || [],
              projects
            );

      setValue(selectedFile || nullFile);
      setInFlight(true);

      try {
        const result = await setParameter(parameterArg);
        if (result?.success) {
          const updatedItem = result.data?.updated_item ?? {
            ...result.data,
            _objectPath: objectPath,
          };
          onChange?.(updatedItem);
        }
      } catch (error) {
        console.error("Error setting parameter:", error);
        alert(`Error: ${error}`);
      } finally {
        setInFlight(false);
        await Promise.all([mutateContainer(), mutateContent(), mutateDigest()]);
      }
    },
    [
      item?._objectPath,
      projects,
      projectJobs,
      fileItemToParameterArg,
      setParameter,
      onChange,
      setInFlight,
      mutateContainer,
      mutateContent,
      mutateDigest,
    ]
  );

  const handleFileChange = useCallback(
    (event: ChangeEvent<HTMLInputElement>) => {
      setFiles?.(event.currentTarget.files);
    },
    [setFiles]
  );

  const handleMenuClick = useCallback(
    (event: React.MouseEvent<HTMLButtonElement>) => {
      event.stopPropagation();
      event.preventDefault();
      setFileMenuAnchorEl(event.currentTarget);
      setFile(value);
    },
    [setFileMenuAnchorEl, setFile, value]
  );

  const handleClear = useCallback(
    (event: React.MouseEvent) => {
      event.stopPropagation();
      handleFileSelect(event as unknown as SyntheticEvent, nullFile, "clear");
    },
    [handleFileSelect]
  );

  const handleToggle = useCallback(() => {
    if (!hasValidationError) setIsManuallyExpanded((prev) => !prev);
  }, [hasValidationError]);

  const getOptionLabel = useCallback(
    (option: CCP4i2File) => {
      const fileJob = projectJobs?.find((job) => job.id === option.job);
      return fileJob
        ? `${fileJob.number}: ${option.annotation}`
        : option.annotation;
    },
    [projectJobs]
  );

  // Loading and visibility checks
  if (!projectFiles || !projectJobs) return <LinearProgress />;
  if (!isVisible) return null;

  const canUpload = job.status === 1;
  const canFetch = qualifiers?.downloadModes?.length > 0 && job.status === 1;

  return (
    <Box
      sx={{
        mx: FIELD_SPACING.marginLeft,
        my: 0.5,
      }}
    >
      <Box
        ref={setNodeRef}
        sx={{
          border: 1,
          borderColor: hasError ? "error.main" : hasFile ? "primary.main" : "divider",
          borderLeftWidth: hasError ? 4 : hasFile ? 3 : 1,
          borderRadius: 1,
          bgcolor: isOver
            ? isValidDrop
              ? "success.light"
              : "error.light"
            : "background.paper",
          transition: "all 0.2s",
          "&:hover": {
            borderColor: hasError ? "error.main" : "primary.main",
          },
        }}
      >
        {/* Main row */}
        <Stack direction="row" alignItems="center" sx={{ p: 1 }}>
          {/* File type icon */}
          <Avatar
            src={`/qticons/${item?._class?.slice(1)}.png`}
            alt={item?._class || "File type"}
            sx={{
              width: 40,
              height: 40,
              mr: 1.5,
              flexShrink: 0,
              bgcolor: hasFile ? "primary.light" : "action.hover",
            }}
          />

          {/* File selector */}
          <Box sx={{ flex: 1, minWidth: 0 }}>
            <Autocomplete
              disabled={isDisabled}
              sx={{ ...sx }}
              size="small"
              value={value}
              onChange={handleFileSelect}
              options={[...fileOptions, nullFile]}
              getOptionLabel={getOptionLabel}
              getOptionKey={(option: CCP4i2File) => option.uuid}
              freeSolo={false}
              renderInput={(params) => (
                <TextField
                  {...params}
                  error={hasError}
                  slotProps={{
                    inputLabel: { shrink: true, disableAnimation: true },
                  }}
                  label={guiLabel}
                  size="small"
                  placeholder={`Select ${fileTypeLabel} file...`}
                />
              )}
              title={item?._objectPath || item?._className || "File selector"}
            />
          </Box>

          {/* Action buttons */}
          <Stack direction="row" spacing={0.5} sx={{ ml: 1, flexShrink: 0 }}>
            {canUpload && (
              <InputFileUpload
                sx={{ minWidth: "auto" }}
                disabled={isDisabled}
                accept={fileConfig.acceptedExtensions}
                handleFileChange={handleFileChange}
              />
            )}

            {canFetch && (
              <InputFileFetch
                sx={{ minWidth: "auto" }}
                disabled={isDisabled}
                modes={qualifiers.downloadModes}
                handleFileChange={handleFileChange}
                onChange={onChange}
                item={item}
              />
            )}

            {hasFile && (
              <>
                <Tooltip title="File options">
                  <IconButton
                    size="small"
                    onClick={handleMenuClick}
                    aria-label="Open file menu"
                  >
                    <MenuIcon fontSize="small" />
                  </IconButton>
                </Tooltip>
                {!isDisabled && (
                  <Tooltip title="Clear selection">
                    <IconButton
                      size="small"
                      onClick={handleClear}
                      aria-label="Clear file selection"
                      color="default"
                    >
                      <Clear fontSize="small" />
                    </IconButton>
                  </Tooltip>
                )}
              </>
            )}

            {hasChildren && (
              <Tooltip
                title={
                  hasValidationError
                    ? "Options expanded due to validation error"
                    : isExpanded
                      ? "Collapse options"
                      : "Expand options"
                }
              >
                <IconButton
                  onClick={handleToggle}
                  size="small"
                  disabled={hasValidationError}
                  sx={{
                    transition: "transform 0.2s ease-in-out",
                    transform: isExpanded ? "rotate(90deg)" : "rotate(0deg)",
                    opacity: hasValidationError ? 0.6 : 1,
                  }}
                  aria-label={isExpanded ? "Collapse options" : "Expand options"}
                >
                  <ChevronRightIcon fontSize="small" />
                </IconButton>
              </Tooltip>
            )}

            <ErrorTrigger item={item} job={job} />
          </Stack>
        </Stack>

        {/* Expandable children */}
        {hasChildren && (
          <ExpandableSection
            expanded={isExpanded}
            onToggle={(expanded) => setIsManuallyExpanded(expanded)}
            forceExpanded={hasValidationError || forceExpanded}
            hasError={hasValidationError}
            hideTitle
            forceExpandedTitle="Required Options (Error)"
          >
            <Box sx={{ px: 1, pb: 1 }}>{children}</Box>
          </ExpandableSection>
        )}
      </Box>
    </Box>
  );
};
