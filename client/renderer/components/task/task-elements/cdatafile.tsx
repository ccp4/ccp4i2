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
  Button,
  Collapse,
  IconButton,
  LinearProgress,
  Stack,
  TextField,
  Typography,
} from "@mui/material";
import {
  Menu as MenuIcon,
  ExpandMore as ExpandMoreIcon,
  ChevronRight as ChevronRightIcon,
} from "@mui/icons-material";
import { useDndContext, useDroppable } from "@dnd-kit/core";

import { useApi } from "../../../api";
import { useJob, useProject } from "../../../utils";
import { CCP4i2TaskElementProps } from "./task-element";
import { File as CCP4i2File, nullFile, Project } from "../../../types/models";
import { useTaskInterface } from "../../../providers/task-provider";
import { useFileMenu } from "../../../providers/file-context-menu";
import { ErrorTrigger } from "./error-info";
import { InputFileFetch } from "./input-file-fetch";
import { InputFileUpload } from "./input-file-upload";
import { FIELD_SPACING } from "./field-sizes";

const BORDER_RADIUS_STYLES = {
  none: { borderRadius: 0 },
  left: { borderTopLeftRadius: "0.5rem", borderBottomLeftRadius: "0.5rem" },
  right: { borderTopRightRadius: "0.5rem", borderBottomRightRadius: "0.5rem" },
  full: { borderRadius: "0.5rem" },
} as const;

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
  const { files: projectFiles, jobs: projectJobs } = useProject(job.project);
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
    return { allowedTypes, acceptedExtensions };
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
        return isValidType && isNotParentJob;
      })
      .sort((a, b) => b.job - a.job);
  }, [projectFiles, projectJobs, fileConfig.allowedTypes]);

  const borderColor = getValidationColor(item);
  const computedValidationError = borderColor === "error.light";
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

  const backgroundColor = isOver
    ? isValidDrop
      ? "success.light"
      : "error.light"
    : "background.paper";

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
        // Call onChange on success - pass updated_item if available (from upload_file_param),
        // otherwise pass the result.data (from set_parameter) with the objectPath added
        if (result?.success) {
          const updatedItem = result.data?.updated_item ?? { ...result.data, _objectPath: objectPath };
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

  const showMenuButton = value && value !== nullFile;
  const canUpload = job.status === 1;
  const canFetch = qualifiers?.downloadModes?.length > 0 && job.status === 1;

  return (
    <Stack
      sx={{
        border: "3px solid",
        borderColor,
        backgroundColor,
        borderRadius: "0.5rem",
        mx: FIELD_SPACING.marginLeft,
        my: 0.5,
      }}
      direction="column"
    >
      <Stack ref={setNodeRef} direction="row" alignItems="center">
        <Avatar
          src={`/api/proxy/djangostatic/qticons/${item?._class?.slice(1)}.png`}
          alt={item?._class || "File type"}
        />

        <Autocomplete
          disabled={isDisabled}
          sx={{ m: 1, flex: 1, minWidth: 0, ...sx }}
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
              error={borderColor === "error.light"}
              slotProps={{
                inputLabel: { shrink: true, disableAnimation: true },
              }}
              label={guiLabel}
              size="small"
            />
          )}
          title={item?._objectPath || item?._className || "File selector"}
        />

        <Stack direction="row">
          {canUpload && (
            <InputFileUpload
              sx={{
                my: 2,
                ...BORDER_RADIUS_STYLES.none,
                "&:first-of-type": BORDER_RADIUS_STYLES.left,
                "&:last-of-type": BORDER_RADIUS_STYLES.right,
              }}
              disabled={isDisabled}
              accept={fileConfig.acceptedExtensions}
              handleFileChange={handleFileChange}
            />
          )}

          {canFetch && (
            <InputFileFetch
              sx={{
                my: 2,
                ...BORDER_RADIUS_STYLES.none,
                "&:first-of-type": BORDER_RADIUS_STYLES.left,
                "&:last-of-type": BORDER_RADIUS_STYLES.right,
              }}
              disabled={isDisabled}
              modes={qualifiers.downloadModes}
              handleFileChange={handleFileChange}
              onChange={onChange}
              item={item}
            />
          )}

          {showMenuButton && (
            <Button
              disabled={false}
              role="button"
              variant="outlined"
              tabIndex={-1}
              size="small"
              startIcon={<MenuIcon fontSize="small" />}
              sx={{
                my: 2,
                ...BORDER_RADIUS_STYLES.none,
                "&:first-of-type": BORDER_RADIUS_STYLES.left,
                "&:last-of-type": BORDER_RADIUS_STYLES.right,
              }}
              onClick={handleMenuClick}
              aria-label="Open file menu"
            />
          )}
        </Stack>

        {hasChildren && (
          <IconButton
            onClick={handleToggle}
            size="small"
            disabled={hasValidationError}
            sx={{
              ml: 1,
              transition: "transform 0.2s ease-in-out",
              transform: isExpanded ? "rotate(90deg)" : "rotate(0deg)",
              opacity: hasValidationError ? 0.6 : 1,
            }}
            aria-label={
              hasValidationError
                ? "Options expanded due to validation error"
                : isExpanded
                  ? "Collapse options"
                  : "Expand options"
            }
          >
            {isExpanded ? (
              <ExpandMoreIcon fontSize="small" />
            ) : (
              <ChevronRightIcon fontSize="small" />
            )}
          </IconButton>
        )}

        <ErrorTrigger item={item} job={job} />
      </Stack>

      {hasChildren && (
        <Collapse in={isExpanded} timeout={200}>
          <Stack
            sx={{
              px: 2,
              pb: 1,
              pt: 0,
              backgroundColor: hasValidationError
                ? "error.lighter"
                : "background.paper",
              borderTop: "1px solid",
              borderTopColor: hasValidationError ? "error.light" : "divider",
              borderBottomLeftRadius: "0.4rem",
              borderBottomRightRadius: "0.4rem",
            }}
            spacing={0.5}
          >
            <Typography
              variant="caption"
              color={hasValidationError ? "error.main" : "text.secondary"}
              sx={{
                fontWeight: 500,
                textTransform: "uppercase",
                letterSpacing: 0.5,
                mb: 0.5,
              }}
            >
              {hasValidationError
                ? "Required Options (Error)"
                : "Additional Options"}
            </Typography>
            {children}
          </Stack>
        </Collapse>
      )}
    </Stack>
  );
};
