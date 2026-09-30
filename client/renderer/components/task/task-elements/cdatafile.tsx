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
  Divider,
  IconButton,
  LinearProgress,
  ListItemIcon,
  ListItemText,
  Menu,
  MenuItem,
  Stack,
  TextField,
  Tooltip,
} from "@mui/material";
import { alpha, SxProps, Theme } from "@mui/material/styles";
import {
  AccountTree as AccountTreeIcon,
  MoreVert as MoreVertIcon,
  ChevronRight as ChevronRightIcon,
  ContentCopy,
  ContentPaste,
  DeleteOutline,
  Download,
  HelpOutline,
  Preview,
  SaveAlt,
} from "@mui/icons-material";
import { useDndContext, useDroppable } from "@dnd-kit/core";

import { doDownload, useApi } from "../../../api";
import { ACTIVE_JOB_STATUSES, useJob, useProject, useProjectFiles } from "../../../utils";
import { CCP4i2TaskElementProps } from "./task-element";
import { File as CCP4i2File, nullFile, Project } from "../../../types/models";
import { FileMenuExtraItem, useFileMenu } from "../../../providers/file-context-menu";
import { ErrorTrigger } from "./error-info";
import { InputFileFetch } from "./input-file-fetch";
import { InputFileUpload } from "./input-file-upload";
import { FIELD_SPACING } from "./field-sizes";
import { ExpandableSection } from "./expandable-section";
import { BrowseProjectFilesDialog } from "./browse-project-files-dialog";
import { useContainerField } from "./hooks/useContainerField";

/**
 * Content flag conversion capability mapping.
 * Defines which content types can be converted to which target types.
 * Mirrors CAN_CONVERT_TO from server/ccp4i2/core/CCP4XtalData.py
 *
 * Content flags: IPAIR=1, FPAIR=2, IMEAN=3, FMEAN=4
 */
const CAN_CONVERT_TO: Record<number, number[]> = {
  1: [1, 2, 3, 4], // IPAIR can convert to IPAIR, FPAIR, IMEAN, FMEAN
  2: [2, 4],       // FPAIR can convert to FPAIR, FMEAN
  3: [3, 4],       // IMEAN can convert to IMEAN, FMEAN
  4: [4],          // FMEAN can convert to FMEAN only
};

/**
 * Check if a file's content type can be converted to any of the required types.
 * @param fileContent The file's content flag (1-4)
 * @param requiredFlags Array of required content flags
 * @returns true if the file can provide data in any required format
 */
const canConvertToRequired = (
  fileContent: number | null | undefined,
  requiredFlags: number[]
): boolean => {
  if (fileContent == null) return false;
  const convertibleTo = CAN_CONVERT_TO[fileContent];
  if (!convertibleTo) return false;
  return requiredFlags.some((required) => convertibleTo.includes(required));
};

/**
 * Shared look for the boxed group of file-source buttons (file system /
 * project database / internet): one outline round the lot, hairline dividers
 * between the buttons, square corners on the buttons themselves.
 */
const FILE_SOURCE_GROUP_SX: SxProps<Theme> = {
  display: "flex",
  alignItems: "stretch",
  // Match the outline weight of an outlined TextField/Button (23% ink). The
  // theme's `divider` (12%) is too faint to read as a group at this size.
  border: "1px solid",
  borderColor: (theme) => alpha(theme.palette.text.primary, 0.23),
  borderRadius: 1,
  overflow: "hidden",
  "& .MuiIconButton-root": { borderRadius: 0 },
  // Hairlines between the buttons, same ink as the surrounding outline
  "& .MuiDivider-root": { borderColor: "inherit" },
};

/** Extra menu items injected by subtype widgets (e.g. CPdbDataFile "Select atoms") */
export interface IconMenuItem {
  label: string;
  icon?: ReactNode;
  onClick: () => void;
  disabled?: boolean;
  /** Render a divider before this item */
  divider?: boolean;
  checkable?: boolean;
  checked?: boolean;
}

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
  /** Additional icon-menu items from subtype widgets (inserted before separator) */
  iconMenuItems?: IconMenuItem[];
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
  iconMenuItems,
}) => {
  const api = useApi();
  const { fileItemToParameterArg, mutateContainer, mutateFileDigest, mutateFileContent } =
    useJob(job.id);

  const {
    item,
    isDisabled,
    isVisible,
    validationColor,
    commit,
  } = useContainerField<any>({
    job,
    itemName,
    visibility,
    disabled: disabledProp,
    onChange,
  });
  const { setFileMenuAnchorEl, setFile, setExtraMenuItems } = useFileMenu();

  // Data and state
  // Poll for files only while the job is active, so task widgets see newly
  // created output files without polling forever once the job is done.
  const { files: projectFiles } = useProjectFiles(
    job.project,
    ACTIVE_JOB_STATUSES.includes(job.status)
  );
  const { jobs: projectJobs } = useProject(job.project);
  const { data: projects } = api.get<Project[]>("projects");
  // Invalidate, never subscribe: this element renders once per file in a
  // task, and subscribing here fetched a digest and the file's content for
  // every one of them on mount. The helpers revalidate only the mounted
  // subscribers (the elements that show a digest or content).
  const mutateDigest = useCallback(
    () => mutateFileDigest(item?._objectPath ?? ""),
    [mutateFileDigest, item?._objectPath]
  );
  const mutateContent = mutateFileContent;
  const [value, setValue] = useState<CCP4i2File>(nullFile);
  const [isManuallyExpanded, setIsManuallyExpanded] = useState(false);
  const [browseDialogOpen, setBrowseDialogOpen] = useState(false);
  // iconMenuAnchorEl removed — avatar now opens the shared FileMenu

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

  // Compare job numbers hierarchically (e.g., "1.2.3" < "1.2.4" < "2")
  // Returns positive if a should come after b (higher job numbers first)
  const compareJobNumbers = useCallback(
    (jobIdA: number, jobIdB: number): number => {
      const jobA = projectJobs?.find((j) => j.id === jobIdA);
      const jobB = projectJobs?.find((j) => j.id === jobIdB);
      if (!jobA || !jobB) return 0;

      const aParts = jobA.number.split(".").map(Number);
      const bParts = jobB.number.split(".").map(Number);

      for (let i = 0; i < Math.max(aParts.length, bParts.length); i++) {
        const aVal = aParts[i] || 0;
        const bVal = bParts[i] || 0;
        if (aVal !== bVal) return bVal - aVal; // Descending order (higher first)
      }
      return 0;
    },
    [projectJobs]
  );

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
        // Check if file's content type can be CONVERTED to any required type
        // (e.g., IPAIR file can satisfy FMEAN requirement via conversion)
        // null, undefined, or empty array means no filtering
        const hasValidContentFlag =
          !fileConfig.requiredContentFlag ||
          fileConfig.requiredContentFlag.length === 0 ||
          canConvertToRequired(file.content, fileConfig.requiredContentFlag);
        return isValidType && isNotParentJob && hasValidContentFlag;
      })
      .sort((a, b) => compareJobNumbers(a.job, b.job));
  }, [projectFiles, projectJobs, fileConfig.allowedTypes, fileConfig.requiredContentFlag, compareJobNumbers]);

  const borderColor = validationColor;
  const hasError = borderColor === "error.light";
  const computedValidationError = hasError;
  const hasValidationError = overrideValidationError ?? computedValidationError;
  const hasChildren = React.Children.toArray(children).length > 0;
  const isExpanded = hasValidationError || isManuallyExpanded || forceExpanded;

  const guiLabel =
    qualifiers?.guiLabel ?? item?._objectPath?.split(".").at(-1) ?? "";

  // Drag and drop — drop target (only for pending job inputs)
  const { isOver, setNodeRef } = useDroppable({
    id: `job_${job.uuid}_${itemName}`,
    data: { job, item },
  });

  const fileIsSet = value && value !== nullFile;

  // Native HTML5 drag for the file icon — enables OS drag-out and within-window drops
  const handleFileDragStart = useCallback(
    (e: React.DragEvent) => {
      if (!fileIsSet || !value) return;

      const ref = {
        ccp4i2_file: true,
        uuid: value.uuid,
        id: value.id,
        name: value.name,
        type: value.type,
        sub_type: value.sub_type,
        content: value.content,
        annotation: value.annotation,
        job: value.job,
        job_param_name: value.job_param_name,
      };
      e.dataTransfer.setData("application/ccp4i2-file", JSON.stringify(ref));
      e.dataTransfer.effectAllowed = "copy";

      // Write to clipboard for cross-window paste
      navigator.clipboard.writeText(JSON.stringify(ref)).catch(() => {});

    },
    [fileIsSet, value, projectJobs]
  );
  const { active } = useDndContext();
  const isValidDrop =
    active?.data?.current?.file &&
    (fileConfig.allowedTypes?.includes(active.data.current.file.type) ||
      false) &&
    job.status === 1;

  const fileTypeLabel = getFileTypeLabel(item?._class);
  const hasFile = value && value !== nullFile;

  // Get the current dbFileId from the item
  const dbFileId = item?._value?.dbFileId?._value?.trim() || null;

  // Check if the selected file exists in fileOptions
  const selectedFileInOptions = useMemo(() => {
    if (!dbFileId || !fileOptions) return null;
    return fileOptions.find(
      (file) => file.uuid.replace(/-/g, "") === dbFileId.replace(/-/g, "")
    );
  }, [dbFileId, fileOptions]);

  // A file from a subjob or another project is not in fileOptions, and the
  // Autocomplete needs an option object to render its value. Everything that
  // needs is already in the parameter -- dbFileId, baseName, annotation --
  // and for a file outside this project getOptionLabel falls back to
  // `annotation || name` anyway, because its job is not among projectJobs.
  //
  // So the display costs nothing. This matters at scale: a pandda_campaign
  // job carries 146 such references, and fetching each one took 173 ms, so
  // the interface spent 25 seconds of round trips filling in fields that the
  // parameters already described.
  const optionFromParameter = useMemo<CCP4i2File | null>(() => {
    if (!dbFileId || selectedFileInOptions) return null;
    const contents = item?._value ?? {};
    return {
      uuid: dbFileId,
      name: contents?.baseName?._value ?? "",
      annotation: contents?.annotation?._value ?? "",
    } as CCP4i2File;
  }, [dbFileId, selectedFileInOptions, item]);

  // The full record is only needed by things the user does to a file -- the
  // drag payload below wants id, type, sub_type, content and job, which the
  // parameter does not carry. dragstart cannot await, so this is armed when
  // the pointer reaches the file, which always precedes a drag.
  const [needsFileRecord, setNeedsFileRecord] = useState(false);
  const { data: fetchedFile } = api.get<CCP4i2File>(
    dbFileId && !selectedFileInOptions && needsFileRecord
      ? `files_by_uuid/${dbFileId}/`
      : null
  );

  // Combine fileOptions with the selected file when it is not one of them:
  // the record once fetched, otherwise the one described by the parameter.
  const displayOptions = useMemo(() => {
    const selected = fetchedFile ?? optionFromParameter;
    if (!selected || selectedFileInOptions) return fileOptions;
    return [selected, ...fileOptions];
  }, [fileOptions, fetchedFile, optionFromParameter, selectedFileInOptions]);

  // Update value when item changes
  useEffect(() => {
    if (!item?._objectPath) return;
    if (!dbFileId) {
      setValue(nullFile);
      return;
    }
    // fileOptions first, then the fetched record, then what the parameter
    // itself says the file is.
    const selectedFile =
      selectedFileInOptions ||
      (fetchedFile &&
        fetchedFile.uuid?.replace(/-/g, "") === dbFileId.replace(/-/g, "")
        ? fetchedFile
        : optionFromParameter);
    setValue(selectedFile || nullFile);
  }, [item, dbFileId, selectedFileInOptions, fetchedFile, optionFromParameter]);

  // Reset expansion when children disappear
  useEffect(() => {
    if (!hasChildren) setIsManuallyExpanded(false);
  }, [hasChildren]);

  // Event handlers
  const handleFileSelect = useCallback(
    async (
      _event: SyntheticEvent,
      selectedFile: CCP4i2File | null,
      reason: AutocompleteChangeReason
    ) => {
      if (!item?._objectPath || !projects) return;

      const writeValue =
        reason === "clear" || selectedFile === nullFile
          ? null
          : fileItemToParameterArg(
              selectedFile!,
              item._objectPath,
              projectJobs || [],
              projects
            ).value;

      const previous = value;
      setValue(selectedFile || nullFile);

      const result = await commit(writeValue);
      if (result && !result.success) {
        setValue(previous);
      }
      await Promise.all([mutateContainer(), mutateContent(), mutateDigest()]);
    },
    [
      item?._objectPath,
      projects,
      projectJobs,
      fileItemToParameterArg,
      commit,
      value,
      mutateContainer,
      mutateContent,
      mutateDigest,
    ]
  );

  const handleFileChange = useCallback(
    (event: ChangeEvent<HTMLInputElement>) => {
      const input = event.currentTarget;
      const picked = input.files;
      // Copy the selection into an independent FileList before touching the
      // input, so resetting its value below can't empty what we hand on.
      let files: FileList | null = picked;
      if (picked && picked.length && typeof DataTransfer !== "undefined") {
        const dt = new DataTransfer();
        for (let i = 0; i < picked.length; i++) dt.items.add(picked[i]);
        files = dt.files;
      }
      setFiles?.(files);
      // Clear the input's value so the SAME file can be chosen again after a
      // Clear: an <input type="file"> fires no change event when its value is
      // unchanged, which otherwise leaves a cleared field impossible to
      // repopulate by re-picking the same path.
      input.value = "";
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

  const handleBrowseFileSelect = useCallback(
    (file: CCP4i2File) => {
      // Use the same path as the autocomplete selection
      handleFileSelect(
        {} as React.SyntheticEvent,
        file,
        "selectOption"
      );
    },
    [handleFileSelect]
  );

  const handleToggle = useCallback(() => {
    if (!hasValidationError) setIsManuallyExpanded((prev) => !prev);
  }, [hasValidationError]);

  // Open the shared file context menu with cdatafile-specific extra items
  const handleIconContextMenu = useCallback(
    (event: React.MouseEvent<HTMLElement>) => {
      event.preventDefault();
      event.stopPropagation();

      // Build extra items specific to this cdatafile widget
      const extras: FileMenuExtraItem[] = [];

      if (!isDisabled) {
        extras.push({
          key: "clear",
          label: "Clear",
          icon: <DeleteOutline fontSize="small" />,
          onClick: () => handleFileSelect({} as React.SyntheticEvent, nullFile, "clear"),
        });
      }

      // Subtype-specific items (e.g. "Select atoms")
      if (iconMenuItems) {
        for (const extra of iconMenuItems) {
          if (extra.divider) {
            extras.push({ key: `div-${extra.label}`, label: "", onClick: () => {}, divider: true });
          }
          extras.push({
            key: extra.label,
            label: extra.label,
            icon: extra.icon,
            onClick: extra.onClick,
            disabled: extra.disabled,
          });
        }
      }

      // Paste
      if (!isDisabled) {
        extras.push({
          key: "paste",
          label: "Paste",
          icon: <ContentPaste fontSize="small" />,
          onClick: async () => {
            try {
              const text = await navigator.clipboard.readText();
              const parsed = JSON.parse(text);
              if (parsed?.ccp4i2_file && parsed.uuid && item?._objectPath && projectJobs && projects) {
                const fileRef = { ...parsed, exports: [], fileimport: -1, file_uses: [] } as CCP4i2File;
                const arg = fileItemToParameterArg(fileRef, item._objectPath, projectJobs, projects);
                await commit(arg.value);
                await mutateContainer();
              }
            } catch { /* ignore invalid clipboard */ }
          },
        });
      }

      // Help
      if (qualifiers?.helpFile) {
        extras.push({
          key: "help",
          label: "Help",
          icon: <HelpOutline fontSize="small" />,
          onClick: () => window.open(`https://www.ccp4.ac.uk/html/${qualifiers.helpFile}.html`, "_blank", "noopener,noreferrer"),
        });
      }

      setExtraMenuItems(extras);
      setFile(hasFile ? value : null);
      setFileMenuAnchorEl(event.currentTarget);
    },
    [isDisabled, hasFile, value, iconMenuItems, qualifiers, item, projectJobs, projects,
     handleFileSelect, fileItemToParameterArg, commit, mutateContainer,
     setExtraMenuItems, setFile, setFileMenuAnchorEl]
  );

  // Native HTML5 drag and drop for filesystem files
  const [nativeDragOver, setNativeDragOver] = useState(false);

  const handleNativeDragOver = useCallback(
    (e: React.DragEvent) => {
      // Handle filesystem files or ccp4i2 file references, but not internal @dnd-kit drags
      const hasCcp4i2File = e.dataTransfer.types.includes("application/ccp4i2-file");
      const hasFiles = e.dataTransfer.types.includes("Files");
      if ((hasCcp4i2File || hasFiles) && !active && !isDisabled) {
        e.preventDefault();
        e.stopPropagation();
        setNativeDragOver(true);
      }
    },
    [active, isDisabled]
  );

  const handleNativeDragLeave = useCallback((e: React.DragEvent) => {
    e.preventDefault();
    setNativeDragOver(false);
  }, []);

  const handleNativeDrop = useCallback(
    async (e: React.DragEvent) => {
      e.preventDefault();
      e.stopPropagation();
      setNativeDragOver(false);

      // Handle ccp4i2 file reference drop (from another file element)
      const ccp4i2Data = e.dataTransfer.getData("application/ccp4i2-file");
      if (ccp4i2Data && !isDisabled && item?._objectPath && projectJobs && projects) {
        try {
          const parsed = JSON.parse(ccp4i2Data);
          if (parsed?.ccp4i2_file && parsed.uuid) {
            const fileRef = {
              ...parsed,
              exports: [],
              fileimport: -1,
              file_uses: [],
            } as CCP4i2File;
            const arg = fileItemToParameterArg(fileRef, item._objectPath, projectJobs, projects);
            await commit(arg.value);
            await mutateContainer();
            return;
          }
        } catch {
          // Invalid JSON — fall through to filesystem drop
        }
      }

      // Handle filesystem file drop
      if (e.dataTransfer.files?.length > 0 && setFiles && !isDisabled) {
        setFiles(e.dataTransfer.files);
      }
    },
    [setFiles, isDisabled, item, projectJobs, projects, fileItemToParameterArg, commit, mutateContainer]
  );

  // A file registered with no annotation (dimple's final.pdb, for one) still
  // has a name; labelling it by annotation alone rendered it as nothing, which
  // looked like an empty picker.
  const getOptionLabel = useCallback(
    (option: CCP4i2File) => {
      const fileJob = projectJobs?.find((job) => job.id === option.job);
      return fileJob
        ? `${fileJob.number}: ${option.annotation || option.name}`
        : option.annotation || option.name;
    },
    [projectJobs]
  );

  // Loading and visibility checks
  if (!projectFiles || !projectJobs) return <LinearProgress />;
  if (!isVisible) return null;

  const canUpload = job.status === 1;
  const canFetch = qualifiers?.downloadModes?.length > 0 && job.status === 1;

  // The three logically equivalent ways of choosing a file. Fetch-from-internet
  // is only offered for file types that declare downloadModes, so the group can
  // hold two or three buttons.
  const fileSourceButtons: { key: string; node: ReactNode }[] = [];
  if (canUpload) {
    fileSourceButtons.push({
      key: "filesystem",
      node: (
        <InputFileUpload
          disabled={isDisabled}
          accept={fileConfig.acceptedExtensions}
          handleFileChange={handleFileChange}
        />
      ),
    });
    fileSourceButtons.push({
      key: "database",
      node: (
        <Tooltip title="Browse the project hierarchy">
          <span>
            <IconButton
              size="small"
              onClick={() => setBrowseDialogOpen(true)}
              disabled={isDisabled}
              aria-label="Browse project files"
            >
              <AccountTreeIcon fontSize="small" />
            </IconButton>
          </span>
        </Tooltip>
      ),
    });
  }
  if (canFetch) {
    fileSourceButtons.push({
      key: "internet",
      node: (
        <InputFileFetch
          disabled={isDisabled}
          modes={qualifiers.downloadModes}
          onChange={onChange}
          item={item}
        />
      ),
    });
  }

  return (
    <Box
      ref={setNodeRef}
      onDragOver={handleNativeDragOver}
      onDragLeave={handleNativeDragLeave}
      onDrop={handleNativeDrop}
      sx={{
        mx: FIELD_SPACING.marginLeft,
        my: 0,
        bgcolor: nativeDragOver
          ? "success.light"
          : isOver
            ? isValidDrop
              ? "success.light"
              : "error.light"
            : "transparent",
        borderRadius: 1,
        transition: "background-color 0.2s",
      }}
    >
      {/* Main row */}
      <Stack direction="row" alignItems="center">
        {/* File type icon — native HTML5 draggable + right-click context menu */}
        <Avatar
          draggable={!!hasFile}
          // The drag payload needs the file's full record, and dragstart
          // cannot wait for a fetch. Arming it when the pointer arrives (or
          // the icon takes focus) means the one file being touched is
          // resolved, rather than every file on the page at mount.
          onPointerEnter={() => setNeedsFileRecord(true)}
          onFocus={() => setNeedsFileRecord(true)}
          onDragStart={handleFileDragStart}
          src={`/svgicons/${item?._class?.slice(1)}.svg`}
          alt={item?._class || "File type"}
          onContextMenu={handleIconContextMenu}
          sx={{
            width: 32,
            height: 32,
            mr: 1,
            flexShrink: 0,
            bgcolor: hasFile ? "primary.light" : "action.hover",
            cursor: hasFile ? "grab" : "context-menu",
            transition: "box-shadow 0.2s ease",
            "&:hover": hasFile
              ? {
                  boxShadow: (theme) =>
                    `0 0 0 3px ${alpha(theme.palette.primary.main, 0.5)}`,
                }
              : {},
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
            options={[...displayOptions, nullFile]}
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

        {/* Action buttons: the three ways of choosing a file (file system,
            project database, internet) are boxed together as one group; the
            per-file utility icons that follow are unboxed. */}
        <Stack
          direction="row"
          spacing={0.5}
          alignItems="center"
          sx={{ ml: 1, flexShrink: 0 }}
        >
          {fileSourceButtons.length > 0 && (
            <Box sx={FILE_SOURCE_GROUP_SX}>
              {fileSourceButtons.map(({ key, node }, index) => (
                <React.Fragment key={key}>
                  {index > 0 && <Divider orientation="vertical" flexItem />}
                  {node}
                </React.Fragment>
              ))}
            </Box>
          )}

          {hasFile && (
            <Tooltip title="File options">
              <IconButton
                size="small"
                onClick={handleMenuClick}
                aria-label="Open file menu"
              >
                <MoreVertIcon fontSize="small" />
              </IconButton>
            </Tooltip>
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
          sx={{
            ml: 5,
            borderTop: "none",
            borderLeft: "2px solid",
            borderLeftColor: hasValidationError ? "error.main" : "divider",
            borderRadius: 0,
            borderBottomLeftRadius: 0,
            borderBottomRightRadius: 0,
            px: 1.5,
          }}
        >
          {children}
        </ExpandableSection>
      )}

      {/* Browse files from other projects dialog */}
      <BrowseProjectFilesDialog
        open={browseDialogOpen}
        onClose={() => setBrowseDialogOpen(false)}
        onFileSelect={handleBrowseFileSelect}
        allowedTypes={fileConfig.allowedTypes}
        requiredContentFlag={fileConfig.requiredContentFlag}
        fileTypeLabel={fileTypeLabel}
      />

      {/* File context menu is now the shared FileMenu from file-context-menu.tsx */}
      {/* Extra items (Clear, Copy, Paste, Help, subtype items) are injected via setExtraMenuItems */}
    </Box>
  );
};
