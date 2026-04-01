/*
 * Copyright (C) 2025-2026 Newcastle University
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
import React, {
  forwardRef,
  useCallback,
  useMemo,
  useState,
  createContext,
  useContext,
} from "react";
import {
  Box,
  Button,
  Checkbox,
  Chip,
  IconButton,
  List,
  ListItem,
  Paper,
  Skeleton,
  Stack,
  Toolbar,
  Tooltip,
  Typography,
} from "@mui/material";
import SearchField from "./search-field";
import {
  RichTreeView,
  TreeItem2Content,
  TreeItem2GroupTransition,
  TreeItem2Icon,
  TreeItem2IconContainer,
  TreeItem2Label,
  TreeItem2LabelInput,
  TreeItem2Props,
  TreeItem2Root,
  useTreeItem2,
} from "@mui/x-tree-view";
import { CheckBoxOutlined, Clear, Delete, Menu as MenuIcon } from "@mui/icons-material";
import { useDraggable } from "@dnd-kit/core";
import { useRouter } from "next/navigation";

import { EndpointFetch, useApi } from "../api";
import {
  Job,
  JobTreeNode,
  JobTreeResponse,
  File as DjangoFile,
} from "../types/models";
import { CCP4i2JobAvatar } from "./job-avatar";
import { FileAvatar } from "./file-avatar";
import { useJobMenu } from "../providers/job-context-menu";
import { useFileMenu } from "../providers/file-context-menu";
import { useRecentlyStartedJobs } from "../providers/recently-started-jobs-context";
import { useDeleteDialog } from "../providers/delete-dialog";
import { useSet } from "../hooks";

// =============================================================================
// Types
// =============================================================================

interface ClassicJobListProps {
  projectId: number;
  parent?: number;
  withSubtitles?: boolean;
}

/** Tree item type for RichTreeView - combines JobTreeNode with children */
interface TreeViewItem extends Omit<JobTreeNode, 'children'> {
  children: (TreeViewItem | DjangoFile)[];
}

interface TreeItemData {
  job?: JobTreeNode;
  file?: DjangoFile;
  isJob: boolean;
  displayLabel: string;
  timestamp?: string;
}

// =============================================================================
// Constants
// =============================================================================

/** Polling interval when jobs are actively running (2, 3, or 7) */
const ACTIVE_POLL_INTERVAL = 3000;
/** Polling interval when all jobs are idle */
const IDLE_POLL_INTERVAL = 30000;
/** Job statuses that indicate active work: Queued (2), Running (3), Running remotely (7) */
const ACTIVE_JOB_STATUSES = [2, 3, 7];

const TIME_DISPLAY_STYLE = { fontSize: "75%" };

// =============================================================================
// Context for sharing job tree data with tree items
// =============================================================================

interface JobTreeContextValue {
  jobsByUuid: Map<string, JobTreeNode>;
  jobsById: Map<number, JobTreeNode>;
  filesByUuid: Map<string, DjangoFile>;
  taskShortTitles: Map<string, string>;
  selectMode: boolean;
  selectedJobIds: Set<number>;
  toggleJobSelection: (jobId: number) => void;
}

const JobTreeContext = createContext<JobTreeContextValue>({
  jobsByUuid: new Map(),
  jobsById: new Map(),
  filesByUuid: new Map(),
  taskShortTitles: new Map(),
  selectMode: false,
  selectedJobIds: new Set(),
  toggleJobSelection: () => {},
});

// =============================================================================
// Custom hooks
// =============================================================================

/**
 * Check if any job in the tree has an active status (queued, running, running remotely)
 */
const hasActiveJobs = (jobTree: JobTreeNode[] | undefined): boolean => {
  if (!jobTree) return false;

  const checkNode = (node: JobTreeNode): boolean => {
    if (ACTIVE_JOB_STATUSES.includes(node.status)) return true;
    return node.children.some(checkNode);
  };

  return jobTree.some(checkNode);
};

/**
 * Fetch job tree data using consolidated endpoint.
 * Replaces 4 separate calls (jobs, files)
 * with a single optimized request.
 *
 * Uses adaptive polling: faster when jobs are running OR recently started,
 * slower when all jobs are idle.
 */
const useJobTree = (projectId: number) => {
  const api = useApi();
  // Track whether jobs are active to determine polling interval
  const [jobsActive, setJobsActive] = useState(false);

  // Check if any jobs were recently started (grace period for DB latency)
  const { hasRecentlyStartedJobsInProject } = useRecentlyStartedJobs();
  const hasRecentlyStarted = hasRecentlyStartedJobsInProject(projectId);

  const endpointFetch: EndpointFetch = {
    type: "projects",
    id: projectId,
    endpoint: "job_tree",
  };

  // Use adaptive polling interval based on job activity OR recent starts
  // This handles the race condition where DB status hasn't updated yet
  const pollInterval = (jobsActive || hasRecentlyStarted) ? ACTIVE_POLL_INTERVAL : IDLE_POLL_INTERVAL;

  const { data, isLoading, error, mutate } = api.get_endpoint<JobTreeResponse>(
    endpointFetch,
    pollInterval
  );

  // Update active state when data changes
  // Using useEffect for side effect (updating state based on data)
  React.useEffect(() => {
    if (data?.job_tree) {
      const active = hasActiveJobs(data.job_tree);
      setJobsActive(active);
    }
  }, [data?.job_tree]);

  return {
    jobTree: data?.job_tree,
    totalJobs: data?.total_jobs ?? 0,
    totalFiles: data?.total_files ?? 0,
    isLoading,
    error,
    mutate,
    jobsActive: jobsActive || hasRecentlyStarted, // Include recently started in active state
  };
};

/**
 * Build lookup maps for O(1) access to jobs and files by UUID.
 * This allows tree items to efficiently find their data.
 */
const useJobTreeLookups = (jobTree: JobTreeNode[] | undefined) => {
  return useMemo(() => {
    const jobsByUuid = new Map<string, JobTreeNode>();
    const jobsById = new Map<number, JobTreeNode>();
    const filesByUuid = new Map<string, DjangoFile>();

    if (!jobTree) {
      return { jobsByUuid, jobsById, filesByUuid };
    }

    // Recursively index all jobs and files
    const indexNode = (node: JobTreeNode) => {
      jobsByUuid.set(node.uuid, node);
      jobsById.set(node.id, node);
      node.files.forEach((file) => filesByUuid.set(file.uuid, file));
      node.children.forEach(indexNode);
    };

    jobTree.forEach(indexNode);

    return { jobsByUuid, jobsById, filesByUuid };
  }, [jobTree]);
};

/**
 * Check if a job node matches the search text.
 * Matches against job title, task_name, job number, and finish time.
 */
const jobMatchesFilter = (node: JobTreeNode, searchUpper: string): boolean => {
  if (node.title?.toUpperCase().includes(searchUpper)) return true;
  if (node.task_name?.toUpperCase().includes(searchUpper)) return true;
  if (node.number?.toUpperCase().includes(searchUpper)) return true;
  if (
    node.finish_time &&
    new Date(node.finish_time).getFullYear() > 1970 &&
    new Date(node.finish_time).toLocaleString().toUpperCase().includes(searchUpper)
  )
    return true;
  return false;
};

/**
 * Recursively filter a job tree, keeping nodes that match or have matching descendants.
 */
const filterJobTree = (
  nodes: JobTreeNode[],
  searchUpper: string
): JobTreeNode[] => {
  return nodes.reduce<JobTreeNode[]>((acc, node) => {
    const filteredChildren = filterJobTree(node.children, searchUpper);
    if (jobMatchesFilter(node, searchUpper) || filteredChildren.length > 0) {
      acc.push({ ...node, children: filteredChildren });
    }
    return acc;
  }, []);
};

/**
 * Transform job tree to format expected by RichTreeView.
 * Merges files into children array for each job.
 * Optionally filters by search text.
 */
const useTreeViewItems = (
  jobTree: JobTreeNode[] | undefined,
  filterText: string | null
): TreeViewItem[] => {
  return useMemo(() => {
    if (!jobTree) return [];

    // Apply filter first if search text is provided
    const filtered =
      filterText && filterText.trim().length > 0
        ? filterJobTree(jobTree, filterText.trim().toUpperCase())
        : jobTree;

    const transformNode = (node: JobTreeNode): TreeViewItem => {
      // Combine child jobs and files into children array
      const transformedChildren: (TreeViewItem | DjangoFile)[] = [
        ...node.children.map(transformNode),
        ...node.files,
      ];

      return {
        ...node,
        children: transformedChildren,
      };
    };

    return filtered.map(transformNode);
  }, [jobTree, filterText]);
};

const useItemLabel = (taskShortTitles: Map<string, string>) => {
  return useCallback((jobOrFile: TreeViewItem | DjangoFile): string => {
    const isJob = "parent" in jobOrFile;

    if (isJob) {
      const job = jobOrFile as TreeViewItem;
      const title = job.title?.trim()
        ? job.title
        : taskShortTitles.get(job.task_name) || job.task_name;
      return `${job.number}: ${title}`;
    }

    const file = jobOrFile as DjangoFile;
    return file.annotation.trim().length > 0
      ? file.annotation
      : file.job_param_name;
  }, [taskShortTitles]);
};

const useTreeItemData = (itemId: string): TreeItemData => {
  const { jobsByUuid, filesByUuid, taskShortTitles } = useContext(JobTreeContext);

  return useMemo(() => {
    const job = jobsByUuid.get(itemId);
    const file = filesByUuid.get(itemId);
    const isJob = Boolean(job);

    let displayLabel = "";
    let timestamp: string | undefined;

    if (job) {
      const title = job.title?.trim()
        ? job.title
        : taskShortTitles.get(job.task_name) || job.task_name;
      displayLabel = `${job.number}: ${title}`;
      timestamp =
        job.finish_time && new Date(job.finish_time).getFullYear() > 1970
          ? `Finished ${formatFriendlyDate(job.finish_time)}`
          : ACTIVE_JOB_STATUSES.includes(job.status) && job.creation_time
            ? `Started ${formatFriendlyDate(job.creation_time)}`
            : undefined;
    } else if (file) {
      displayLabel =
        file.annotation.trim().length > 0
          ? file.annotation
          : file.job_param_name;
    }

    return {
      job,
      file,
      isJob,
      displayLabel,
      timestamp,
    };
  }, [itemId, jobsByUuid, filesByUuid]);
};

// =============================================================================
// Utility functions
// =============================================================================

/**
 * Format a date string as a friendly relative/short form:
 * - Today: "14:32"
 * - Yesterday: "Yesterday 14:32"
 * - This year: "12 Mar 14:32"
 * - Older: "12 Mar 2024"
 */
const formatFriendlyDate = (isoString: string): string => {
  const date = new Date(isoString);
  const now = new Date();
  const time = date.toLocaleTimeString(undefined, { hour: "2-digit", minute: "2-digit" });

  const isToday =
    date.getDate() === now.getDate() &&
    date.getMonth() === now.getMonth() &&
    date.getFullYear() === now.getFullYear();
  if (isToday) return time;

  const yesterday = new Date(now);
  yesterday.setDate(yesterday.getDate() - 1);
  const isYesterday =
    date.getDate() === yesterday.getDate() &&
    date.getMonth() === yesterday.getMonth() &&
    date.getFullYear() === yesterday.getFullYear();
  if (isYesterday) return `Yesterday ${time}`;

  const dayMonth = date.toLocaleDateString(undefined, { day: "numeric", month: "short" });
  if (date.getFullYear() === now.getFullYear()) return `${dayMonth} ${time}`;

  return `${dayMonth} ${date.getFullYear()}`;
};

const formatFloatValue = (value: number): string => {
  return value.toPrecision(3);
};

// =============================================================================
// Main component
// =============================================================================

export const ClassicJobList: React.FC<ClassicJobListProps> = ({
  projectId,
  parent = null,
  withSubtitles = false,
}) => {
  const [selectedItems, setSelectedItems] = useState<string | null>(null);
  const navigate = useRouter();
  const api = useApi();
  const deleteDialog = useDeleteDialog();

  // Search/filter
  const [filterText, setFilterText] = useState<string | null>(null);

  // Multi-select mode for bulk delete
  const [selectMode, setSelectMode] = useState(false);
  const selectedJobIds = useSet<number>();

  const exitSelectMode = useCallback(() => {
    setSelectMode(false);
    selectedJobIds.clear();
  }, [selectedJobIds]);

  const toggleJobSelection = useCallback(
    (jobId: number) => {
      if (selectedJobIds.has(jobId)) {
        selectedJobIds.delete(jobId);
      } else {
        selectedJobIds.add(jobId);
      }
    },
    [selectedJobIds]
  );

  // Single consolidated API call
  const { jobTree, isLoading, mutate: mutateJobTree } = useJobTree(projectId);

  // Fetch task tree for shortTitle fallback (SWR deduplicates if already fetched)
  const { data: taskTreeResult } = api.get<any>(`task_tree/`);
  const taskShortTitles = useMemo(() => {
    const map = new Map<string, string>();
    const raw = taskTreeResult?.success
      ? taskTreeResult?.data?.task_tree
      : taskTreeResult?.task_tree;
    if (raw?.lookup) {
      for (const [taskName, meta] of Object.entries<any>(raw.lookup)) {
        if (meta?.shortTitle) {
          map.set(taskName, meta.shortTitle);
        }
      }
    }
    return map;
  }, [taskTreeResult]);

  // Build lookup maps for tree items
  const lookups = useJobTreeLookups(jobTree);

  // Context value including selection state
  const contextValue = useMemo<JobTreeContextValue>(
    () => ({
      ...lookups,
      taskShortTitles,
      selectMode,
      selectedJobIds,
      toggleJobSelection,
    }),
    [lookups, taskShortTitles, selectMode, selectedJobIds, toggleJobSelection]
  );

  // Transform to tree view format (with optional filtering)
  const treeViewItems = useTreeViewItems(jobTree, filterText);

  const getItemLabel = useItemLabel(taskShortTitles);

  const handleSelectedItemsChange = useCallback(
    (event: React.SyntheticEvent, ids: string | null) => {
      if (!ids) return;

      // In select mode, clicking a top-level job toggles its selection
      if (selectMode) {
        const job = lookups.jobsByUuid.get(ids);
        if (job && !job.number.includes(".")) {
          toggleJobSelection(job.id);
        }
        return;
      }

      const job = lookups.jobsByUuid.get(ids);
      if (job) {
        navigate.push(`/ccp4i2/project/${job.project}/job/${job.id}`);
      }
      setSelectedItems(ids);
    },
    [lookups.jobsByUuid, navigate, selectMode, toggleJobSelection]
  );

  const handleTreeSelection = useCallback(
    (event: React.SyntheticEvent, ids: string | null) => {
      // Prevent selection when clicking the expand/collapse icon
      const isIconClick = (event.target as Element).closest(
        '[class*="-MuiTreeItem2-iconContainer"]'
      );

      if (event.type !== "click" || !isIconClick) {
        handleSelectedItemsChange(event, ids);
      }
    },
    [handleSelectedItemsChange]
  );

  // Bulk delete handler
  const handleBulkDelete = useCallback(async () => {
    const jobIds = Array.from(selectedJobIds);
    if (jobIds.length === 0) return;

    try {
      const response: any = await api.post("jobs/bulk_dependent_jobs/", {
        job_ids: jobIds,
      });
      const { additional_dependents, total_to_delete, has_active_dependents } =
        response.data;

      // Filter to top-level dependents for display
      const topLevelDependents = (additional_dependents || []).filter(
        (job: Job) => job.parent === null
      );

      if (deleteDialog) {
        deleteDialog({
          type: "show",
          what: `${jobIds.length} selected job${jobIds.length !== 1 ? "s" : ""}`,
          onDelete: async () => {
            await api.post("jobs/bulk_delete/", { job_ids: jobIds });
            mutateJobTree();
            exitSelectMode();
          },
          onCancel: () => {},
          children:
            topLevelDependents.length > 0
              ? [
                  <Paper
                    key="dependentJobs"
                    sx={{ maxHeight: "10rem", overflowY: "auto" }}
                  >
                    <Typography variant="body2" sx={{ mb: 1 }}>
                      The following {topLevelDependents.length} dependent job
                      {topLevelDependents.length !== 1 ? "s" : ""} would also be
                      deleted:
                    </Typography>
                    <List dense>
                      {topLevelDependents.map((dependentJob: Job) => (
                        <ListItem key={dependentJob.uuid}>
                          <Toolbar>
                            <CCP4i2JobAvatar job={dependentJob} />
                            {`${dependentJob.number}: ${dependentJob.title}`}
                          </Toolbar>
                        </ListItem>
                      ))}
                    </List>
                  </Paper>,
                ]
              : undefined,
          deleteDisabled: has_active_dependents,
        });
      }
    } catch (error) {
      console.error("Failed to fetch bulk dependencies:", error);
    }
  }, [selectedJobIds, api, deleteDialog, mutateJobTree]);

  if (isLoading) {
    return <Skeleton variant="rectangular" width="100%" height={200} />;
  }

  if (!jobTree?.length) {
    return null;
  }

  return (
    <JobTreeContext.Provider value={contextValue}>
      {selectMode ? (
        <Paper
          elevation={2}
          sx={{ p: 1, mb: 1, bgcolor: "action.selected" }}
        >
          <Stack direction="row" alignItems="center" spacing={1}>
            <Button
              size="small"
              onClick={exitSelectMode}
              startIcon={<Clear fontSize="small" />}
            >
              Cancel
            </Button>
            <Typography
              color="primary.main"
              variant="subtitle2"
              sx={{ fontWeight: 600, flexGrow: 1 }}
            >
              {selectedJobIds.size > 0
                ? `${selectedJobIds.size} job${selectedJobIds.size !== 1 ? "s" : ""} selected`
                : "Select jobs to delete"}
            </Typography>
            <Tooltip title="Delete selected jobs">
              <span>
                <IconButton
                  onClick={handleBulkDelete}
                  size="small"
                  color="error"
                  disabled={selectedJobIds.size === 0}
                >
                  <Delete fontSize="small" />
                </IconButton>
              </span>
            </Tooltip>
          </Stack>
        </Paper>
      ) : (
        <Stack direction="row" alignItems="center" spacing={2} sx={{ p: 2 }}>
          <SearchField
            value={filterText || ""}
            onChange={setFilterText}
            placeholder="Search jobs…"
            size="small"
          />
          <Tooltip title="Select jobs to delete">
            <IconButton
              size="small"
              onClick={() => setSelectMode(true)}
            >
              <CheckBoxOutlined fontSize="small" />
            </IconButton>
          </Tooltip>
        </Stack>
      )}
      <RichTreeView
        items={treeViewItems}
        isItemEditable={() => true}
        experimentalFeatures={{ labelEditing: true }}
        getItemId={(jobOrFile) => jobOrFile.uuid}
        getItemLabel={getItemLabel}
        slots={{ item: CustomTreeItem }}
        onSelectedItemsChange={handleTreeSelection}
        selectedItems={selectedItems}
        sx={{ flex: "auto", overflowY: "auto", scrollbarWidth: "thin" }}
      />
    </JobTreeContext.Provider>
  );
};

// =============================================================================
// Custom Tree Item Component
// =============================================================================

const CustomTreeItem = forwardRef<HTMLLIElement, TreeItem2Props>(
  function CustomTreeItem({ id, itemId, label, disabled, children }, ref) {
    const { job, file, isJob, timestamp } = useTreeItemData(itemId);
    const { selectMode, selectedJobIds, toggleJobSelection } = useContext(JobTreeContext);

    const { setJobMenuAnchorEl, setJob } = useJobMenu();
    const { setFileMenuAnchorEl, setFile } = useFileMenu();

    // @dnd-kit drag — used for jobs only (context setting)
    const { attributes, listeners, setNodeRef } = useDraggable({
      id: job ? `job_${itemId}` : `file_${itemId}`,
      data: { job, file },
      disabled: !!file, // Disable @dnd-kit for files — they use native HTML5 drag
    });

    // Native HTML5 drag for files — enables OS drag-out to Finder/Explorer
    const handleFileDragStart = useCallback(
      (e: React.DragEvent) => {
        if (!file) return;

        // Set ccp4i2 file reference for within-window drops
        const ref = {
          ccp4i2_file: true,
          uuid: file.uuid,
          id: file.id,
          name: file.name,
          type: file.type,
          sub_type: file.sub_type,
          content: file.content,
          annotation: file.annotation,
          job: file.job,
          job_param_name: file.job_param_name,
        };
        e.dataTransfer.setData("application/ccp4i2-file", JSON.stringify(ref));
        e.dataTransfer.effectAllowed = "copy";

        // Write to clipboard for cross-window paste
        navigator.clipboard.writeText(JSON.stringify(ref)).catch(() => {});

      },
      [file]
    );

    // Tree item hooks
    const {
      getRootProps,
      getContentProps,
      getLabelProps,
      getGroupTransitionProps,
      getIconContainerProps,
      getLabelInputProps,
      status,
    } = useTreeItem2({ id, itemId, label, disabled, children, rootRef: ref });

    // KPI content for jobs - now uses embedded data, no filtering needed
    const kpiContent = useMemo(() => {
      if (!job) return null;

      const { kpis } = job;
      if (!kpis) return null;

      const charChips = Object.entries(kpis.char_values || {}).map(
        ([key, value]) => (
          <Chip key={`char_${key}`} label={`${key}: ${value}`} size="small" />
        )
      );

      const floatChips = Object.entries(kpis.float_values || {}).map(
        ([key, value]) => (
          <Chip
            key={`float_${key}`}
            label={`${key}: ${formatFloatValue(value)}`}
            size="small"
          />
        )
      );

      if (charChips.length === 0 && floatChips.length === 0) {
        return null;
      }

      return (
        <Stack direction="row" spacing={0.5} flexWrap="wrap">
          {charChips}
          {floatChips}
        </Stack>
      );
    }, [job]);

    const handleMenuClick = useCallback(
      (event: React.MouseEvent<HTMLElement>) => {
        event.stopPropagation();
        event.preventDefault(); // Prevent default context menu

        if (job) {
          setJobMenuAnchorEl(event.currentTarget);
          // Cast to Job since setJob expects Job type but we have JobTreeNode
          // They share the same properties except files array type
          setJob(job as unknown as Job);
        } else if (file) {
          setFileMenuAnchorEl(event.currentTarget);
          setFile(file);
        }
      },
      [job, file, setJobMenuAnchorEl, setJob, setFileMenuAnchorEl, setFile]
    );

    // Handle right-click context menu
    const handleContextMenu = useCallback(
      (event: React.MouseEvent<HTMLElement>) => {
        handleMenuClick(event);
      },
      [handleMenuClick]
    );

    // Handle button click (left-click on menu button)
    const handleButtonClick = useCallback(
      (event: React.MouseEvent<HTMLButtonElement>) => {
        handleMenuClick(event);
      },
      [handleMenuClick]
    );

    // Double click handler for jobs
    const handleDoubleClick = useCallback(
      (event: React.MouseEvent<HTMLElement>) => {
        if (job) {
          const path = `/job/${job.id}`;
          window.open(path, "_blank", "noopener,noreferrer");
        }
      },
      [job]
    );

    const renderAvatar = () => {
      if (job) {
        return (
          <CCP4i2JobAvatar
            job={job as unknown as Job}
            ref={setNodeRef}
            {...listeners}
            {...attributes}
          />
        );
      }

      if (file) {
        return (
          <div
            draggable
            onDragStart={handleFileDragStart}
            style={{ display: "inline-flex", cursor: "grab" }}
          >
            <FileAvatar
              file={file}
              sx={{ pointerEvents: "none" }}
            />
          </div>
        );
      }

      return null;
    };

    const labelStyle = useMemo(() => {
      if (!isJob) {
        // Files: smaller font, italic
        return { fontSize: "0.8rem", fontStyle: "italic", fontWeight: 400 };
      }
      if (job && !job.number.includes(".")) {
        // Top-level jobs: bold
        return { fontWeight: 700 };
      }
      // Subjobs: smaller, non-bold
      return { fontSize: "0.85rem", fontWeight: 400 };
    }, [isJob, job]);

    const renderContent = () => {
      if (status.editing) {
        return <TreeItem2LabelInput {...getLabelInputProps()} />;
      }

      return (
        <Stack direction="column" sx={{ flexGrow: 1 }}>
          <TreeItem2Label {...getLabelProps()} sx={labelStyle} />
          {timestamp && (
            <Typography
              variant="body2"
              sx={TIME_DISPLAY_STYLE}
              color="text.secondary"
              noWrap
            >
              {timestamp}
            </Typography>
          )}
          {kpiContent}
        </Stack>
      );
    };

    return (
      <TreeItem2Root
        {...getRootProps()}
        sx={undefined}
      >
        <TreeItem2Content
          {...getContentProps()}
          onContextMenu={handleContextMenu}
          onDoubleClick={handleDoubleClick}
          sx={{
            cursor: "context-menu",
            "&:hover": {
              backgroundColor: "action.hover",
            },
          }}
        >
          <TreeItem2IconContainer
            {...getIconContainerProps()}
            sx={{
              width: 24,
              height: 24,
              minWidth: 24,
              minHeight: 24,
              display: "flex",
              alignItems: "center",
              justifyContent: "center",
              borderRadius: 4,
              backgroundColor: children ? "rgba(0,0,0,0.04)" : "transparent",
              boxSizing: "border-box",
              mr: 1,
            }}
          >
            <TreeItem2Icon status={status} />
          </TreeItem2IconContainer>

          {selectMode && job && !job.number.includes(".") && (
            <Checkbox
              size="small"
              checked={selectedJobIds.has(job.id)}
              onChange={(e) => {
                e.stopPropagation();
                toggleJobSelection(job.id);
              }}
              onClick={(e) => e.stopPropagation()}
              sx={{ p: 0.5 }}
            />
          )}

          <Stack
            direction="row"
            alignItems="center"
            spacing={1}
            sx={{ flexGrow: 1 }}
          >
            {renderAvatar()}
            {renderContent()}
          </Stack>

          <Button
            size="small"
            variant="outlined"
            sx={{
              p: 0.5,
              minWidth: "auto",
              ml: 1,
            }}
            onClick={handleButtonClick}
            aria-label={`Open ${isJob ? "job" : "file"} menu`}
          >
            <MenuIcon fontSize="small" />
          </Button>
        </TreeItem2Content>

        {children && (
          <TreeItem2GroupTransition {...getGroupTransitionProps()} />
        )}
      </TreeItem2Root>
    );
  }
);

// Set display name for debugging
CustomTreeItem.displayName = "CustomTreeItem";
