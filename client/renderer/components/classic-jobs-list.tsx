import React, {
  forwardRef,
  useCallback,
  useMemo,
  useState,
  createContext,
  useContext,
} from "react";
import { Button, Chip, Skeleton, Stack, Typography } from "@mui/material";
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
import { Menu as MenuIcon } from "@mui/icons-material";
import { useDraggable } from "@dnd-kit/core";
import { useRouter } from "next/navigation";

import { EndpointFetch, useApi } from "../api";
import {
  Job,
  File as DjangoFile,
} from "../types/models";
import { CCP4i2JobAvatar } from "./job-avatar";
import { FileAvatar } from "./file-avatar";
import { useCCP4i2Window } from "../app-context";
import { JobWithChildren, useJobMenu } from "../providers/job-context-menu";
import { useFileMenu } from "../providers/file-context-menu";
import { useRecentlyStartedJobs } from "../providers/recently-started-jobs-context";

// =============================================================================
// Types
// =============================================================================

interface ClassicJobListProps {
  projectId: number;
  parent?: number;
  withSubtitles?: boolean;
}

/** Job node with embedded files and KPIs from job_tree endpoint
 * We use Omit to override the files property type (Job has number[] for IDs,
 * JobTreeNode has full DjangoFile objects).
 */
interface JobTreeNode extends Omit<Job, 'files'> {
  files: DjangoFile[];
  kpis: {
    float_values: Record<string, number>;
    char_values: Record<string, string>;
  };
  children: JobTreeNode[];
}

/** Response from job_tree endpoint */
interface JobTreeResponse {
  job_tree: JobTreeNode[];
  total_jobs: number;
  total_files: number;
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

const TREE_ITEM_BORDER_STYLE = { border: "1px solid #999" };
const TIME_DISPLAY_STYLE = { fontSize: "75%" };

// =============================================================================
// Context for sharing job tree data with tree items
// =============================================================================

interface JobTreeContextValue {
  jobsByUuid: Map<string, JobTreeNode>;
  filesByUuid: Map<string, DjangoFile>;
}

const JobTreeContext = createContext<JobTreeContextValue>({
  jobsByUuid: new Map(),
  filesByUuid: new Map(),
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
 * Replaces 4 separate calls (jobs, files, job_char_values, job_float_values)
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

  const { data, isLoading, error } = api.get_endpoint<JobTreeResponse>(
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
    const filesByUuid = new Map<string, DjangoFile>();

    if (!jobTree) {
      return { jobsByUuid, filesByUuid };
    }

    // Recursively index all jobs and files
    const indexNode = (node: JobTreeNode) => {
      jobsByUuid.set(node.uuid, node);
      node.files.forEach((file) => filesByUuid.set(file.uuid, file));
      node.children.forEach(indexNode);
    };

    jobTree.forEach(indexNode);

    return { jobsByUuid, filesByUuid };
  }, [jobTree]);
};

/**
 * Transform job tree to format expected by RichTreeView.
 * Merges files into children array for each job.
 */
const useTreeViewItems = (jobTree: JobTreeNode[] | undefined): TreeViewItem[] => {
  return useMemo(() => {
    if (!jobTree) return [];

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

    return jobTree.map(transformNode);
  }, [jobTree]);
};

const useItemLabel = () => {
  return useCallback((jobOrFile: TreeViewItem | DjangoFile): string => {
    const isJob = "parent" in jobOrFile;

    if (isJob) {
      const job = jobOrFile as TreeViewItem;
      return `${job.number}: ${job.title}`;
    }

    const file = jobOrFile as DjangoFile;
    return file.annotation.trim().length > 0
      ? file.annotation
      : file.job_param_name;
  }, []);
};

const useTreeItemData = (itemId: string): TreeItemData => {
  const { jobsByUuid, filesByUuid } = useContext(JobTreeContext);

  return useMemo(() => {
    const job = jobsByUuid.get(itemId);
    const file = filesByUuid.get(itemId);
    const isJob = Boolean(job);

    let displayLabel = "";
    let timestamp: string | undefined;

    if (job) {
      displayLabel = `${job.number}: ${job.title}`;
      timestamp =
        job.finish_time && new Date(job.finish_time).getFullYear() > 1970
          ? `Finished ${new Date(job.finish_time).toLocaleString()}`
          : job.creation_time
            ? `Modified ${new Date(job.creation_time).toLocaleString()}`
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

const shouldShowTreeItemBorder = (job?: JobTreeNode): boolean => {
  return Boolean(job && !job.number.includes("."));
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

  // Single consolidated API call
  const { jobTree, isLoading } = useJobTree(projectId);

  // Build lookup maps for tree items
  const lookups = useJobTreeLookups(jobTree);

  // Transform to tree view format
  const treeViewItems = useTreeViewItems(jobTree);

  const getItemLabel = useItemLabel();

  const handleSelectedItemsChange = useCallback(
    (event: React.SyntheticEvent, ids: string | null) => {
      if (!ids) return;

      const job = lookups.jobsByUuid.get(ids);
      if (job) {
        navigate.push(`/project/${job.project}/job/${job.id}`);
      }
      setSelectedItems(ids);
    },
    [lookups.jobsByUuid, navigate]
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

  if (isLoading) {
    return <Skeleton variant="rectangular" width="100%" height={200} />;
  }

  if (!treeViewItems.length) {
    return null;
  }

  return (
    <JobTreeContext.Provider value={lookups}>
      <RichTreeView
        items={treeViewItems}
        isItemEditable={() => true}
        experimentalFeatures={{ labelEditing: true }}
        getItemId={(jobOrFile) => jobOrFile.uuid}
        getItemLabel={getItemLabel}
        slots={{ item: CustomTreeItem }}
        onSelectedItemsChange={handleTreeSelection}
        selectedItems={selectedItems}
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

    const { setJobMenuAnchorEl, setJob } = useJobMenu();
    const { setFileMenuAnchorEl, setFile } = useFileMenu();

    // Drag and drop setup
    const { attributes, listeners, setNodeRef } = useDraggable({
      id: job ? `job_${itemId}` : `file_${itemId}`,
      data: { job, file },
    });

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

      // Debug: log KPI data with actual counts
      const floatCount = Object.keys(kpis?.float_values || {}).length;
      const charCount = Object.keys(kpis?.char_values || {}).length;
      if (floatCount > 0 || charCount > 0) {
        console.log(`Job ${job.number} has ${floatCount} float KPIs, ${charCount} char KPIs:`, kpis);
      }

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
          <FileAvatar
            file={file}
            ref={setNodeRef}
            {...listeners}
            {...attributes}
          />
        );
      }

      return null;
    };

    const renderContent = () => {
      if (status.editing) {
        return <TreeItem2LabelInput {...getLabelInputProps()} />;
      }

      return (
        <Stack direction="column" sx={{ flexGrow: 1 }}>
          <TreeItem2Label {...getLabelProps()} />
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
        sx={shouldShowTreeItemBorder(job) ? TREE_ITEM_BORDER_STYLE : undefined}
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
