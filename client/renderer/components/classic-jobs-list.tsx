import React, {
  forwardRef,
  useCallback,
  useContext,
  useMemo,
  useState,
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
  JobCharValue,
  JobFloatValue,
} from "../types/models";
import { CCP4i2JobAvatar } from "./job-avatar";
import { FileAvatar } from "./file-avatar";
import { useCCP4i2Window } from "../app-context";
import { JobWithChildren, useJobMenu } from "../providers/job-context-menu";
import { useFileMenu } from "../providers/file-context-menu";

// Types
interface ClassicJobListProps {
  projectId: number;
  parent?: number;
  withSubtitles?: boolean;
}

interface TreeItemData {
  job?: Job;
  file?: DjangoFile;
  isJob: boolean;
  displayLabel: string;
  timestamp?: string;
}

// Constants
const CACHE_TIME = 10000;
const TREE_ITEM_BORDER_STYLE = { border: "1px solid #999" };
const TIME_DISPLAY_STYLE = { fontSize: "75%" };

// Custom hooks
const useProjectData = (projectId: number) => {
  const api = useApi();

  const endpointFetch: EndpointFetch = {
    type: "projects",
    id: projectId,
    endpoint: "jobs",
  };

  const { data: jobs } = api.get_endpoint<Job[]>(endpointFetch, CACHE_TIME);
  const { data: files } = api.get_endpoint<DjangoFile[]>(
    {
      type: "projects",
      id: projectId,
      endpoint: "files",
    },
    CACHE_TIME
  );
  const { data: jobCharValues } = api.get_endpoint<JobCharValue[]>(
    {
      type: "projects",
      id: projectId,
      endpoint: "job_char_values/",
    },
    CACHE_TIME
  );
  const { data: jobFloatValues } = api.get_endpoint<JobFloatValue[]>(
    {
      type: "projects",
      id: projectId,
      endpoint: "job_float_values/",
    },
    CACHE_TIME
  );

  return {
    jobs,
    files,
    jobCharValues,
    jobFloatValues,
    isLoading: !jobs || !files,
  };
};

const useDecoratedJobs = (
  jobs: Job[] | undefined,
  files: DjangoFile[] | undefined,
  parent: number | null
) => {
  return useMemo<(JobWithChildren | DjangoFile)[] | undefined>(() => {
    if (!jobs?.filter) return [];

    // Helper function to parse job number into comparable array
    const parseJobNumber = (jobNumber: string): number[] => {
      return jobNumber.split(".").map((num) => parseInt(num, 10) || 0);
    };

    // Helper function to compare job numbers (highest ordinal first)
    const compareJobNumbers = (a: string, b: string): number => {
      const aParts = parseJobNumber(a);
      const bParts = parseJobNumber(b);

      // Compare each part
      for (let i = 0; i < Math.max(aParts.length, bParts.length); i++) {
        const aPart = aParts[i] || 0;
        const bPart = bParts[i] || 0;

        if (aPart !== bPart) {
          // Higher numbers first (descending order)
          return bPart - aPart;
        }
      }

      // If all parts are equal, maintain original order
      return 0;
    };

    return jobs
      .filter((job) => job.parent === parent)
      .map((job) => {
        const childJobs: (Job | DjangoFile)[] = jobs.filter(
          (childJob) => childJob.parent === job.id
        );
        const childFiles = files?.filter((file) => file.job === job.id) || [];
        return {
          ...job,
          children: [...childJobs, ...childFiles],
        };
      })
      .sort((a, b) => compareJobNumbers(a.number, b.number));
  }, [jobs, files, parent]);
};

const useItemLabel = () => {
  return useCallback((jobOrFile: JobWithChildren | DjangoFile): string => {
    const isJob = "parent" in jobOrFile;

    if (isJob) {
      const job = jobOrFile as Job;
      return `${job.number}: ${job.title}`;
    }

    const file = jobOrFile as DjangoFile;
    return file.annotation.trim().length > 0
      ? file.annotation
      : file.job_param_name;
  }, []);
};

const useTreeItemData = (
  itemId: string,
  jobs?: Job[],
  files?: DjangoFile[]
): TreeItemData => {
  return useMemo(() => {
    const job = jobs?.find((j) => j.uuid === itemId);
    const file = files?.find((f) => f.uuid === itemId);
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
  }, [itemId, jobs, files]);
};

// Utility functions
const shouldShowTreeItemBorder = (job?: Job): boolean => {
  return Boolean(job && !job.number.includes("."));
};

const formatFloatValue = (value: number): string => {
  return value.toPrecision(3);
};

// Main component
export const ClassicJobList: React.FC<ClassicJobListProps> = ({
  projectId,
  parent = null,
  withSubtitles = false,
}) => {
  const [selectedItems, setSelectedItems] = useState<string | null>(null);
  const navigate = useRouter();

  const { jobs, files, isLoading } = useProjectData(projectId);
  const decoratedJobs = useDecoratedJobs(jobs, files, parent);
  const getItemLabel = useItemLabel();

  const handleSelectedItemsChange = useCallback(
    (event: React.SyntheticEvent, ids: string | null) => {
      if (!ids || !jobs) return;

      const job = jobs.find((job) => job.uuid === ids);
      if (job) {
        navigate.push(`/project/${job.project}/job/${job.id}`);
      }
      setSelectedItems(ids);
    },
    [jobs, navigate]
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

  if (!decoratedJobs) {
    return null;
  }

  return (
    <RichTreeView
      items={decoratedJobs}
      isItemEditable={() => true}
      experimentalFeatures={{ labelEditing: true }}
      getItemId={(jobOrFile) => jobOrFile.uuid}
      getItemLabel={getItemLabel}
      slots={{ item: CustomTreeItem }}
      onSelectedItemsChange={handleTreeSelection}
      selectedItems={selectedItems}
    />
  );
};

// Custom Tree Item Component
const CustomTreeItem = forwardRef<HTMLLIElement, TreeItem2Props>(
  function CustomTreeItem({ id, itemId, label, disabled, children }, ref) {
    const { projectId } = useCCP4i2Window();
    const { jobs, files, jobCharValues, jobFloatValues } = useProjectData(
      projectId ?? 0
    );
    const { job, file, isJob, timestamp } = useTreeItemData(
      itemId,
      jobs,
      files
    );

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

    // KPI content for jobs
    const kpiContent = useMemo(() => {
      if (!job || !jobCharValues || !jobFloatValues) return null;

      const charValues = jobCharValues
        .filter((item) => item.job === job.id)
        .map((item) => (
          <Chip
            key={item.key}
            label={`${item.key}: ${item.value}`}
            size="small"
          />
        ));

      const floatValues = jobFloatValues
        .filter((item) => item.job === job.id)
        .map((item) => (
          <Chip
            key={item.key}
            label={`${item.key}: ${formatFloatValue(item.value)}`}
            size="small"
          />
        ));

      return (
        <Stack direction="row" spacing={0.5} flexWrap="wrap">
          {charValues}
          {floatValues}
        </Stack>
      );
    }, [job, jobCharValues, jobFloatValues]);

    const handleMenuClick = useCallback(
      (event: React.MouseEvent<HTMLElement>) => {
        event.stopPropagation();
        event.preventDefault(); // Prevent default context menu

        if (job) {
          setJobMenuAnchorEl(event.currentTarget);
          setJob(job);
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
            job={job}
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
              backgroundColor: children ? "rgba(0,0,0,0.04)" : "transparent", // Only show background if there are children
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
