"use client";

import { useCallback, useMemo, useRef, useState } from "react";
import { useRouter } from "next/navigation";
import {
  Avatar,
  Box,
  Chip,
  IconButton,
  ListItemIcon,
  ListItemText,
  Menu,
  MenuItem,
  Skeleton,
  Stack,
  Table,
  TableBody,
  TableCell,
  TableContainer,
  TableHead,
  TableRow,
  Tooltip,
  Typography,
} from "@mui/material";
import {
  MoreHoriz as MoreIcon,
  Delete as DeleteIcon,
  PlayArrow as RerunIcon,
} from "@mui/icons-material";
import { useVirtualizer } from "@tanstack/react-virtual";
import { apiPost } from "../../api-fetch";
import {
  MemberProjectWithSummary,
  CampaignJobInfo,
  parseDatasetFilename,
} from "../../types/campaigns";
import { SmilesView } from "./smiles-view";
import { ProjectTagChips } from "../project-tag-chips";

// Status ID to color mapping (matching legacy CCP4i2)
const STATUS_COLORS: Record<number, string> = {
  0: "#AAA", // Unknown
  1: "#FFF", // Created
  2: "#FFA", // Pending
  3: "#AAF", // Running
  4: "#FDA", // Interrupted
  5: "#FAA", // Failed
  6: "#AFA", // Finished
};

interface VirtualizedMemberProjectsTableProps {
  projects: MemberProjectWithSummary[];
  smilesMap: Record<number, string>;
  showSubJobs: boolean;
  latestCoordsFileId?: number;
  onRefresh: () => void;
  onDelete: (project: MemberProjectWithSummary) => void;
  onProjectClick: (project: MemberProjectWithSummary) => void;
  /** Maximum height of the table container (default: 500) */
  maxHeight?: number;
}

export function VirtualizedMemberProjectsTable({
  projects,
  smilesMap,
  showSubJobs,
  latestCoordsFileId,
  onRefresh,
  onDelete,
  onProjectClick,
  maxHeight = 500,
}: VirtualizedMemberProjectsTableProps) {
  const parentRef = useRef<HTMLDivElement>(null);

  // Set up virtualizer for windowed rendering
  const rowVirtualizer = useVirtualizer({
    count: projects.length,
    getScrollElement: () => parentRef.current,
    estimateSize: () => 85, // Estimated row height (accounts for SMILES and tags)
    overscan: 5, // Render 5 extra rows above/below viewport
  });

  const virtualItems = rowVirtualizer.getVirtualItems();

  if (projects.length === 0) {
    return (
      <Box sx={{ textAlign: "center", py: 4 }}>
        <Typography color="text.secondary">
          No member projects yet. Use Batch Import to add datasets.
        </Typography>
      </Box>
    );
  }

  return (
    <TableContainer
      ref={parentRef}
      sx={{
        maxHeight,
        overflow: "auto",
      }}
    >
      <Table stickyHeader size="small" sx={{ tableLayout: "fixed" }}>
        <TableHead>
          <TableRow>
            <TableCell width={60}>Actions</TableCell>
            <TableCell>Project Name</TableCell>
            <TableCell width={120}>Ligand</TableCell>
            <TableCell align="center" width={80}>Resolution</TableCell>
            <TableCell align="center" width={80}>R-Factor</TableCell>
            <TableCell align="center" width={80}>R-Free</TableCell>
            <TableCell>Jobs (click for CCP4i2, shift-click for Moorhen)</TableCell>
          </TableRow>
        </TableHead>
        <TableBody>
          {/* Spacer for rows above viewport */}
          {virtualItems.length > 0 && virtualItems[0].start > 0 && (
            <TableRow>
              <TableCell
                colSpan={7}
                sx={{
                  height: virtualItems[0].start,
                  padding: 0,
                  border: "none",
                }}
              />
            </TableRow>
          )}

          {/* Virtualized rows */}
          {virtualItems.map((virtualRow) => {
            const project = projects[virtualRow.index];
            return (
              <MemberProjectRow
                key={project.id}
                project={project}
                smilesMap={smilesMap}
                showSubJobs={showSubJobs}
                latestCoordsFileId={latestCoordsFileId}
                onRefresh={onRefresh}
                onDelete={() => onDelete(project)}
                onProjectClick={() => onProjectClick(project)}
                virtualIndex={virtualRow.index}
                measureElement={rowVirtualizer.measureElement}
              />
            );
          })}

          {/* Spacer for rows below viewport */}
          {virtualItems.length > 0 && (
            <TableRow>
              <TableCell
                colSpan={7}
                sx={{
                  height:
                    rowVirtualizer.getTotalSize() -
                    (virtualItems[virtualItems.length - 1]?.end ?? 0),
                  padding: 0,
                  border: "none",
                }}
              />
            </TableRow>
          )}
        </TableBody>
      </Table>
    </TableContainer>
  );
}

// Member project row component
interface MemberProjectRowProps {
  project: MemberProjectWithSummary;
  smilesMap: Record<number, string>;
  showSubJobs: boolean;
  latestCoordsFileId?: number;
  onRefresh: () => void;
  onDelete: () => void;
  onProjectClick: () => void;
  virtualIndex: number;
  measureElement: (element: HTMLElement | null) => void;
}

function MemberProjectRow({
  project,
  smilesMap,
  showSubJobs,
  latestCoordsFileId,
  onRefresh,
  onDelete,
  onProjectClick,
  virtualIndex,
  measureElement,
}: MemberProjectRowProps) {
  const router = useRouter();
  const [menuAnchor, setMenuAnchor] = useState<HTMLElement | null>(null);

  // Parse project name for metadata
  const parsed = parseDatasetFilename(project.name);
  const regId = parsed.nclId ? parseInt(parsed.nclId) : null;
  const smiles = regId && regId !== 0 ? smilesMap[regId] : null;

  // Filter jobs - show only top-level unless showSubJobs is true
  const visibleJobs = useMemo(() => {
    if (!project.jobs) return [];
    const filtered = showSubJobs
      ? project.jobs
      : project.jobs.filter((j) => !j.number.includes("."));
    // Sort by job number (natural order)
    return [...filtered].sort((a, b) => {
      const aParts = a.number.split(".").map(Number);
      const bParts = b.number.split(".").map(Number);
      for (let i = 0; i < Math.max(aParts.length, bParts.length); i++) {
        const aVal = aParts[i] || 0;
        const bVal = bParts[i] || 0;
        if (aVal !== bVal) return aVal - bVal;
      }
      return 0;
    });
  }, [project.jobs, showSubJobs]);

  // Handle job click - navigate to job view, shift-click for Moorhen
  const handleJobClick = useCallback(
    (job: CampaignJobInfo, event: React.MouseEvent) => {
      event.stopPropagation();
      event.preventDefault();
      if (event.shiftKey) {
        // Shift-click opens Moorhen in new tab
        window.open(`/ccp4i2/moorhen-page/job-by-id/${job.id}`, "_blank");
      } else {
        // Regular click opens job in project view
        router.push(`/ccp4i2/project/${project.id}/job/${job.id}`);
      }
    },
    [router, project.id]
  );

  // Handle rerun initial job
  const handleRerun = useCallback(async () => {
    setMenuAnchor(null);
    try {
      const slJob = project.jobs?.find((j) => j.task_name === "SubstituteLigand");
      if (!slJob) return;
      const cloneResult = await apiPost<{ new_job: { id: number } }>(
        `jobs/${slJob.id}/clone/`,
        {}
      );
      await apiPost(`jobs/${cloneResult.new_job.id}/run/`, { queue: "batch" });
      onRefresh();
    } catch (err) {
      console.error("Failed to rerun job:", err);
    }
  }, [project.jobs, onRefresh]);

  return (
    <TableRow
      hover
      data-index={virtualIndex}
      ref={(el) => measureElement(el)}
      sx={{ cursor: "pointer" }}
      onClick={onProjectClick}
    >
      {/* Actions menu */}
      <TableCell onClick={(e) => e.stopPropagation()}>
        <IconButton size="small" onClick={(e) => setMenuAnchor(e.currentTarget)}>
          <MoreIcon fontSize="small" />
        </IconButton>
        <Menu
          anchorEl={menuAnchor}
          open={Boolean(menuAnchor)}
          onClose={() => setMenuAnchor(null)}
        >
          <MenuItem
            onClick={() => {
              setMenuAnchor(null);
              onDelete();
            }}
          >
            <ListItemIcon>
              <DeleteIcon fontSize="small" />
            </ListItemIcon>
            <ListItemText>Delete project</ListItemText>
          </MenuItem>
          <MenuItem onClick={handleRerun}>
            <ListItemIcon>
              <RerunIcon fontSize="small" />
            </ListItemIcon>
            <ListItemText>Rerun initial job</ListItemText>
          </MenuItem>
        </Menu>
      </TableCell>

      {/* Project name with tags */}
      <TableCell>
        <Stack spacing={0.5}>
          <Typography variant="body2">{project.name}</Typography>
          <ProjectTagChips tags={project.tags} maxVisible={4} size="small" hideEmpty />
        </Stack>
      </TableCell>

      {/* Ligand */}
      <TableCell>
        {regId === 0 ? (
          <Typography variant="body2" color="text.secondary">
            Presumed Apo
          </Typography>
        ) : smiles ? (
          <Tooltip title={`NCL-${parsed.nclId} (reg#${regId})`} placement="bottom">
            <Box>
              <SmilesView smiles={smiles} width={100} height={75} />
            </Box>
          </Tooltip>
        ) : parsed.nclId ? (
          <Box
            sx={{
              display: "flex",
              flexDirection: "column",
              alignItems: "center",
            }}
          >
            <Skeleton variant="rectangular" width={100} height={55} />
            <Chip
              label={`NCL-${parsed.nclId}`}
              size="small"
              color="primary"
              variant="outlined"
              sx={{ mt: 0.5 }}
            />
          </Box>
        ) : (
          <Typography color="text.secondary">-</Typography>
        )}
      </TableCell>

      {/* Resolution */}
      <TableCell align="center">
        <Chip
          label={
            project.kpis.highResLimit !== undefined
              ? Number(project.kpis.highResLimit).toFixed(2)
              : "-"
          }
          size="small"
          variant="outlined"
        />
      </TableCell>

      {/* R-Factor */}
      <TableCell align="center">
        <Chip
          label={
            project.kpis.RFactor !== undefined
              ? Number(project.kpis.RFactor).toFixed(3)
              : "-"
          }
          size="small"
          variant="outlined"
        />
      </TableCell>

      {/* R-Free */}
      <TableCell align="center">
        <Chip
          label={
            project.kpis.RFree !== undefined
              ? Number(project.kpis.RFree).toFixed(3)
              : "-"
          }
          size="small"
          variant="outlined"
        />
      </TableCell>

      {/* Jobs - clickable icons matching legacy style */}
      <TableCell onClick={(e) => e.stopPropagation()}>
        <Stack direction="row" spacing={0.5} flexWrap="wrap" useFlexGap>
          {visibleJobs.map((job) => (
            <Tooltip
              key={job.id}
              title={
                <Box>
                  <Typography variant="body2">{job.task_name}</Typography>
                  <Typography variant="caption">Job {job.number}</Typography>
                </Box>
              }
            >
              <Box
                sx={{
                  display: "flex",
                  flexDirection: "column",
                  alignItems: "center",
                  cursor: "pointer",
                  "&:hover": { opacity: 0.8 },
                }}
                onClick={(e) => handleJobClick(job, e)}
              >
                <Avatar
                  sx={{
                    width: 32,
                    height: 32,
                    bgcolor: STATUS_COLORS[job.status] || "#AAA",
                    border: "1px solid rgba(0,0,0,0.1)",
                  }}
                  src={`/svgicons/${job.task_name}.svg`}
                  imgProps={{
                    onError: (e: React.SyntheticEvent<HTMLImageElement>) => {
                      e.currentTarget.src = `/qticons/${job.task_name}.png`;
                    },
                  }}
                >
                  {job.task_name?.[0]?.toUpperCase()}
                </Avatar>
                <Typography variant="caption" sx={{ fontSize: "0.65rem" }}>
                  {job.number}
                </Typography>
              </Box>
            </Tooltip>
          ))}
          {visibleJobs.length === 0 && (
            <Typography color="text.secondary" variant="body2">
              No jobs
            </Typography>
          )}
        </Stack>
      </TableCell>
    </TableRow>
  );
}
