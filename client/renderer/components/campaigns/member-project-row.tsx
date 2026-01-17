"use client";

import { useCallback, useState, useMemo } from "react";
import { useRouter } from "next/navigation";
import {
  Avatar,
  Box,
  Button,
  Chip,
  IconButton,
  ListItemIcon,
  ListItemText,
  Menu,
  MenuItem,
  Skeleton,
  Stack,
  TableCell,
  TableRow,
  Tooltip,
  Typography,
} from "@mui/material";
import {
  MoreHoriz as MoreIcon,
  Delete as DeleteIcon,
  PlayArrow as RerunIcon,
  Refresh as RefreshCoordsIcon,
} from "@mui/icons-material";
import { useApi } from "../../api";
import { apiPost, apiDelete } from "../../api-fetch";
import { parseDatasetFilename } from "../../types/campaigns";
import { Job } from "../../types/models";

// Status ID to color mapping (matching legacy)
const STATUS_COLORS: Record<number, string> = {
  0: "#AAA", // Unknown
  1: "#FFF", // Created
  2: "#FFA", // Pending
  3: "#AAF", // Running
  4: "#FDA", // Interrupted
  5: "#FAA", // Failed
  6: "#AFA", // Finished
};

interface MemberProjectRowProps {
  project: {
    id: number;
    name: string;
    uuid: string;
  };
  jobs: Job[];
  smilesMap: Record<number, string>;
  showSubJobs: boolean;
  latestCoordsFileId?: number;
  onRefresh: () => void;
  onDelete: (project: { id: number; name: string }) => void;
}

export function MemberProjectRow({
  project,
  jobs,
  smilesMap,
  showSubJobs,
  latestCoordsFileId,
  onRefresh,
  onDelete,
}: MemberProjectRowProps) {
  const router = useRouter();
  const api = useApi();
  const [menuAnchor, setMenuAnchor] = useState<HTMLElement | null>(null);

  // Parse project name for metadata
  const parsed = parseDatasetFilename(project.name);
  const regId = parsed.nclId ? parseInt(parsed.nclId) : null;
  const smiles = regId && regId !== 0 ? smilesMap[regId] : null;

  // Filter and sort jobs for this project
  const projectJobs = useMemo(() => {
    return jobs
      .filter((j) => j.project === project.id)
      .sort((a, b) => {
        // Sort by job number (highest first)
        const aParts = a.number.split(".").map(Number);
        const bParts = b.number.split(".").map(Number);
        for (let i = 0; i < Math.max(aParts.length, bParts.length); i++) {
          const aVal = aParts[i] || 0;
          const bVal = bParts[i] || 0;
          if (aVal !== bVal) return bVal - aVal;
        }
        return 0;
      });
  }, [jobs, project.id]);

  // Get visible jobs (filter sub-jobs if needed)
  const visibleJobs = useMemo(() => {
    if (showSubJobs) return projectJobs;
    return projectJobs.filter((j) => !j.number.includes("."));
  }, [projectJobs, showSubJobs]);

  // Extract KPIs from most recent finished job
  const kpis = useMemo(() => {
    // This will need to be fetched from job_tree endpoint which includes kpis
    // For now, return empty - we'll enhance the API to include this
    return {
      highResLimit: undefined as number | undefined,
      RFactor: undefined as number | undefined,
      RFree: undefined as number | undefined,
    };
  }, [projectJobs]);

  // Handle job click
  const handleJobClick = useCallback(
    (job: Job, event: React.MouseEvent) => {
      event.stopPropagation();
      if (event.shiftKey) {
        // Shift-click opens Moorhen in new tab (needs separate window for cross-origin isolation)
        window.open(`/ccp4i2/moorhen-page/job-by-id/${job.id}`, '_blank');
      } else {
        // Regular click opens job in project view
        router.push(`/ccp4i2/project/${project.id}/job/${job.id}`);
      }
    },
    [router, project.id]
  );

  // Handle rerun
  const handleRerun = useCallback(async () => {
    setMenuAnchor(null);
    try {
      // Find initial SubstituteLigand job
      const slJob = projectJobs.find((j) => j.task_name === "SubstituteLigand");
      if (!slJob) {
        console.error("No SubstituteLigand job found");
        return;
      }

      // Clone and run
      const cloneResult = await apiPost<{ new_job: { id: number } }>(
        `jobs/${slJob.id}/clone/`,
        {}
      );
      await apiPost(`jobs/${cloneResult.new_job.id}/run/`, {});
      onRefresh();
    } catch (err) {
      console.error("Failed to rerun job:", err);
    }
  }, [projectJobs, onRefresh]);

  // Handle rerun with new coords
  const handleRerunWithNewCoords = useCallback(async () => {
    setMenuAnchor(null);
    if (!latestCoordsFileId) return;

    try {
      const slJob = projectJobs.find((j) => j.task_name === "SubstituteLigand");
      if (!slJob) return;

      // Clone job
      const cloneResult = await apiPost<{ new_job: { id: number; uuid: string } }>(
        `jobs/${slJob.id}/clone/`,
        {}
      );

      // TODO: Upload new coords to the cloned job
      // This requires the file upload workflow

      await apiPost(`jobs/${cloneResult.new_job.id}/run/`, {});
      onRefresh();
    } catch (err) {
      console.error("Failed to rerun with new coords:", err);
    }
  }, [projectJobs, latestCoordsFileId, onRefresh]);

  const handleRowClick = () => {
    router.push(`/ccp4i2/project/${project.id}`);
  };

  return (
    <TableRow
      hover
      sx={{ cursor: "pointer" }}
      onClick={handleRowClick}
    >
      {/* Actions */}
      <TableCell onClick={(e) => e.stopPropagation()}>
        <IconButton
          size="small"
          onClick={(e) => setMenuAnchor(e.currentTarget)}
        >
          <MoreIcon />
        </IconButton>
        <Menu
          anchorEl={menuAnchor}
          open={Boolean(menuAnchor)}
          onClose={() => setMenuAnchor(null)}
        >
          <MenuItem onClick={() => { setMenuAnchor(null); onDelete(project); }}>
            <ListItemIcon><DeleteIcon fontSize="small" /></ListItemIcon>
            <ListItemText>Delete project</ListItemText>
          </MenuItem>
          <MenuItem onClick={handleRerun}>
            <ListItemIcon><RerunIcon fontSize="small" /></ListItemIcon>
            <ListItemText>Rerun initial job</ListItemText>
          </MenuItem>
          {latestCoordsFileId && (
            <MenuItem onClick={handleRerunWithNewCoords}>
              <ListItemIcon><RefreshCoordsIcon fontSize="small" /></ListItemIcon>
              <ListItemText>Rerun with new coords</ListItemText>
            </MenuItem>
          )}
        </Menu>
      </TableCell>

      {/* Project Name */}
      <TableCell>
        <Typography variant="body2">{project.name}</Typography>
      </TableCell>

      {/* Ligand (SMILES or Apo) */}
      <TableCell>
        {regId === 0 ? (
          <Typography variant="body2" color="text.secondary">
            Presumed Apo
          </Typography>
        ) : smiles ? (
          <SmilesViewer smiles={smiles} width={100} height={75} />
        ) : parsed.nclId ? (
          <Skeleton variant="rectangular" width={100} height={75} />
        ) : (
          <Typography color="text.secondary">-</Typography>
        )}
      </TableCell>

      {/* Resolution */}
      <TableCell align="center">
        {kpis.highResLimit !== undefined ? (
          <Chip label={kpis.highResLimit.toFixed(2)} size="small" />
        ) : (
          "-"
        )}
      </TableCell>

      {/* R-Factor */}
      <TableCell align="center">
        {kpis.RFactor !== undefined ? (
          <Chip label={kpis.RFactor.toFixed(3)} size="small" />
        ) : (
          "-"
        )}
      </TableCell>

      {/* R-Free */}
      <TableCell align="center">
        {kpis.RFree !== undefined ? (
          <Chip label={kpis.RFree.toFixed(3)} size="small" />
        ) : (
          "-"
        )}
      </TableCell>

      {/* Jobs - clickable icons */}
      <TableCell onClick={(e) => e.stopPropagation()}>
        <Stack direction="row" spacing={0.5} flexWrap="wrap" useFlexGap>
          {visibleJobs.map((job) => (
            <Tooltip
              key={job.id}
              title={
                <Box>
                  <Typography variant="body2">{job.task_name}</Typography>
                  <Typography variant="caption">
                    Job {job.number} - Click to view, Shift+Click for Moorhen
                  </Typography>
                </Box>
              }
            >
              <Box
                sx={{
                  display: "flex",
                  flexDirection: "column",
                  alignItems: "center",
                  cursor: "pointer",
                }}
                onClick={(e) => handleJobClick(job, e)}
              >
                <Avatar
                  sx={{
                    width: 32,
                    height: 32,
                    bgcolor: STATUS_COLORS[job.status] || "#AAA",
                  }}
                  src={`/svgicons/${job.task_name}.svg`}
                  imgProps={{
                    onError: (e: any) => {
                      e.target.src = `/qticons/${job.task_name}.png`;
                    },
                  }}
                >
                  {/* Fallback to first letter if no icon */}
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

// Simple SMILES viewer using RDKit (if available) or fallback
function SmilesViewer({
  smiles,
  width,
  height,
}: {
  smiles: string;
  width: number;
  height: number;
}) {
  // Try to use RDKit if available via window
  // This requires RDKit to be loaded (which CCP4i2 does have)
  const [svg, setSvg] = useState<string | null>(null);
  const [error, setError] = useState(false);

  // Attempt to render with RDKit
  useMemo(() => {
    if (typeof window !== "undefined" && (window as any).RDKit) {
      try {
        const RDKit = (window as any).RDKit;
        const mol = RDKit.get_mol(smiles);
        if (mol) {
          const svgString = mol.get_svg(width, height);
          setSvg(svgString);
          mol.delete();
        } else {
          setError(true);
        }
      } catch (e) {
        setError(true);
      }
    } else {
      // RDKit not loaded yet, show skeleton
      setError(true);
    }
  }, [smiles, width, height]);

  if (error) {
    // Fallback: show SMILES as tooltip on a chip
    return (
      <Tooltip title={smiles}>
        <Chip
          label={`NCL-${smiles.substring(0, 8)}...`}
          size="small"
          color="primary"
          variant="outlined"
        />
      </Tooltip>
    );
  }

  if (!svg) {
    return <Skeleton variant="rectangular" width={width} height={height} />;
  }

  return (
    <Box
      sx={{ width, height }}
      dangerouslySetInnerHTML={{ __html: svg }}
    />
  );
}
