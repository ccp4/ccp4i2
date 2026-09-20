"use client";

/**
 * A dataset's Sites cell in the campaign overview.
 *
 * Shows where something was found, and how much of the dataset is still
 * unlooked at. What it deliberately does not show is the empties: in a
 * campaign with 30-40 sites most are empty for most datasets, and chipping
 * them would bury the one or two results somebody is scanning the column for.
 * A dataset nobody has looked at and one examined and found empty both render
 * as nothing here; the count below tells them apart, and the per-dataset view
 * in Moorhen has the detail.
 */

import { Box, Chip, Stack, Tooltip, Typography } from "@mui/material";
import {
  MemberProjectWithSummary,
  SiteEvaluationSummary,
} from "../../types/campaigns";
import {
  evaluationProgress,
  preferredJobId,
  siteViewUrl,
  verdictChips,
} from "../../lib/site-verdicts";

interface SiteVerdictChipsProps {
  project: MemberProjectWithSummary;
  /** Absent outside a campaign, in which case the chips are not links. */
  campaignId?: number;
  /** How many chips before the rest collapse into a "+N". */
  maxChips?: number;
}

export function SiteVerdictChips({
  project,
  campaignId,
  maxChips = 3,
}: SiteVerdictChipsProps) {
  const { shown, overflow } = verdictChips(project.site_evaluations, maxChips);
  const progress = evaluationProgress(
    project.sites_evaluated,
    project.sites_total
  );
  const jobId = preferredJobId(project.jobs);

  if (shown.length === 0 && !progress) return null;

  const open = (evaluation: SiteEvaluationSummary) => {
    const url = siteViewUrl(campaignId, jobId, evaluation.site_id);
    if (url) window.open(url, "_blank");
  };

  return (
    <Stack spacing={0.25} alignItems="flex-start">
      <Stack direction="row" spacing={0.5} flexWrap="wrap" useFlexGap>
        {shown.map((evaluation) => {
          const clickable = Boolean(siteViewUrl(campaignId, jobId, 0));
          const isHit = evaluation.verdict === "hit";
          return (
            <Tooltip
              key={evaluation.site_id}
              title={
                clickable
                  ? `${evaluation.site_name} — ${
                      isHit ? "ligand present" : "unclear"
                    }. Click to view.`
                  : `${evaluation.site_name} — ${
                      isHit ? "ligand present" : "unclear"
                    }`
              }
            >
              <Chip
                label={evaluation.site_name}
                size="small"
                // A hit is the result; an unclear is deliberately quieter, so
                // a column of maybes does not read like a column of findings.
                color={isHit ? "success" : "warning"}
                variant={isHit ? "filled" : "outlined"}
                onClick={
                  clickable
                    ? (event) => {
                        event.stopPropagation();
                        open(evaluation);
                      }
                    : undefined
                }
                sx={{
                  maxWidth: 140,
                  cursor: clickable ? "pointer" : "default",
                }}
              />
            </Tooltip>
          );
        })}
        {overflow > 0 && (
          <Tooltip title={`${overflow} more`}>
            <Chip
              label={`+${overflow}`}
              size="small"
              variant="outlined"
              sx={{ cursor: "default" }}
            />
          </Tooltip>
        )}
      </Stack>
      {progress && (
        <Tooltip title="Sites evaluated in this dataset, of the campaign's sites">
          <Box>
            <Typography variant="caption" color="text.secondary">
              {progress}
            </Typography>
          </Box>
        </Tooltip>
      )}
    </Stack>
  );
}
