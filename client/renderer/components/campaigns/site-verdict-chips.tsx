"use client";

/**
 * A dataset's Sites cell in the campaign overview.
 *
 * Shows where something was found, and — in one quiet line underneath — how
 * far the evaluation of this dataset has got. What it deliberately does not
 * chip is the empties: in a campaign with 30-40 sites most are empty for most
 * datasets, and chipping them would bury the one or two results somebody is
 * scanning the column for.
 *
 * That makes three different datasets look alike from the chips alone: one
 * nobody has opened, one part-way through with nothing found so far, and one
 * examined throughout and found empty. The note is what tells them apart, and
 * it says least on the rows the chips already explain.
 */

import { Box, Chip, Stack, Tooltip, Typography } from "@mui/material";
import {
  MemberProjectWithSummary,
  SiteEvaluationSummary,
} from "../../types/campaigns";
import {
  EvaluationTone,
  evaluationNote,
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

/**
 * The note is a footnote to the chips, never a competitor to them: an
 * unevaluated dataset is the quietest thing in the column, because a column of
 * them is the normal state of a campaign early on.
 */
const NOTE_SX: Record<EvaluationTone, object> = {
  progress: { color: "text.secondary" },
  clear: { color: "text.secondary" },
  untouched: { color: "text.disabled", fontStyle: "italic" },
};

export function SiteVerdictChips({
  project,
  campaignId,
  maxChips = 3,
}: SiteVerdictChipsProps) {
  const { shown, overflow } = verdictChips(project.site_evaluations, maxChips);
  const note = evaluationNote(
    project.sites_evaluated,
    project.sites_total,
    (project.site_evaluations?.length ?? 0) > 0
  );
  const jobId = preferredJobId(project.jobs);

  if (shown.length === 0 && !note) return null;

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
      {note && (
        <Tooltip title={note.detail}>
          <Box>
            <Typography variant="caption" sx={NOTE_SX[note.tone]}>
              {note.text}
            </Typography>
          </Box>
        </Tooltip>
      )}
    </Stack>
  );
}
