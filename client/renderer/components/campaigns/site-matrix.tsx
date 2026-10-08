"use client";

/**
 * The pieces of the campaign overview's site matrix: a column header per
 * site, one box per dataset per site, and the legend that explains them.
 * The rules (which box, which colour, what the tooltip says, where a click
 * goes) live in lib/site-matrix and are tested there.
 */

import { Box, Stack, Theme, Tooltip, Typography } from "@mui/material";
import {
  CampaignSite,
  MemberProjectWithSummary,
  SiteVerdict,
} from "../../types/campaigns";
import {
  SiteCellKind,
  VERDICT_LEGEND,
  VerdictColourKey,
  cellFor,
  siteCellKind,
  siteCellTooltip,
  siteCellUrl,
  verdictColourKey,
} from "../../lib/site-matrix";

/** Width of a site column, in px. Narrow, because there may be 30 or more. */
export const SITE_COLUMN_WIDTH = 30;
const BOX_SIZE = 16;

/**
 * A palette colour that reads in both modes: the theme's own status colours
 * (which MUI already adjusts for dark mode), and a mid grey for "not
 * evaluated" that is visible on a light or a dark background alike.
 */
function paletteColour(theme: Theme, key: VerdictColourKey): string {
  return key === "grey" ? theme.palette.grey[500] : theme.palette[key].main;
}

function VerdictSwatch({
  kind,
  verdict,
}: {
  kind: "filled" | "outlined";
  verdict: SiteVerdict | null;
}) {
  const key = verdictColourKey(verdict);
  return (
    <Box
      data-testid="site-box"
      data-kind={kind}
      data-colour={key}
      sx={(theme) => {
        const colour = paletteColour(theme, key);
        return {
          width: BOX_SIZE,
          height: BOX_SIZE,
          boxSizing: "border-box",
          borderRadius: "3px",
          border: `2px solid ${colour}`,
          bgcolor: kind === "filled" ? colour : "transparent",
        };
      }}
    />
  );
}

export function SiteMatrixLegend() {
  return (
    <Stack
      direction="row"
      spacing={2}
      useFlexGap
      flexWrap="wrap"
      alignItems="center"
      sx={{ mb: 1 }}
    >
      <Stack direction="row" spacing={0.75} alignItems="center">
        <VerdictSwatch kind="filled" verdict={null} />
        <Typography variant="caption">PanDDA event</Typography>
      </Stack>
      <Stack direction="row" spacing={0.75} alignItems="center">
        <VerdictSwatch kind="outlined" verdict={null} />
        <Typography variant="caption">Verdict, no event</Typography>
      </Stack>
      <Typography variant="caption" color="text.secondary">
        Colour is the verdict:
      </Typography>
      {VERDICT_LEGEND.map(({ verdict, label }) => (
        <Stack
          key={label}
          direction="row"
          spacing={0.75}
          alignItems="center"
        >
          <VerdictSwatch kind="filled" verdict={verdict} />
          <Typography variant="caption">{label}</Typography>
        </Stack>
      ))}
    </Stack>
  );
}

/**
 * A site's column header: the name written upwards, so a column stays as
 * narrow as its boxes, cut off when long with the whole name on hover.
 */
export function SiteColumnHeaderContent({ site }: { site: CampaignSite }) {
  return (
    <Tooltip title={site.name}>
      <Box
        sx={{
          writingMode: "vertical-rl",
          transform: "rotate(180deg)",
          maxHeight: 96,
          overflow: "hidden",
          textOverflow: "ellipsis",
          whiteSpace: "nowrap",
          mx: "auto",
          fontSize: "0.75rem",
          lineHeight: 1.2,
        }}
      >
        {site.name}
      </Box>
    </Tooltip>
  );
}

interface SiteMatrixCellContentProps {
  project: MemberProjectWithSummary;
  site: CampaignSite;
  campaignId?: number;
}

/** What goes in one dataset's cell at one site. */
export function SiteMatrixCellContent({
  project,
  site,
  campaignId,
}: SiteMatrixCellContentProps) {
  const cell = cellFor(project, site);
  const kind: SiteCellKind = siteCellKind(cell, project.frame_mismatch);
  const url = siteCellUrl(campaignId, project.current_model_job, site);
  const lines = siteCellTooltip(
    site.name,
    cell,
    project.frame_mismatch,
    Boolean(project.current_model_job)
  );
  const title = (
    <Box>
      {lines.map((line, i) => (
        <Typography
          key={i}
          variant={i === 0 ? "body2" : "caption"}
          fontWeight={i === 0 ? "bold" : undefined}
          display="block"
        >
          {line}
        </Typography>
      ))}
    </Box>
  );

  if (kind === "mismatch") {
    return (
      <Tooltip title={title}>
        <Typography
          variant="body2"
          color="text.disabled"
          aria-label={lines.join(". ")}
          sx={{ textAlign: "center", cursor: "default" }}
        >
          {"–"}
        </Typography>
      </Tooltip>
    );
  }

  // An empty cell still carries the tooltip, so hovering anywhere in a column
  // says which site it is.
  const clickable = kind !== "none" && Boolean(url);
  return (
    <Tooltip title={title}>
      <Box
        role={clickable ? "link" : undefined}
        aria-label={lines.join(". ")}
        onClick={
          clickable
            ? (event) => {
                event.stopPropagation();
                window.open(url as string, "_blank");
              }
            : undefined
        }
        sx={{
          display: "flex",
          justifyContent: "center",
          alignItems: "center",
          minHeight: BOX_SIZE,
          cursor: clickable ? "pointer" : "default",
          "&:hover": clickable ? { opacity: 0.75 } : undefined,
        }}
      >
        {kind !== "none" && (
          <VerdictSwatch kind={kind} verdict={cell?.verdict ?? null} />
        )}
      </Box>
    </Tooltip>
  );
}
