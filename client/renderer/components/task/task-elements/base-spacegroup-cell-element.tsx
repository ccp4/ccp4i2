import { Box, Chip, Stack, Typography } from "@mui/material";
import GridOnIcon from "@mui/icons-material/GridOn";

interface CCell {
  a: number;
  b: number;
  c: number;
  alpha: number;
  beta: number;
  gamma: number;
}

interface CResolutionRange {
  low: number;
  high: number;
}

/**
 * Digest data from MTZ/reflection files.
 * Supports both formats:
 * - Legacy: resolutionRange.low / resolutionRange.high
 * - Current: lowRes / highRes (from digest API)
 */
interface CObsData {
  cell?: CCell;
  spaceGroup?: string;
  resolutionRange?: CResolutionRange;
  // Alternative resolution format from digest API
  lowRes?: number;
  highRes?: number;
}

interface BaseSpacegroupCellElementProps {
  data?: CObsData;
}

/** Format cell parameter value with appropriate precision */
const formatCellParam = (value: number, isAngle: boolean): string => {
  if (isAngle) {
    return value.toFixed(1) + "°";
  }
  return value.toFixed(2) + " Å";
};

/** Format resolution value */
const formatResolution = (value: number | undefined): string => {
  if (typeof value !== "number") return "—";
  return value.toFixed(2) + " Å";
};

export const BaseSpacegroupCellElement: React.FC<
  BaseSpacegroupCellElementProps
> = ({ data }) => {
  if (!data?.cell) {
    return null;
  }

  const { cell, spaceGroup } = data;
  const lowRes = data.resolutionRange?.low ?? data.lowRes;
  const highRes = data.resolutionRange?.high ?? data.highRes;

  return (
    <Box
      sx={{
        p: 1.5,
        borderRadius: 1,
        bgcolor: "action.hover",
        border: 1,
        borderColor: "divider",
      }}
    >
      <Stack spacing={1.5}>
        {/* Spacegroup */}
        <Stack direction="row" alignItems="center" spacing={1}>
          <GridOnIcon fontSize="small" color="action" />
          <Typography variant="body2" color="text.secondary">
            Spacegroup:
          </Typography>
          <Chip
            label={spaceGroup || "Unknown"}
            size="small"
            variant="outlined"
            color="primary"
          />
        </Stack>

        {/* Cell parameters */}
        <Box>
          <Typography
            variant="caption"
            color="text.secondary"
            sx={{ mb: 0.5, display: "block" }}
          >
            Unit cell
          </Typography>
          <Stack direction="row" spacing={1} flexWrap="wrap" useFlexGap>
            <Chip
              label={`a = ${formatCellParam(cell.a, false)}`}
              size="small"
              variant="outlined"
            />
            <Chip
              label={`b = ${formatCellParam(cell.b, false)}`}
              size="small"
              variant="outlined"
            />
            <Chip
              label={`c = ${formatCellParam(cell.c, false)}`}
              size="small"
              variant="outlined"
            />
            <Chip
              label={`α = ${formatCellParam(cell.alpha, true)}`}
              size="small"
              variant="outlined"
            />
            <Chip
              label={`β = ${formatCellParam(cell.beta, true)}`}
              size="small"
              variant="outlined"
            />
            <Chip
              label={`γ = ${formatCellParam(cell.gamma, true)}`}
              size="small"
              variant="outlined"
            />
          </Stack>
        </Box>

        {/* Resolution */}
        <Stack direction="row" alignItems="center" spacing={1}>
          <Typography variant="body2" color="text.secondary">
            Resolution:
          </Typography>
          <Typography variant="body2">
            {formatResolution(lowRes)} – {formatResolution(highRes)}
          </Typography>
        </Stack>
      </Stack>
    </Box>
  );
};
