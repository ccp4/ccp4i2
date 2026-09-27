import { useMemo } from "react";
import { Alert, AlertTitle, Box, Typography } from "@mui/material";

import { useApi } from "../../../api";

/**
 * What resolution this campaign will actually be processed at, said before
 * anyone pays for a run.
 *
 * PanDDA processes a dataset at the WORST resolution among its comparators,
 * and takes comparators regardless of resolution until it has
 * max_shell_datasets of them (default 60). So a campaign with fewer datasets
 * than that is processed at the resolution of its worst crystal however good
 * the rest are -- and nothing said so until the run was over and the report
 * printed a median processing resolution as though it were a fact about the
 * data.
 *
 * DDU's 50-dataset CDK4 campaign was processed at 6.71 A, the resolution of
 * one crystal, while 48 of its 50 datasets were better than 4 A. No fragment
 * is findable there. Removing 14 low-resolution datasets and running the same
 * software over the same crystals gave 2.97 A, interpretable density, and
 * poses that sit in it.
 *
 * So this offers the two remedies with their arithmetic done, rather than a
 * warning the reader has to act on blind.
 */
interface ResolutionOutlook {
  datasets: { dtag: string; resolution: number | null }[];
  best: number | null;
  worst: number | null;
  processing: number | null;
  dragged_by_comparator_floor: boolean;
  max_shell_datasets: number;
}

/** The resolution a set of datasets would be processed at with a given shell. */
export function shellResolution(resolutions: number[], maxShellDatasets: number): number | null {
  if (resolutions.length === 0) return null;
  const sorted = [...resolutions].sort((a, b) => a - b);
  return sorted[Math.min(sorted.length, maxShellDatasets) - 1];
}

export const PanddaResolutionNotice: React.FC<{ jobId: number }> = ({ jobId }) => {
  const api = useApi();
  const { data } = api.get<{ data?: ResolutionOutlook } | ResolutionOutlook>(
    jobId ? `jobs/${jobId}/pandda_resolution/` : null
  );
  const outlook = (data as { data?: ResolutionOutlook })?.data ?? (data as ResolutionOutlook);

  const remedies = useMemo(() => {
    if (!outlook) return null;
    const resolutions = outlook.datasets
      .map((d) => d.resolution)
      .filter((r): r is number => typeof r === "number");
    if (resolutions.length === 0 || outlook.best == null) return null;

    // What a cut would give: the worst dataset kept sets the resolution.
    const cuts = [4.0, 3.5, 3.0, 2.8]
      .map((cut) => {
        const kept = resolutions.filter((r) => r <= cut);
        return { cut, kept: kept.length, resolution: kept.length ? Math.max(...kept) : null };
      })
      // Only cuts that lose something, keep enough to characterise a ground
      // state, and actually improve on what we have.
      .filter(
        (c) =>
          c.kept >= 25 &&
          c.kept < resolutions.length &&
          c.resolution != null &&
          outlook.processing != null &&
          c.resolution < outlook.processing - 0.1
      );

    // What a smaller shell would give, keeping every dataset.
    const shells = [40, 30, 25]
      .map((n) => ({ n, resolution: shellResolution(resolutions, n) }))
      .filter(
        (s) =>
          s.n < resolutions.length &&
          s.resolution != null &&
          outlook.processing != null &&
          s.resolution < outlook.processing - 0.1
      );

    return { cuts, shells, count: resolutions.length };
  }, [outlook]);

  if (!outlook || !outlook.dragged_by_comparator_floor || !remedies) return null;
  const { processing, best, worst, max_shell_datasets: maxShell } = outlook;
  if (processing == null || best == null || processing <= best + 0.5) return null;

  const worstDatasets = outlook.datasets
    .filter((d) => d.resolution != null && worst != null && d.resolution >= worst - 0.01)
    .map((d) => d.dtag);

  return (
    <Alert severity="warning" sx={{ mt: 1 }}>
      <AlertTitle>
        This will be processed at {processing.toFixed(2)} A, not {best.toFixed(2)} A
      </AlertTitle>
      <Typography variant="body2">
        PanDDA processes each dataset at the worst resolution among its comparators, and takes
        comparators regardless of resolution until it has {maxShell} of them. With{" "}
        {remedies.count} datasets that floor is never reached, so every dataset is a comparator and
        the worst crystal
        {worstDatasets.length ? ` (${worstDatasets.slice(0, 2).join(", ")})` : ""} sets the
        resolution for the whole run.
      </Typography>
      {remedies.cuts.length > 0 && (
        <Box sx={{ mt: 1 }}>
          <Typography variant="body2" sx={{ fontWeight: 600 }}>
            Remove the low-resolution datasets:
          </Typography>
          {remedies.cuts.map((c) => (
            <Typography variant="body2" key={c.cut} sx={{ ml: 1 }}>
              keep the {c.kept} better than {c.cut.toFixed(1)} A &rarr; processed at{" "}
              {c.resolution!.toFixed(2)} A
            </Typography>
          ))}
        </Box>
      )}
      {remedies.shells.length > 0 && (
        <Box sx={{ mt: 1 }}>
          <Typography variant="body2" sx={{ fontWeight: 600 }}>
            Or keep them all and lower MAX_SHELL_DATASETS (on the Run tab):
          </Typography>
          {remedies.shells.map((s) => (
            <Typography variant="body2" key={s.n} sx={{ ml: 1 }}>
              {s.n} &rarr; processed at {s.resolution!.toFixed(2)} A
            </Typography>
          ))}
        </Box>
      )}
      <Typography variant="caption" color="text.secondary" sx={{ display: "block", mt: 1 }}>
        Comparators must stay above MIN_CHARACTERISATION_DATASETS (25) either way.
      </Typography>
    </Alert>
  );
};
