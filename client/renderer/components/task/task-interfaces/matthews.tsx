/*
 * Copyright (C) 2026 Newcastle University
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
import {
  Alert,
  Box,
  Grid2,
  LinearProgress,
  Paper,
  Stack,
  Typography,
} from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "./task-container";
import { CCP4i2TaskElement } from "../task-elements/task-element";
import { CCP4i2ContainerElement } from "../task-elements/ccontainer";
import { useJob } from "../../../utils";
import { useApi } from "../../../api";
import { BaseSpacegroupCellElement } from "../task-elements/base-spacegroup-cell-element";

/**
 * Renders a Matthews probability table showing nmol/ASU, solvent %, and
 * probability bars.  Reused from ProvideAsuContents.
 */
const MatthewsResultsBox: React.FC<{ results: any[] }> = ({ results }) => (
  <Box
    sx={{
      p: 1.5,
      borderRadius: 1,
      bgcolor: "action.hover",
      border: 1,
      borderColor: "divider",
    }}
  >
    <Typography
      variant="caption"
      color="text.secondary"
      sx={{ mb: 1, display: "block" }}
    >
      Matthews Analysis
    </Typography>
    <Stack spacing={1}>
      {results.map(
        (result: {
          nmol_in_asu: number;
          percent_solvent: number;
          prob_matth: number;
        }) => {
          const probability = result.prob_matth;
          const isLikely = probability > 0.5;
          return (
            <Box
              key={result.nmol_in_asu}
              sx={{
                p: 1,
                borderRadius: 0.5,
                bgcolor: isLikely ? "success.main" : "background.paper",
                color: isLikely ? "success.contrastText" : "text.primary",
                border: 1,
                borderColor: isLikely ? "success.main" : "divider",
              }}
            >
              <Stack
                direction="row"
                justifyContent="space-between"
                alignItems="center"
              >
                <Typography variant="body2" fontWeight="medium">
                  {result.nmol_in_asu} mol/ASU
                </Typography>
                <Typography variant="body2" fontWeight="bold">
                  {(probability * 100).toFixed(0)}%
                </Typography>
              </Stack>
              <Stack direction="row" spacing={2} sx={{ mt: 0.5 }}>
                <Typography
                  variant="caption"
                  sx={{ opacity: isLikely ? 0.9 : 0.7 }}
                >
                  {result.percent_solvent.toFixed(1)}% solvent
                </Typography>
                <LinearProgress
                  variant="determinate"
                  value={probability * 100}
                  sx={{
                    flex: 1,
                    alignSelf: "center",
                    height: 4,
                    borderRadius: 2,
                    bgcolor: isLikely
                      ? "success.light"
                      : "action.disabledBackground",
                    "& .MuiLinearProgress-bar": {
                      bgcolor: isLikely
                        ? "success.contrastText"
                        : "primary.main",
                    },
                  }}
                />
              </Stack>
            </Box>
          );
        }
      )}
    </Stack>
  </Box>
);

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { job } = props;
  const api = useApi();
  const { useTaskItem, useFileDigest } = useJob(job.id);

  // Track input values
  const { value: modeValue } = useTaskItem("MODE");
  const { value: molwtValue } = useTaskItem("MOLWT");
  const { value: nresValue } = useTaskItem("NRES");
  const { value: hklinValue } = useTaskItem("HKLIN");

  const mode = modeValue ?? "asu_components";

  // File digest for HKLIN — gives us cell parameters and space group
  const hasHKLINFile = Boolean(hklinValue?.dbFileId);
  const hklinDigestPath = hasHKLINFile ? "matthews.inputData.HKLIN" : "";
  const { data: HKLINDigest, mutate: mutateHKLINDigest } =
    useFileDigest(hklinDigestPath);

  // For ASU mode: get molecular weight from the ASU file
  const { item: asuItem } = useTaskItem("ASUIN");
  const hasAsuFile = Boolean(asuItem?._value?.dbFileId);
  // Fetch digest to trigger re-fetch when ASU file changes
  const asuDigestPath = hasAsuFile ? "matthews.inputData.ASUIN" : "";
  const { data: asuDigest } = useFileDigest(asuDigestPath);

  const { data: molWeight, mutate: mutateMolWeight } = api.objectMethod<any>(
    job.id,
    "matthews.inputData.ASUIN.fileContent.seqList",
    "molecularWeight",
    {},
    [asuDigest],
    hasAsuFile && mode === "asu_components"
  );

  // Determine the molWt/nRes to pass to matthewsCoeff
  const effectiveMolWt =
    mode === "asu_components"
      ? molWeight?.data?.result
      : mode === "molwt"
        ? molwtValue
        : undefined;
  const effectiveNres = mode === "nres" ? nresValue : undefined;

  // Build kwargs for matthewsCoeff
  const matthewsKwargs: Record<string, any> = {};
  if (effectiveNres) matthewsKwargs.nRes = effectiveNres;
  else if (effectiveMolWt) matthewsKwargs.molWt = effectiveMolWt;

  const canFetchMatthews =
    hasHKLINFile &&
    HKLINDigest &&
    (effectiveMolWt || effectiveNres);

  const { data: matthewsAnalysis, mutate: mutateMatthews } =
    api.objectMethod<any>(
      job.id,
      "matthews.inputData.HKLIN.fileContent",
      "matthewsCoeff",
      matthewsKwargs,
      [effectiveMolWt, effectiveNres, HKLINDigest],
      !!canFetchMatthews
    );

  const matthewsResults = matthewsAnalysis?.data?.result?.results;

  return (
    <Paper sx={{ display: "flex", flexDirection: "column", gap: 1, p: 1 }}>
      <CCP4i2ContainerElement
        {...props}
        itemName=""
        qualifiers={{ guiLabel: "Cell parameters from reflection data" }}
        containerHint="FolderLevel"
        initiallyOpen={true}
      >
        <CCP4i2TaskElement
          itemName="HKLIN"
          {...props}
          onChange={() => {
            mutateHKLINDigest();
            mutateMatthews();
          }}
        />
        {HKLINDigest && (
          <Stack spacing={1} sx={{ mt: 1 }}>
            <BaseSpacegroupCellElement data={HKLINDigest} />
          </Stack>
        )}
      </CCP4i2ContainerElement>

      <CCP4i2ContainerElement
        {...props}
        itemName=""
        qualifiers={{ guiLabel: "Calculate molecular weight from..." }}
        containerHint="FolderLevel"
        initiallyOpen={true}
      >
        <CCP4i2TaskElement
          itemName="MODE"
          {...props}
          qualifiers={{ guiLabel: "Mode" }}
          onChange={() => mutateMatthews()}
        />

        {mode === "asu_components" && (
          <CCP4i2TaskElement
            itemName="ASUIN"
            {...props}
            qualifiers={{ guiLabel: "Contents of biological unit" }}
            onChange={async () => {
              await mutateMolWeight();
              mutateMatthews();
            }}
          />
        )}
        {mode === "nres" && (
          <CCP4i2TaskElement
            itemName="NRES"
            {...props}
            qualifiers={{ guiLabel: "Number of residues" }}
            onChange={() => mutateMatthews()}
          />
        )}
        {mode === "molwt" && (
          <CCP4i2TaskElement
            itemName="MOLWT"
            {...props}
            qualifiers={{ guiLabel: "Molecular weight (Da)" }}
            onChange={() => mutateMatthews()}
          />
        )}

      </CCP4i2ContainerElement>

      {/* Molecular weight feedback + Matthews results side by side */}
      <Grid2 container spacing={2}>
        <Grid2 size={{ xs: 12, sm: 4 }}>
          {mode === "asu_components" && molWeight?.data?.result && (
            <Typography variant="body2" color="text.secondary">
              Molecular weight: {molWeight.data.result.toFixed(0)} Da
            </Typography>
          )}
          {mode === "nres" && effectiveNres && (
            <Typography variant="body2" color="text.secondary">
              Est. MW: {(effectiveNres * 112.5).toFixed(0)} Da
            </Typography>
          )}
          {mode === "molwt" && effectiveMolWt && (
            <Typography variant="body2" color="text.secondary">
              Molecular weight: {Number(effectiveMolWt).toFixed(0)} Da
            </Typography>
          )}
          {!canFetchMatthews && (
            <Alert severity="info" sx={{ py: 0, mt: 1 }}>
              Provide MTZ + molecular weight to calculate Matthews coefficient
            </Alert>
          )}
        </Grid2>
        <Grid2 size={{ xs: 12, sm: 8 }}>
          {matthewsResults?.length > 0 ? (
            <MatthewsResultsBox results={matthewsResults} />
          ) : canFetchMatthews ? (
            <Typography variant="body2" color="text.secondary">
              Calculating...
            </Typography>
          ) : null}
        </Grid2>
      </Grid2>
    </Paper>
  );
};

export default TaskInterface;
