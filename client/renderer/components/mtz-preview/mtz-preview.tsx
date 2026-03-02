"use client";
import { useMemo, useRef, useState } from "react";
import {
  Box,
  Chip,
  Slider,
  Stack,
  ToggleButton,
  ToggleButtonGroup,
  Typography,
  CircularProgress,
} from "@mui/material";
import { Panel, PanelGroup, PanelResizeHandle } from "react-resizable-panels";
import { parseMtzHeader, parseMtzReflections } from "../../lib/mtz-parser";
import {
  reciprocalLatticeVectors,
  computeSectionBasis,
  indexReflectionsBySection,
  getFixedAxisRange,
  getSectionAxes,
  expandReflectionsBySymmetry,
  type SectionPlane,
} from "../../lib/reciprocal-lattice";
import { spaceGroups } from "../../spacegroups";
import { MtzHeaderSummary } from "./mtz-header-summary";
import { ReciprocalLatticeCanvas } from "./reciprocal-lattice-canvas";

interface MtzPreviewProps {
  data: ArrayBuffer;
}

/** Format a number for display in the color range labels */
function formatValue(v: number): string {
  if (v === 0) return "0";
  const abs = Math.abs(v);
  if (abs >= 1000) return v.toExponential(1);
  if (abs >= 1) return v.toFixed(1);
  if (abs >= 0.01) return v.toFixed(3);
  return v.toExponential(1);
}

export const MtzPreview: React.FC<MtzPreviewProps> = ({ data }) => {
  const header = useMemo(() => parseMtzHeader(data), [data]);

  // Color column: null means auto-detect, "" means none/off
  const [colorColumn, setColorColumn] = useState<string | null>(null);
  const [sectionPlane, setSectionPlane] = useState<SectionPlane>("hk");
  const [sectionIndex, setSectionIndex] = useState(0);

  // Determine the effective intensity column to pass to the parser
  const effectiveColorColumn = colorColumn === "" ? undefined : (colorColumn ?? undefined);

  const reflectionData = useMemo(
    () => parseMtzReflections(data, header, effectiveColorColumn),
    [data, header, effectiveColorColumn]
  );

  // Color range: user's manual override (null = derive from data automatically)
  const [userColorRange, setUserColorRange] = useState<[number, number] | null>(null);

  // Reset user override when the intensity data source changes
  const prevIntensityRange = useRef(reflectionData.intensityRange);
  if (reflectionData.intensityRange !== prevIntensityRange.current) {
    prevIntensityRange.current = reflectionData.intensityRange;
    setUserColorRange(null);
  }

  // Effective range: user override, or auto [0, max] from data
  const colorRange: [number, number] = userColorRange ??
    (reflectionData.intensityRange ? [0, reflectionData.intensityRange[1]] : [0, 1]);

  // Whether coloring is active
  const coloringActive = colorColumn !== "" && reflectionData.intensityLabel !== undefined;

  // Compute reciprocal lattice vectors from cell parameters
  const recipVecs = useMemo(() => {
    const cell = header.cell;
    if (!cell) return null;
    return reciprocalLatticeVectors(
      cell.a, cell.b, cell.c,
      cell.alpha, cell.beta, cell.gamma
    );
  }, [header]);

  const sectionBasis = useMemo(() => {
    if (!recipVecs) return null;
    return computeSectionBasis(recipVecs, sectionPlane);
  }, [recipVecs, sectionPlane]);

  // Look up symmetry operators from space group number
  const symops = useMemo(() => {
    if (!header.spaceGroupNumber) return null;
    const sg = spaceGroups.find(sg => sg.numbers.includes(header.spaceGroupNumber!));
    return sg?.transformations ?? null;
  }, [header.spaceGroupNumber]);

  // Expand reflections by symmetry (ASU → full reciprocal sphere)
  // Skip for unmerged data — observations are already in the full sphere
  const expandedReflections = useMemo(() => {
    if (!header.isMerged) return reflectionData.reflections;
    if (!symops || symops.length <= 1) return reflectionData.reflections;
    return expandReflectionsBySymmetry(reflectionData.reflections, symops);
  }, [reflectionData.reflections, symops, header.isMerged]);

  const sectionMap = useMemo(() => {
    if (!sectionBasis) return new Map();
    return indexReflectionsBySection(
      expandedReflections,
      sectionPlane,
      sectionBasis
    );
  }, [expandedReflections, sectionPlane, sectionBasis]);

  const sectionData = useMemo(
    () => sectionMap.get(sectionIndex) ?? [],
    [sectionMap, sectionIndex]
  );

  // Compute hkl ranges from expanded reflections for slider range
  const expandedRanges = useMemo(() => {
    let hMin = Infinity, hMax = -Infinity;
    let kMin = Infinity, kMax = -Infinity;
    let lMin = Infinity, lMax = -Infinity;
    for (const r of expandedReflections) {
      if (r.h < hMin) hMin = r.h;
      if (r.h > hMax) hMax = r.h;
      if (r.k < kMin) kMin = r.k;
      if (r.k > kMax) kMax = r.k;
      if (r.l < lMin) lMin = r.l;
      if (r.l > lMax) lMax = r.l;
    }
    return {
      hRange: [hMin, hMax] as [number, number],
      kRange: [kMin, kMax] as [number, number],
      lRange: [lMin, lMax] as [number, number],
    };
  }, [expandedReflections]);

  const fixedRange = useMemo(
    () =>
      getFixedAxisRange(
        expandedRanges.hRange,
        expandedRanges.kRange,
        expandedRanges.lRange,
        sectionPlane
      ),
    [expandedRanges, sectionPlane]
  );

  const [uLabel, vLabel, fixedLabel] = getSectionAxes(sectionPlane);

  const handlePlaneChange = (_: unknown, value: SectionPlane | null) => {
    if (!value) return;
    setSectionPlane(value);
    setSectionIndex(0);
  };

  // Slider range: use the full intensity range for slider min/max
  const sliderMin = reflectionData.intensityRange?.[0] ?? 0;
  const sliderMax = reflectionData.intensityRange?.[1] ?? 1;

  // Log-scale slider: map slider position [0,100] to value via log interpolation
  const logMin = sliderMin > 0 ? Math.log10(sliderMin) : 0;
  const logMax = sliderMax > 0 ? Math.log10(sliderMax) : 1;

  const valueToSlider = (v: number): number => {
    if (v <= 0) return 0;
    const logV = Math.log10(v);
    return ((logV - logMin) / (logMax - logMin)) * 100;
  };

  const sliderToValue = (s: number): number => {
    if (s <= 0) return 0;
    const logV = logMin + (s / 100) * (logMax - logMin);
    return Math.pow(10, logV);
  };

  const handleColorRangeChange = (_: unknown, newValue: number | number[]) => {
    const [lo, hi] = newValue as number[];
    setUserColorRange([sliderToValue(lo), sliderToValue(hi)]);
  };

  if (!header) {
    return (
      <Box sx={{ display: "flex", justifyContent: "center", alignItems: "center", height: "100%" }}>
        <CircularProgress />
      </Box>
    );
  }

  return (
    <Box sx={{ height: "calc(100vh - 20rem)", minHeight: 400 }}>
      <PanelGroup direction="horizontal">
        <Panel defaultSize={35} minSize={20}>
          <Box sx={{ height: "100%", overflowY: "auto", p: 2 }}>
            <MtzHeaderSummary
              header={header}
              onColorColumnChange={setColorColumn}
              selectedColorColumn={colorColumn}
              autoDetectedColumn={reflectionData.intensityLabel ?? null}
            />
          </Box>
        </Panel>
        <PanelResizeHandle
          style={{
            width: 10,
            backgroundColor: "transparent",
            display: "flex",
            alignItems: "center",
            justifyContent: "center",
            cursor: "col-resize",
          }}
        >
          <div
            style={{
              width: 4,
              height: "50%",
              backgroundColor: "gray",
              borderRadius: 2,
            }}
          />
        </PanelResizeHandle>
        <Panel defaultSize={65} minSize={30}>
          <Stack sx={{ height: "100%" }} spacing={0}>
            {/* Controls toolbar */}
            <Stack
              direction="row"
              spacing={2}
              alignItems="center"
              sx={{ px: 2, py: 1, flexShrink: 0 }}
            >
              <ToggleButtonGroup
                value={sectionPlane}
                exclusive
                onChange={handlePlaneChange}
                size="small"
              >
                <ToggleButton value="hk">hk</ToggleButton>
                <ToggleButton value="hl">hl</ToggleButton>
                <ToggleButton value="kl">kl</ToggleButton>
              </ToggleButtonGroup>
              <Typography
                variant="body2"
                sx={{ minWidth: 40, fontFamily: "monospace" }}
              >
                {fixedLabel}={sectionIndex}
              </Typography>
              <Slider
                value={sectionIndex}
                min={fixedRange[0]}
                max={fixedRange[1]}
                step={1}
                onChange={(_, v) => setSectionIndex(v as number)}
                sx={{ maxWidth: 200, minWidth: 100 }}
                valueLabelDisplay="auto"
                valueLabelFormat={(v) => `${fixedLabel}=${v}`}
              />
              <Chip
                label={`${sectionData.length} refl.`}
                size="small"
                variant="outlined"
                title={symops && symops.length > 1
                  ? `${reflectionData.reflections.length} unique → ${expandedReflections.length} expanded`
                  : undefined}
              />
            </Stack>
            {/* Color range slider */}
            {coloringActive && reflectionData.intensityRange && (
              <Stack
                direction="row"
                spacing={1}
                alignItems="center"
                sx={{ px: 2, pb: 0.5, flexShrink: 0 }}
              >
                <Typography variant="caption" color="text.secondary" sx={{ minWidth: 70 }}>
                  {reflectionData.intensityLabel}
                </Typography>
                <Typography variant="caption" sx={{ fontFamily: "monospace", minWidth: 50 }}>
                  {formatValue(colorRange[0])}
                </Typography>
                <Slider
                  value={[valueToSlider(colorRange[0]), valueToSlider(colorRange[1])]}
                  onChange={handleColorRangeChange}
                  min={0}
                  max={100}
                  step={0.5}
                  sx={{ minWidth: 120 }}
                  size="small"
                />
                <Typography variant="caption" sx={{ fontFamily: "monospace", minWidth: 50 }}>
                  {formatValue(colorRange[1])}
                </Typography>
              </Stack>
            )}
            {/* Canvas */}
            <Box sx={{ flex: 1, minHeight: 0 }}>
              <ReciprocalLatticeCanvas
                sectionData={sectionData}
                uLabel={`${uLabel}*`}
                vLabel={`${vLabel}*`}
                showIntensity={coloringActive}
                colorRange={coloringActive ? colorRange : undefined}
                sectionBasis={sectionBasis}
                resolution={header.resolution}
                recipVecs={recipVecs}
                intensityLabel={reflectionData.intensityLabel}
              />
            </Box>
          </Stack>
        </Panel>
      </PanelGroup>
    </Box>
  );
};
