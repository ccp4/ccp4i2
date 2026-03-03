import { useMemo } from "react";
import $ from "jquery";
//@ts-ignore
import CytoscapeComponent from "react-cytoscapejs";
import Cytoscape from "cytoscape";
import COSEBilkent from "cytoscape-cose-bilkent";
import { Box, Typography } from "@mui/material";
import { useTheme } from "@mui/material/styles";
import { CCP4i2ReportElementProps } from "./CCP4i2ReportElement";

Cytoscape.use(COSEBilkent);

/** Color palette for different pipeline step types. */
const STEP_COLORS: Record<string, { bg: string; bgDark: string }> = {
  // Data preparation
  WALK_DIRECTORY: { bg: "#e8f5e9", bgDark: "#1b5e20" },
  MTZ_FILE_READER: { bg: "#e8f5e9", bgDark: "#1b5e20" },
  DATA_PREPARATION: { bg: "#e8f5e9", bgDark: "#1b5e20" },
  BIOLOGICAL_UNIT_BUILDER: { bg: "#e8f5e9", bgDark: "#1b5e20" },
  // Analysis
  ANISOTROPY_CORRECTION: { bg: "#fff3e0", bgDark: "#e65100" },
  SPACE_GROUP_ANALYSIS: { bg: "#fff3e0", bgDark: "#e65100" },
  TNCS_COMMENSURATE_MODULATION_ANALYSIS: { bg: "#fff3e0", bgDark: "#e65100" },
  TNCS_CORRECTION_REFINEMENT: { bg: "#fff3e0", bgDark: "#e65100" },
  TWIN_TEST: { bg: "#fff3e0", bgDark: "#e65100" },
  CELL_CONTENT_ANALYSIS: { bg: "#fff3e0", bgDark: "#e65100" },
  CELL_CONTENT_SCALING: { bg: "#fff3e0", bgDark: "#e65100" },
  // Search
  ENSEMBLE_SINGLE_MODEL: { bg: "#e3f2fd", bgDark: "#0d47a1" },
  FIND_COMPONENTS: { bg: "#e3f2fd", bgDark: "#0d47a1" },
  FAST_ROTATION_FUNCTION: { bg: "#e1f5fe", bgDark: "#01579b" },
  FAST_ROTATION_FUNCTION_RESCORE: { bg: "#e1f5fe", bgDark: "#01579b" },
  GYRE_REFINEMENT: { bg: "#e1f5fe", bgDark: "#01579b" },
  FAST_TRANSLATION_FUNCTION: { bg: "#e8eaf6", bgDark: "#1a237e" },
  FAST_TRANSLATION_FUNCTION_RESCORE: { bg: "#e8eaf6", bgDark: "#1a237e" },
  // Evaluation
  POSE_SCORING: { bg: "#f3e5f5", bgDark: "#4a148c" },
  PACKING: { bg: "#f3e5f5", bgDark: "#4a148c" },
  PERMUTATION_FOR_SEARCH: { bg: "#f3e5f5", bgDark: "#4a148c" },
  // Refinement & scoring
  RIGID_BODY_REFINEMENT: { bg: "#fce4ec", bgDark: "#880e4f" },
  MOVING: { bg: "#fce4ec", bgDark: "#880e4f" },
  RIGID_BODY_MAPS_FOR_INTENSITIES: { bg: "#fce4ec", bgDark: "#880e4f" },
  XRAY_COORDINATE_REFINEMENT: { bg: "#fce4ec", bgDark: "#880e4f" },
  TRANSLATION_FUNCTION_ZSCORE_EQUIVALENT: { bg: "#fff9c4", bgDark: "#f57f17" },
  "R-FACTOR_CALCULATION": { bg: "#fff9c4", bgDark: "#f57f17" },
  // Utility
  EFFECTIVE_AMPLITUDES: { bg: "#f5f5f5", bgDark: "#424242" },
  JOIN_NODES: { bg: "#f5f5f5", bgDark: "#424242" },
  SPAN_SPACE: { bg: "#f5f5f5", bgDark: "#424242" },
};

const DEFAULT_COLOR = { bg: "#f5f5f5", bgDark: "#616161" };

export const CCP4i2ReportDAG: React.FC<CCP4i2ReportElementProps> = ({
  item,
}) => {
  const theme = useTheme();
  const isDark = theme.palette.mode === "dark";

  const { elements, title } = useMemo(() => {
    const titleAttr = $(item).attr("title") || "Pipeline DAG";
    // The JSON is embedded as text content of the XML element
    const jsonText = $(item).text().trim();
    let parsed: any[] = [];
    try {
      parsed = JSON.parse(jsonText);
    } catch {
      // Fall back to empty
    }
    return { elements: parsed, title: titleAttr };
  }, [item]);

  const cytoscapeStyles = useMemo(
    () => [
      {
        selector: "node",
        style: {
          label: "data(label)",
          "text-wrap": "wrap" as const,
          "text-max-width": 160,
          "font-size": 9,
          shape: "round-rectangle" as const,
          width: 180,
          height: 36,
          "text-valign": "center" as const,
          "text-halign": "center" as const,
          "background-color": isDark ? "#424242" : "#f5f5f5",
          "border-width": 1,
          "border-color": isDark ? "#757575" : "#bdbdbd",
          color: isDark ? "#ffffff" : "#000000",
        },
      },
      // Generate per-step-type selectors for coloring
      ...Object.entries(STEP_COLORS).map(([tag, colors]) => ({
        selector: `node[tag = "${tag}"]`,
        style: {
          "background-color": isDark ? colors.bgDark : colors.bg,
          "border-color": isDark ? "#757575" : "#9e9e9e",
        },
      })),
      // Container/grouping nodes (Picard, Space_Groups, etc.)
      {
        selector: "node[?info]",
        style: {
          "font-weight": "bold" as const,
          "border-width": 2,
        },
      },
      {
        selector: "edge",
        style: {
          width: 1.5,
          "curve-style": "bezier" as const,
          "line-color": isDark ? "#616161" : "#bdbdbd",
          "target-arrow-shape": "triangle" as const,
          "target-arrow-color": isDark ? "#616161" : "#bdbdbd",
          "arrow-scale": 0.8,
        },
      },
    ],
    [isDark]
  );

  const layoutOptions = useMemo(
    () => ({
      name: "breadthfirst",
      animate: false,
      directed: true,
      spacingFactor: 0.7,
      fit: true,
      padding: 20,
    }),
    []
  );

  if (elements.length === 0) {
    return null;
  }

  return (
    <Box sx={{ my: 1 }}>
      <Typography variant="subtitle2" sx={{ mb: 1 }}>
        {title}
      </Typography>
      <Box
        sx={{
          border: 1,
          borderColor: "divider",
          borderRadius: 1,
          overflow: "hidden",
        }}
      >
        <CytoscapeComponent
          elements={elements}
          stylesheet={cytoscapeStyles}
          style={{
            width: "100%",
            height: "500px",
          }}
          layout={layoutOptions}
          userZoomingEnabled={true}
          userPanningEnabled={true}
          boxSelectionEnabled={false}
        />
      </Box>
    </Box>
  );
};
