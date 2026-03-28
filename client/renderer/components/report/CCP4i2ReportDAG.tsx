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
import { useEffect, useMemo, useRef, useState } from "react";
import $ from "jquery";
import { Box, Typography, CircularProgress } from "@mui/material";
import { useTheme } from "@mui/material/styles";
import { CCP4i2ReportElementProps } from "./CCP4i2ReportElement";

/** Dark-mode colour swaps for Airlie's named/hex node colours. */
const DARK_COLOR_SWAPS: Record<string, string> = {
  green: "#4caf50",
  gold: "#ffd54f",
  grey: "#9e9e9e",
  darkkhaki: "#bdb76b",
  orange: "#ff9800",
  skyblue: "#4fc3f7",
  // R-factor gradient — lighten for dark backgrounds
  "#0000ff": "#4444ff",
  "#2b00d5": "#5533ff",
  "#5500aa": "#7733cc",
  "#800080": "#aa44aa",
  "#aa0055": "#cc4477",
  "#d5002a": "#ee3355",
  "#ff0000": "#ff4444",
};

interface VisNode {
  id: number | string;
  label: string;
  color: string | { background: string; border: string };
  shape: string;
  value?: number;
  borderWidth?: number;
}

interface VisEdge {
  from: number | string;
  to: number | string;
}

interface VisData {
  nodes: VisNode[];
  edges: VisEdge[];
}

export const CCP4i2ReportDAG: React.FC<CCP4i2ReportElementProps> = ({
  item,
}) => {
  const theme = useTheme();
  const isDark = theme.palette.mode === "dark";
  const containerRef = useRef<HTMLDivElement>(null);
  const networkRef = useRef<any>(null);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);

  const { graphData, title, isLegacy } = useMemo(() => {
    const titleAttr = $(item).attr("title") || "Pipeline DAG";
    const text = $(item).text().trim();

    // DOT strings from old Viz.js format
    if (text.startsWith("digraph")) {
      return { graphData: null, title: titleAttr, isLegacy: true };
    }

    // Try parsing as vis-network JSON
    try {
      const parsed = JSON.parse(text) as VisData;
      if (parsed.nodes && parsed.edges) {
        return { graphData: parsed, title: titleAttr, isLegacy: false };
      }
    } catch {
      // Not JSON
    }

    // Old Cytoscape JSON format
    if (text.startsWith("[") || text.startsWith("{")) {
      return { graphData: null, title: titleAttr, isLegacy: true };
    }

    return { graphData: null, title: titleAttr, isLegacy: false };
  }, [item]);

  // Apply dark mode colour swaps
  const themedData = useMemo(() => {
    if (!graphData) return null;
    if (!isDark) return graphData;

    return {
      ...graphData,
      nodes: graphData.nodes.map((node) => {
        let color: any;
        if (typeof node.color === "string") {
          color = DARK_COLOR_SWAPS[node.color] || node.color;
        } else {
          // Legacy {background, border} format
          color = node.color;
        }
        return {
          ...node,
          color,
          font: { color: "#e0e0e0" },
        };
      }),
    };
  }, [graphData, isDark]);

  // Initialize vis-network
  useEffect(() => {
    if (!themedData || !containerRef.current) return;

    let destroyed = false;

    (async () => {
      try {
        const { Network } = await import("vis-network/standalone");

        if (destroyed || !containerRef.current) return;

        // Options matching Airlie's runTree.py pyvis_display_method (basic)
        const options = {
          layout: {
            hierarchical: {
              enabled: true,
              sortMethod: "hubsize" as const,
              direction: "UD" as const,
              levelSeparation: 115,
              nodeSpacing: 120,
              blockShifting: true,
              edgeMinimization: true,
              parentCentralization: true,
              shakeTowards: "roots" as const,
            },
          },
          nodes: {
            font: {
              size: 31,
              face: "Helvetica, Arial, sans-serif",
              color: isDark ? "#e0e0e0" : "#000000",
            },
            scaling: {
              min: 10,
              max: 60,
              label: { enabled: true, min: 14, max: 31 },
            },
          },
          edges: {
            arrows: { to: { enabled: true } },
            color: { inherit: true },
            smooth: false,
          },
          physics: {
            enabled: false,
          },
          interaction: {
            dragNodes: false,
            zoomView: false,
            dragView: true,
          },
        };

        const network = new Network(
          containerRef.current,
          { nodes: themedData.nodes, edges: themedData.edges },
          options
        );

        networkRef.current = network;
        network.once("afterDrawing", () => {
          if (!destroyed) {
            setLoading(false);
            network.fit({ animation: false });
          }
        });
      } catch (err: any) {
        if (!destroyed) {
          setError(err.message || "Failed to render graph");
          setLoading(false);
        }
      }
    })();

    return () => {
      destroyed = true;
      if (networkRef.current) {
        networkRef.current.destroy();
        networkRef.current = null;
      }
    };
  }, [themedData, isDark]);

  if (isLegacy) {
    return (
      <Box sx={{ my: 1 }}>
        <Typography variant="body2" color="text.secondary">
          Legacy DAG format — regenerate the report to see the updated
          visualization.
        </Typography>
      </Box>
    );
  }

  if (!graphData) return null;

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
          height: "350px",
          position: "relative",
        }}
      >
        {loading && (
          <Box
            sx={{
              position: "absolute",
              inset: 0,
              display: "flex",
              alignItems: "center",
              justifyContent: "center",
              zIndex: 1,
            }}
          >
            <CircularProgress size={20} />
            <Typography variant="body2" sx={{ ml: 1 }}>
              Loading graph…
            </Typography>
          </Box>
        )}
        {error ? (
          <Box sx={{ p: 2 }}>
            <Typography color="error" variant="body2">
              Error rendering graph: {error}
            </Typography>
          </Box>
        ) : (
          <div
            ref={containerRef}
            style={{ width: "100%", height: "100%" }}
          />
        )}
      </Box>
    </Box>
  );
};
