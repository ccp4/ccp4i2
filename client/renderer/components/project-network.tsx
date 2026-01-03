"use client";
import React, { useEffect, useMemo, useRef } from "react";
import { useRouter } from "next/navigation";
//@ts-ignore
import CytoscapeComponent from "react-cytoscapejs";
import Cytoscape from "cytoscape";
import COSEBilkent from "cytoscape-cose-bilkent";
import { useApi } from "../api";
import {
  File as FileInfo,
  FileUse as FileUseInfo,
  Job as JobInfo,
} from "../types/models";
import { FormControlLabel, Radio, RadioGroup } from "@mui/material";
import { useTheme } from "../theme/theme-provider";
// (single useRouter import only)

export interface ProjectNetworkProps {
  projectId: number;
}

Cytoscape.use(COSEBilkent);

export const ProjectNetwork = ({ projectId }: ProjectNetworkProps) => {
  const api = useApi();
  const router = useRouter();
  const { mode } = useTheme();
  const cyRef = useRef<any | null>(null);
  const viewportRef = useRef<
    Record<string, { zoom: number; pan: { x: number; y: number } }>
  >({});

  const [selectedNetwork, setSelectedNetwork] =
    React.useState<string>("fileToFile");
  const { data: fileUses } = api.get<FileUseInfo[]>(
    `projects/${projectId}/file_uses/`
  );
  const { data: files } = api.get<FileInfo[]>(`projects/${projectId}/files/`);
  const { data: jobs } = api.get<JobInfo[]>(`projects/${projectId}/jobs/`);

  const topLevelJobs = useMemo(
    () => (jobs ? jobs.filter((job) => job.parent === null) : []),
    [jobs]
  );

  const topLevelFiles = useMemo(() => {
    const topLevelJobIds = topLevelJobs.map((job) => job.id);
    return files
      ? files.filter((file) => topLevelJobIds.includes(file.job))
      : [];
  }, [files, topLevelJobs]);

  const topLevelFileUses = useMemo(() => {
    const topLevelFileIds = topLevelFiles.map((file) => file.id);
    const topLevelJobIds = topLevelJobs.map((job) => job.id);
    return fileUses
      ? fileUses.filter(
          (fileUse) =>
            topLevelFileIds.includes(fileUse.file) &&
            topLevelJobIds.includes(fileUse.job)
        )
      : [];
  }, [fileUses, topLevelFiles, topLevelJobs]);

  // Helper function to find if there's a path between two jobs through intermediate jobs
  const hasIndirectPath = (
    sourceJobId: number,
    targetJobId: number,
    jobToJobMap: Map<number, Set<number>>
  ): boolean => {
    const visited = new Set<number>();
    const queue = [sourceJobId];
    visited.add(sourceJobId);

    while (queue.length > 0) {
      const currentJobId = queue.shift()!;
      const directTargets = jobToJobMap.get(currentJobId) || new Set();

      for (const nextJobId of directTargets) {
        if (nextJobId === targetJobId) {
          return true; // Found indirect path
        }
        if (!visited.has(nextJobId)) {
          visited.add(nextJobId);
          queue.push(nextJobId);
        }
      }
    }
    return false;
  };

  const networkElements = useMemo(() => {
    if (!topLevelFileUses || !topLevelFiles || !topLevelJobs) return [];

    const fileNodes = topLevelFiles.map((file) => {
      const job = topLevelJobs.find((j) => j.id === file.job);
      return {
        data: {
          id: `file-${file.id}`,
          label: `${job?.number}: ${file.annotation || file.job_param_name}`,
          type: "file",
        },
      };
    });

    const jobNodes = topLevelJobs.map((job) => ({
      data: {
        id: `job-${job.id}`,
        label: `${job.number}: ${job.title || job.task_name}`,
        type: "job",
        jobId: job.id,
        projectId: (job as any)?.project ?? projectId,
      },
    }));

    const fileUseEdges = topLevelFileUses.map((fileUse) => ({
      data: {
        id: `file-use-${fileUse.id}`,
        source: `file-${fileUse.file}`,
        target: `job-${fileUse.job}`,
      },
    }));

    const fileFromEdges = topLevelFiles.map((file) => ({
      data: {
        id: `file-${file.id}-from-job-${file.job}`,
        source: `file-${file.id}`,
        target: `job-${file.job}`,
      },
    }));

    // Calculate job-to-job edges
    const jobToJobEdges = topLevelFileUses
      .map((fileUse) => {
        const file = topLevelFiles.find((f) => f.id === fileUse.file);
        if (!file) return null;
        // Only create edge if file.job !== fileUse.job
        if (file.job !== fileUse.job) {
          return {
            data: {
              id: `job-to-job-${fileUse.id}`,
              source: `job-${file.job}`,
              target: `job-${fileUse.job}`,
              type: "job-to-job",
            },
            sourceJobId: file.job,
            targetJobId: fileUse.job,
          };
        }
        return null;
      })
      .filter(Boolean);

    // Calculate pruned job-to-job edges
    const prunedJobToJobEdges = (() => {
      // Group edges by source-target pair to handle multiple file relationships
      const edgeGroups = new Map<string, (typeof jobToJobEdges)[0][]>();

      jobToJobEdges.forEach((edge) => {
        if (!edge) return;
        const key = `${edge.sourceJobId}-${edge.targetJobId}`;
        if (!edgeGroups.has(key)) {
          edgeGroups.set(key, []);
        }
        edgeGroups.get(key)!.push(edge);
      });

      // Create a map of unique job-to-job connections (one per source-target pair)
      const uniqueConnections = new Map<number, Set<number>>();
      Array.from(edgeGroups.keys()).forEach((key) => {
        const [sourceStr, targetStr] = key.split("-");
        const sourceId = parseInt(sourceStr);
        const targetId = parseInt(targetStr);

        if (!uniqueConnections.has(sourceId)) {
          uniqueConnections.set(sourceId, new Set());
        }
        uniqueConnections.get(sourceId)!.add(targetId);
      });

      // Filter out edge groups that have indirect paths through other jobs
      const keptEdgeGroups = Array.from(edgeGroups.entries()).filter(
        ([key, edges]) => {
          const [sourceStr, targetStr] = key.split("-");
          const sourceId = parseInt(sourceStr);
          const targetId = parseInt(targetStr);

          // Create a temporary map without this direct connection
          const tempMap = new Map<number, Set<number>>();
          uniqueConnections.forEach((targets, source) => {
            const filteredTargets = new Set<number>();
            targets.forEach((target) => {
              // Exclude the current direct connection
              if (!(source === sourceId && target === targetId)) {
                filteredTargets.add(target);
              }
            });
            if (filteredTargets.size > 0) {
              tempMap.set(source, filteredTargets);
            }
          });

          // Check if there's an indirect path without this direct connection
          return !hasIndirectPath(sourceId, targetId, tempMap);
        }
      );

      // Return all edges from the kept groups (preserving multiple file relationships)
      return keptEdgeGroups.flatMap(([key, edges]) => edges);
    })();

    // Calculate file-to-file edges
    const fileToFileEdges = topLevelFileUses
      .filter((fileUse) => fileUse.role === 1)
      .map((fileUse) => {
        const job = topLevelJobs.find((j) => j.id === fileUse.job);
        const sourceFile = topLevelFiles.find((f) => f.job === fileUse.job);
        const targetFile = topLevelFiles.find((f) => f.id === fileUse.file);
        // Only create edge if both files exist and are different
        if (sourceFile && targetFile && sourceFile.id !== targetFile.id) {
          return {
            data: {
              id: `file-to-file-${sourceFile.id}-to-${targetFile.id}`,
              source: `file-${sourceFile.id}`,
              target: `file-${targetFile.id}`,
              type: "file-to-file",
              label: job?.task_name,
            },
          };
        }
        return null;
      })
      .filter(Boolean);

    // Compose elements for the selected network
    let elements: any[];
    switch (selectedNetwork) {
      case "fileToFile":
        elements = [...fileNodes, ...(fileToFileEdges as any[])];
        break;
      case "jobToJob":
        elements = [...jobNodes, ...(jobToJobEdges as any[])];
        break;
      case "prunedJobToJob":
        elements = [...jobNodes, ...(prunedJobToJobEdges as any[])];
        break;
      case "full":
      default:
        elements = [
          ...fileNodes,
          ...jobNodes,
          ...fileUseEdges,
          ...fileFromEdges,
        ];
        break;
    }

    // Filter out unconnected nodes unless there are only nodes (no edges)
    const edgeEndpoints = new Set<string>();
    let edgeCount = 0;
    elements.forEach((el) => {
      const s = el?.data?.source;
      const t = el?.data?.target;
      if (s && t) {
        edgeCount += 1;
        edgeEndpoints.add(String(s));
        edgeEndpoints.add(String(t));
      }
    });

    if (edgeCount === 0) {
      // No edges: show all nodes (original behavior)
      return elements;
    }

    const filtered = elements.filter((el) => {
      const hasEndpoints = el?.data?.source && el?.data?.target;
      if (hasEndpoints) return true; // keep edges
      const id = el?.data?.id;
      return id && edgeEndpoints.has(String(id));
    });

    return filtered;
  }, [topLevelFileUses, topLevelFiles, topLevelJobs, selectedNetwork]);

  // Compute the node id at the end of the longest chain to anchor layout at the top
  const layoutRootId = useMemo(() => {
    if (!networkElements || networkElements.length === 0) return undefined;

    type NodeId = string;
    const nodes = new Set<NodeId>();
    const out = new Map<NodeId, Set<NodeId>>();
    const inDeg = new Map<NodeId, number>();

    // Collect nodes and edges
    networkElements.forEach((el: any) => {
      if (el.data?.id && !el.data?.source && !el.data?.target) {
        nodes.add(el.data.id as NodeId);
        if (!out.has(el.data.id)) out.set(el.data.id, new Set());
        if (!inDeg.has(el.data.id)) inDeg.set(el.data.id, 0);
      }
    });
    networkElements.forEach((el: any) => {
      const s = el.data?.source as NodeId | undefined;
      const t = el.data?.target as NodeId | undefined;
      if (s && t) {
        nodes.add(s);
        nodes.add(t);
        if (!out.has(s)) out.set(s, new Set());
        out.get(s)!.add(t);
        if (!inDeg.has(s)) inDeg.set(s, 0);
        inDeg.set(t, (inDeg.get(t) || 0) + 1);
        if (!inDeg.has(t)) inDeg.set(t, 0);
        if (!out.has(t)) out.set(t, new Set());
      }
    });

    if (nodes.size === 0) return undefined;

    // Kahn's algorithm for toposort (detect cycles)
    const inDegCopy = new Map(inDeg);
    const queue: NodeId[] = [];
    nodes.forEach((n) => {
      if ((inDegCopy.get(n) || 0) === 0) queue.push(n);
    });
    const topo: NodeId[] = [];
    while (queue.length) {
      const u = queue.shift()!;
      topo.push(u);
      (out.get(u) || new Set()).forEach((v) => {
        inDegCopy.set(v, (inDegCopy.get(v) || 0) - 1);
        if ((inDegCopy.get(v) || 0) === 0) queue.push(v);
      });
    }

    const hasCycle = topo.length !== nodes.size;

    if (!hasCycle) {
      // Longest path in DAG ending at sinks
      const dp = new Map<NodeId, number>();
      topo.forEach((n) => dp.set(n, 0));
      topo.forEach((u) => {
        (out.get(u) || new Set()).forEach((v) => {
          dp.set(v, Math.max(dp.get(v) || 0, (dp.get(u) || 0) + 1));
        });
      });
      // Find sinks
      let bestNode: NodeId | undefined;
      let bestDist = -1;
      nodes.forEach((n) => {
        const isSink = (out.get(n)?.size || 0) === 0;
        const d = dp.get(n) || 0;
        if (isSink && d >= bestDist) {
          bestDist = d;
          bestNode = n;
        }
      });
      // Fallback to max dp if no explicit sink (e.g., isolated node)
      if (!bestNode) {
        topo.forEach((n) => {
          const d = dp.get(n) || 0;
          if (d >= bestDist) {
            bestDist = d;
            bestNode = n;
          }
        });
      }
      return bestNode;
    }

    // Cycle fallback: use tree-diameter heuristic on undirected view
    const undirected = new Map<NodeId, Set<NodeId>>();
    nodes.forEach((n) => undirected.set(n, new Set()));
    networkElements.forEach((el: any) => {
      const s = el.data?.source as NodeId | undefined;
      const t = el.data?.target as NodeId | undefined;
      if (s && t) {
        undirected.get(s)!.add(t);
        undirected.get(t)!.add(s);
      }
    });
    const bfsFarthest = (start: NodeId) => {
      const q: NodeId[] = [start];
      const dist = new Map<NodeId, number>([[start, 0]]);
      let far: NodeId = start;
      while (q.length) {
        const u = q.shift()!;
        const du = dist.get(u) || 0;
        (undirected.get(u) || new Set()).forEach((v) => {
          if (!dist.has(v)) {
            dist.set(v, du + 1);
            q.push(v);
            if (du + 1 >= (dist.get(far) || 0)) far = v;
          }
        });
      }
      return far;
    };
    const any = nodes.values().next().value as NodeId;
    const a = bfsFarthest(any);
    const b = bfsFarthest(a);
    return b;
  }, [networkElements]);

  // Cytoscape stylesheet for directional arrows and labels
  const cytoscapeStyles: any = useMemo(
    () => [
      {
        selector: "node",
        style: {
          label: "data(label)",
          "text-wrap": "wrap",
          "text-max-width": 180,
          "font-size": 10,
          "background-color": mode === "dark" ? "#424242" : "#e0e0e0",
          "border-width": 1,
          "border-color": mode === "dark" ? "#757575" : "#999",
          color: mode === "dark" ? "#ffffff" : "#000000",
        },
      },
      {
        selector: 'node[type = "job"]',
        style: {
          shape: "round-rectangle",
          "background-color": mode === "dark" ? "#1e3a5f" : "#d1eaff",
          "border-color": mode === "dark" ? "#4fc3f7" : "#5aa6e8",
        },
      },
      {
        selector: 'node[type = "file"]',
        style: {
          shape: "ellipse",
          "background-color": mode === "dark" ? "#4a2c17" : "#ffe7cc",
          "border-color": mode === "dark" ? "#ffb74d" : "#ffb366",
        },
      },
      {
        selector: "edge",
        style: {
          width: 2,
          "curve-style": "bezier",
          "line-color": mode === "dark" ? "#757575" : "#999",
          // no arrows by default; scoped by type rules below
          "source-arrow-shape": "none",
          "target-arrow-shape": "none",
          "arrow-scale": 1,
          label: "data(label)",
          "font-size": 9,
          "text-rotation": "autorotate",
          "text-margin-y": -6,
          color: mode === "dark" ? "#ffffff" : "#000000",
        },
      },
      {
        selector: 'edge[type = "job-to-job"]',
        style: {
          "line-color": mode === "dark" ? "#9575cd" : "#6a5acd",
          // arrow at the descendent end (target)
          "source-arrow-shape": "none",
          "target-arrow-shape": "triangle",
          "target-arrow-color": mode === "dark" ? "#9575cd" : "#6a5acd",
        },
      },
      {
        selector: 'edge[type = "file-to-file"]',
        style: {
          "line-color": mode === "dark" ? "#4db6ac" : "#2e8b57",
          // arrow at the antecedent end (source)
          "target-arrow-shape": "none",
          "source-arrow-shape": "triangle",
          "source-arrow-color": mode === "dark" ? "#4db6ac" : "#2e8b57",
        },
      },
    ],
    [mode]
  );

  // Layout options: compact for all except "full" view
  const layoutOptions = useMemo(() => {
    const isFull = selectedNetwork === "full";
    return {
      name: "breadthfirst",
      animate: false,
      directed: false,
      spacingFactor: isFull ? 1.1 : 0.5,
      fit: false, // we'll manage viewport separately for stability
      padding: 10,
      roots: layoutRootId ? `#${layoutRootId}` : undefined,
    } as any;
  }, [selectedNetwork, layoutRootId]);

  // Persist current viewport whenever user changes it
  useEffect(() => {
    const cy = cyRef.current;
    if (!cy) return;
    const handler = () => {
      viewportRef.current[selectedNetwork] = {
        zoom: cy.zoom(),
        pan: cy.pan(),
      };
    };
    cy.on("zoom pan", handler);
    return () => {
      if (cy) cy.off("zoom pan", handler);
    };
  }, [selectedNetwork]);

  // Run layout and restore a stable viewport per view
  useEffect(() => {
    const cy = cyRef.current;
    if (!cy) return;
    // Run layout without changing viewport
    const layout = cy.layout(layoutOptions as any);
    layout.run();
    layout.once("layoutstop", () => {
      const saved = viewportRef.current[selectedNetwork];
      if (saved) {
        cy.zoom(saved.zoom);
        cy.pan(saved.pan);
      } else {
        // Initial fit once per view, then persist that viewport
        cy.fit(undefined, 20);
        viewportRef.current[selectedNetwork] = {
          zoom: cy.zoom(),
          pan: cy.pan(),
        };
      }
    });
  }, [networkElements, layoutOptions, selectedNetwork]);

  // Click vs long-press navigation for job nodes in job-to-job views
  useEffect(() => {
    const cy = cyRef.current;
    if (!cy) return;

    let tapStartTime = 0;
    let dragged = false;
    const CLICK_MS = 250; // threshold to consider as a click

    const shouldHandle = () =>
      selectedNetwork === "jobToJob" || selectedNetwork === "prunedJobToJob";

    const onTapStart = (evt: any) => {
      if (!shouldHandle()) return;
      const ele = evt.target;
      if (!ele || ele.group?.() !== "nodes") return;
      if (ele.data?.("type") !== "job") return;
      tapStartTime = performance.now();
      dragged = false;
    };

    const onDrag = (evt: any) => {
      if (!shouldHandle()) return;
      const ele = evt.target;
      if (!ele || ele.group?.() !== "nodes") return;
      if (ele.data?.("type") !== "job") return;
      dragged = true;
    };

    const onTapEnd = (evt: any) => {
      if (!shouldHandle()) return;
      const ele = evt.target;
      if (!ele || ele.group?.() !== "nodes") return;
      if (ele.data?.("type") !== "job") return;
      const dt = performance.now() - tapStartTime;
      if (!dragged && dt < CLICK_MS) {
        const jobId = ele.data("jobId");
        const projId = ele.data("projectId") ?? projectId;
        if (jobId && projId) {
          router.push(`/project/${projId}/job/${jobId}`);
        }
      }
    };

    cy.on("tapstart", 'node[type = "job"]', onTapStart);
    cy.on("drag", 'node[type = "job"]', onDrag);
    cy.on("tapend", 'node[type = "job"]', onTapEnd);

    return () => {
      cy.off("tapstart", 'node[type = "job"]', onTapStart);
      cy.off("drag", 'node[type = "job"]', onDrag);
      cy.off("tapend", 'node[type = "job"]', onTapEnd);
    };
  }, [selectedNetwork, router, projectId]);

  return (
    <>
      <div style={{ marginBottom: 16 }}>
        <RadioGroup
          row
          value={selectedNetwork}
          onChange={(e) => setSelectedNetwork(e.target.value)}
        >
          <FormControlLabel
            value="fileToFile"
            control={<Radio />}
            label="File-to-File"
          />
          <FormControlLabel
            value="jobToJob"
            control={<Radio />}
            label="Job-to-Job"
          />
          <FormControlLabel
            value="prunedJobToJob"
            control={<Radio />}
            label="Pruned Job-to-Job"
          />
          <FormControlLabel
            value="full"
            control={<Radio />}
            label="Full Network"
          />
        </RadioGroup>
      </div>
      {networkElements.length > 0 ? (
        <CytoscapeComponent
          elements={networkElements}
          stylesheet={cytoscapeStyles}
          style={{
            width: "1200px",
            height: "1200px",
            backgroundColor: mode === "dark" ? "#121212" : "#ffffff",
          }}
          layout={layoutOptions}
          cy={(cy) => {
            // Save cy instance and restore viewport if available
            cyRef.current = cy;
            const saved = viewportRef.current[selectedNetwork];
            if (saved) {
              cy.zoom(saved.zoom);
              cy.pan(saved.pan);
            }
          }}
        />
      ) : (
        <p>No data to display</p>
      )}
    </>
  );
};
