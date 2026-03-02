"use client";
import { useRef, useEffect, useCallback, useState } from "react";
import { useTheme } from "@mui/material/styles";
import type { SectionPoint, SectionBasis, ReciprocalVectors } from "../../lib/reciprocal-lattice";

interface ReciprocalLatticeCanvasProps {
  sectionData: SectionPoint[];
  uLabel: string;
  vLabel: string;
  showIntensity: boolean;
  colorRange?: [number, number];
  /** Projected reciprocal-space basis for grid drawing */
  sectionBasis: SectionBasis | null;
  /** Resolution limits [low, high] in Angstroms */
  resolution?: [number, number];
  /** Reciprocal lattice vectors for d-spacing computation */
  recipVecs: ReciprocalVectors | null;
  /** Label for the intensity column (shown in tooltip) */
  intensityLabel?: string;
}

/**
 * Map an intensity value to a normalized [0,1] parameter via log scaling,
 * clamped to the given range.
 */
function intensityToT(value: number, min: number, max: number): number {
  const logMin = min > 0 ? Math.log(min) : 0;
  const logMax = max > 0 ? Math.log(max) : 1;
  const logVal = value > 0 ? Math.log(value) : logMin;
  if (logMax <= logMin) return 0.5;
  return Math.max(0, Math.min(1, (logVal - logMin) / (logMax - logMin)));
}

/** Compute d-spacing from Miller indices and reciprocal vectors */
function computeDSpacing(h: number, k: number, l: number, rv: ReciprocalVectors): number {
  const sx = h * rv.aStar[0] + k * rv.bStar[0] + l * rv.cStar[0];
  const sy = h * rv.aStar[1] + k * rv.bStar[1] + l * rv.cStar[1];
  const sz = h * rv.aStar[2] + k * rv.bStar[2] + l * rv.cStar[2];
  const sMag = Math.sqrt(sx * sx + sy * sy + sz * sz);
  return sMag > 0 ? 1 / sMag : Infinity;
}

/** Format a number for tooltip display */
function fmtTooltip(v: number): string {
  const abs = Math.abs(v);
  if (abs >= 10000) return v.toExponential(2);
  if (abs >= 1) return v.toFixed(1);
  if (abs >= 0.01) return v.toFixed(3);
  return v.toExponential(2);
}

/** Standard d-spacings for resolution circles (Angstroms) */
const RESOLUTION_RINGS = [20, 10, 5, 4, 3, 2.5, 2, 1.5, 1.2, 1.0];

function render(
  ctx: CanvasRenderingContext2D,
  width: number,
  height: number,
  points: SectionPoint[],
  offsetX: number,
  offsetY: number,
  scale: number,
  isDark: boolean,
  showIntensity: boolean,
  colorRange: [number, number] | undefined,
  uLabel: string,
  vLabel: string,
  basis: SectionBasis | null,
  resolution: [number, number] | undefined,
  hoveredPoint: SectionPoint | null,
  recipVecs: ReciprocalVectors | null,
  intensityLabel: string | undefined
) {
  const bgColor = isDark ? "#121212" : "#ffffff";
  const gridColor = isDark ? "rgba(255,255,255,0.08)" : "rgba(0,0,0,0.06)";
  const ringColor = isDark ? "rgba(255,255,255,0.15)" : "rgba(0,0,0,0.12)";
  const ringLabelColor = isDark ? "rgba(255,255,255,0.35)" : "rgba(0,0,0,0.3)";
  const axisColor = isDark ? "#555555" : "#bbbbbb";
  const textColor = isDark ? "#aaaaaa" : "#555555";

  ctx.fillStyle = bgColor;
  ctx.fillRect(0, 0, width, height);

  // Transform: canvas center + user offset, scale = pixels per Å⁻¹
  const centerX = width / 2 + offsetX;
  const centerY = height / 2 + offsetY;
  const toCanvasX = (u: number) => centerX + u * scale;
  const toCanvasY = (v: number) => centerY - v * scale;

  // Draw resolution circles
  const resLow = resolution?.[0] ?? 50;
  const resHigh = resolution?.[1] ?? 1;
  ctx.setLineDash([4, 4]);
  ctx.lineWidth = 0.8;
  ctx.font = "11px sans-serif";
  ctx.textAlign = "left";

  for (const d of RESOLUTION_RINGS) {
    if (d < resHigh * 0.8 || d > resLow * 1.2) continue;
    const radius = (1 / d) * scale; // 1/d in Å⁻¹, scaled to pixels
    if (radius < 10 || radius > Math.max(width, height) * 2) continue;

    ctx.strokeStyle = ringColor;
    ctx.beginPath();
    ctx.arc(centerX, centerY, radius, 0, Math.PI * 2);
    ctx.stroke();

    // Label the ring
    ctx.fillStyle = ringLabelColor;
    const labelAngle = -Math.PI / 4;
    const lx = centerX + radius * Math.cos(labelAngle);
    const ly = centerY + radius * Math.sin(labelAngle);
    ctx.fillText(`${d}\u00C5`, lx + 3, ly - 3);
  }
  ctx.setLineDash([]);

  // Draw reciprocal lattice grid lines (parallelogram grid)
  if (basis) {
    const { uProj, vProj } = basis;

    const det = uProj[0] * vProj[1] - uProj[1] * vProj[0];
    if (Math.abs(det) > 1e-12) {
      const invU: [number, number] = [vProj[1] / det, -vProj[0] / det];
      const invV: [number, number] = [-uProj[1] / det, uProj[0] / det];

      const corners = [
        [-centerX / scale, (centerY) / scale],
        [(width - centerX) / scale, (centerY) / scale],
        [(width - centerX) / scale, (centerY - height) / scale],
        [-centerX / scale, (centerY - height) / scale],
      ];

      let minH = Infinity, maxH = -Infinity;
      let minK = Infinity, maxK = -Infinity;
      for (const [cu, cv] of corners) {
        const hf = cu * invU[0] + cv * invU[1];
        const kf = cu * invV[0] + cv * invV[1];
        if (hf < minH) minH = hf;
        if (hf > maxH) maxH = hf;
        if (kf < minK) minK = kf;
        if (kf > maxK) maxK = kf;
      }

      const hStart = Math.floor(minH) - 1;
      const hEnd = Math.ceil(maxH) + 1;
      const kStart = Math.floor(minK) - 1;
      const kEnd = Math.ceil(maxK) + 1;

      const maxLines = 100;

      ctx.strokeStyle = gridColor;
      ctx.lineWidth = 0.5;

      if (hEnd - hStart < maxLines) {
        for (let h = hStart; h <= hEnd; h++) {
          const x1 = toCanvasX(h * uProj[0] + kStart * vProj[0]);
          const y1 = toCanvasY(h * uProj[1] + kStart * vProj[1]);
          const x2 = toCanvasX(h * uProj[0] + kEnd * vProj[0]);
          const y2 = toCanvasY(h * uProj[1] + kEnd * vProj[1]);
          ctx.beginPath();
          ctx.moveTo(x1, y1);
          ctx.lineTo(x2, y2);
          ctx.stroke();
        }
      }

      if (kEnd - kStart < maxLines) {
        for (let k = kStart; k <= kEnd; k++) {
          const x1 = toCanvasX(hStart * uProj[0] + k * vProj[0]);
          const y1 = toCanvasY(hStart * uProj[1] + k * vProj[1]);
          const x2 = toCanvasX(hEnd * uProj[0] + k * vProj[0]);
          const y2 = toCanvasY(hEnd * uProj[1] + k * vProj[1]);
          ctx.beginPath();
          ctx.moveTo(x1, y1);
          ctx.lineTo(x2, y2);
          ctx.stroke();
        }
      }
    }
  }

  // Draw origin crosshair
  ctx.strokeStyle = axisColor;
  ctx.lineWidth = 1;
  ctx.beginPath();
  ctx.moveTo(centerX, 0);
  ctx.lineTo(centerX, height);
  ctx.stroke();
  ctx.beginPath();
  ctx.moveTo(0, centerY);
  ctx.lineTo(width, centerY);
  ctx.stroke();

  // Draw basis vector arrows from origin
  if (basis) {
    const { uProj, vProj } = basis;

    const uEndX = toCanvasX(uProj[0]);
    const uEndY = toCanvasY(uProj[1]);
    ctx.strokeStyle = isDark ? "#ff8a65" : "#e65100";
    ctx.lineWidth = 1.5;
    ctx.beginPath();
    ctx.moveTo(centerX, centerY);
    ctx.lineTo(uEndX, uEndY);
    ctx.stroke();

    const vEndX = toCanvasX(vProj[0]);
    const vEndY = toCanvasY(vProj[1]);
    ctx.strokeStyle = isDark ? "#81c784" : "#2e7d32";
    ctx.beginPath();
    ctx.moveTo(centerX, centerY);
    ctx.lineTo(vEndX, vEndY);
    ctx.stroke();

    ctx.font = "bold 13px sans-serif";
    ctx.textAlign = "center";
    ctx.fillStyle = isDark ? "#ff8a65" : "#e65100";
    ctx.fillText(uLabel, uEndX + (uEndX > centerX ? 12 : -12), uEndY - 5);
    ctx.fillStyle = isDark ? "#81c784" : "#2e7d32";
    ctx.fillText(vLabel, vEndX + 5, vEndY + (vEndY < centerY ? -8 : 15));
  }

  // Draw reflection points — size + opacity encoding for intensity
  const baseRadius = Math.max(1.5, Math.min(4, scale * 0.003));
  const intMin = colorRange?.[0] ?? 0;
  const intMax = colorRange?.[1] ?? 1;
  // Spot colours: ASU = blue, symmetry-expanded = amber/orange
  const asuR = isDark ? 100 : 21;
  const asuG = isDark ? 181 : 101;
  const asuB = isDark ? 246 : 192;
  const symR = isDark ? 255 : 230;
  const symG = isDark ? 183 : 126;
  const symB = isDark ? 77 : 34;
  const asuDotColor = isDark ? "#64b5f6" : "#1565c0";
  const symDotColor = isDark ? "#ffb74d" : "#e65100";
  const rMin = baseRadius * 0.6;
  const rMax = baseRadius * 2.5;
  const alphaMin = 0.2;
  const alphaMax = 1.0;

  // Draw symmetry-expanded spots first (behind), then ASU spots on top
  for (let pass = 0; pass < 2; pass++) {
    const drawASU = pass === 1;
    for (const pt of points) {
      const isASU = pt.isASU !== false;
      if (isASU !== drawASU) continue;

      const cx = toCanvasX(pt.u);
      const cy = toCanvasY(pt.v);

      let r = baseRadius;
      let alpha = 1.0;

      if (showIntensity && pt.intensity !== undefined && colorRange) {
        const t = intensityToT(pt.intensity, intMin, intMax);
        r = rMin + (rMax - rMin) * t;
        alpha = alphaMin + (alphaMax - alphaMin) * t;
      }

      if (cx < -r || cx > width + r || cy < -r || cy > height + r) continue;

      const sR = isASU ? asuR : symR;
      const sG = isASU ? asuG : symG;
      const sB = isASU ? asuB : symB;

      if (showIntensity && pt.intensity !== undefined && colorRange) {
        ctx.fillStyle = `rgba(${sR},${sG},${sB},${alpha.toFixed(3)})`;
      } else {
        ctx.fillStyle = isASU ? asuDotColor : symDotColor;
      }
      ctx.beginPath();
      ctx.arc(cx, cy, r, 0, Math.PI * 2);
      ctx.fill();
    }
  }

  // Info text
  ctx.fillStyle = textColor;
  ctx.font = "12px sans-serif";
  ctx.textAlign = "left";
  ctx.fillText(`${points.length} reflections`, 8, height - 8);

  // Draw hovered reflection highlight + tooltip
  if (hoveredPoint) {
    const hx = toCanvasX(hoveredPoint.u);
    const hy = toCanvasY(hoveredPoint.v);

    // Highlight ring
    ctx.strokeStyle = isDark ? "#ffffff" : "#000000";
    ctx.lineWidth = 2;
    ctx.beginPath();
    ctx.arc(hx, hy, baseRadius * 3.5, 0, Math.PI * 2);
    ctx.stroke();

    // Build tooltip lines
    const lines: string[] = [];
    lines.push(`(${hoveredPoint.h}, ${hoveredPoint.k}, ${hoveredPoint.l})`);

    if (recipVecs) {
      const d = computeDSpacing(hoveredPoint.h, hoveredPoint.k, hoveredPoint.l, recipVecs);
      if (d < 1000) {
        lines.push(`d = ${d.toFixed(2)} \u00C5`);
      }
    }

    if (hoveredPoint.intensity !== undefined) {
      const label = intensityLabel ?? "I";
      lines.push(`${label} = ${fmtTooltip(hoveredPoint.intensity)}`);
    }

    if (hoveredPoint.isASU === false) {
      lines.push("(symmetry)");
    }

    // Measure text
    ctx.font = "12px monospace";
    const lineHeight = 16;
    const pad = 6;
    const textWidths = lines.map((l) => ctx.measureText(l).width);
    const boxW = Math.max(...textWidths) + 2 * pad;
    const boxH = lines.length * lineHeight + 2 * pad;

    // Position tooltip near the point, keeping it on-screen
    let bx = hx + 14;
    let by = hy - boxH - 6;
    if (bx + boxW > width - 4) bx = hx - boxW - 14;
    if (by < 4) by = hy + 14;
    if (bx < 4) bx = 4;

    // Draw rounded rectangle background
    const cr = 4;
    ctx.fillStyle = isDark ? "rgba(30,30,30,0.92)" : "rgba(255,255,255,0.92)";
    ctx.strokeStyle = isDark ? "rgba(255,255,255,0.25)" : "rgba(0,0,0,0.15)";
    ctx.lineWidth = 1;
    ctx.beginPath();
    ctx.moveTo(bx + cr, by);
    ctx.lineTo(bx + boxW - cr, by);
    ctx.arcTo(bx + boxW, by, bx + boxW, by + cr, cr);
    ctx.lineTo(bx + boxW, by + boxH - cr);
    ctx.arcTo(bx + boxW, by + boxH, bx + boxW - cr, by + boxH, cr);
    ctx.lineTo(bx + cr, by + boxH);
    ctx.arcTo(bx, by + boxH, bx, by + boxH - cr, cr);
    ctx.lineTo(bx, by + cr);
    ctx.arcTo(bx, by, bx + cr, by, cr);
    ctx.closePath();
    ctx.fill();
    ctx.stroke();

    // Draw text
    ctx.fillStyle = isDark ? "#e0e0e0" : "#333333";
    ctx.textAlign = "left";
    for (let i = 0; i < lines.length; i++) {
      ctx.fillText(lines[i], bx + pad, by + pad + lineHeight * (i + 0.8));
    }
  }
}

export const ReciprocalLatticeCanvas: React.FC<ReciprocalLatticeCanvasProps> = ({
  sectionData,
  uLabel,
  vLabel,
  showIntensity,
  colorRange,
  sectionBasis,
  resolution,
  recipVecs,
  intensityLabel,
}) => {
  const canvasRef = useRef<HTMLCanvasElement | null>(null);
  const containerRef = useRef<HTMLDivElement | null>(null);
  const theme = useTheme();
  const isDark = theme.palette.mode === "dark";

  const [offset, setOffset] = useState({ x: 0, y: 0 });
  const [scale, setScale] = useState(500); // pixels per Å⁻¹
  const [dragging, setDragging] = useState(false);
  const [dragStart, setDragStart] = useState({ x: 0, y: 0 });
  const [canvasSize, setCanvasSize] = useState({ width: 600, height: 400 });
  const [hoveredPoint, setHoveredPoint] = useState<SectionPoint | null>(null);
  const hoveredRef = useRef<SectionPoint | null>(null);

  // Clear hover when section data changes
  useEffect(() => {
    hoveredRef.current = null;
    setHoveredPoint(null);
  }, [sectionData]);

  // Auto-fit when section data changes
  useEffect(() => {
    if (sectionData.length === 0) {
      setOffset({ x: 0, y: 0 });
      setScale(500);
      return;
    }

    let minU = Infinity,
      maxU = -Infinity,
      minV = Infinity,
      maxV = -Infinity;
    for (const pt of sectionData) {
      if (pt.u < minU) minU = pt.u;
      if (pt.u > maxU) maxU = pt.u;
      if (pt.v < minV) minV = pt.v;
      if (pt.v > maxV) maxV = pt.v;
    }

    const rangeU = maxU - minU || 0.1;
    const rangeV = maxV - minV || 0.1;
    const padding = 0.15;

    const scaleX = canvasSize.width / (rangeU * (1 + 2 * padding));
    const scaleY = canvasSize.height / (rangeV * (1 + 2 * padding));
    const fitScale = Math.min(scaleX, scaleY);

    setScale(fitScale);
    setOffset({
      x: -((minU + maxU) / 2) * fitScale,
      y: ((minV + maxV) / 2) * fitScale,
    });
  }, [sectionData, canvasSize]);

  // Resize observer
  useEffect(() => {
    const container = containerRef.current;
    if (!container) return;

    const observer = new ResizeObserver((entries) => {
      for (const entry of entries) {
        const { width, height } = entry.contentRect;
        if (width > 0 && height > 0) {
          setCanvasSize({ width: Math.floor(width), height: Math.floor(height) });
        }
      }
    });
    observer.observe(container);
    return () => observer.disconnect();
  }, []);

  // Render
  useEffect(() => {
    const canvas = canvasRef.current;
    if (!canvas) return;

    canvas.width = canvasSize.width;
    canvas.height = canvasSize.height;

    const ctx = canvas.getContext("2d");
    if (!ctx) return;

    render(
      ctx,
      canvasSize.width,
      canvasSize.height,
      sectionData,
      offset.x,
      offset.y,
      scale,
      isDark,
      showIntensity,
      colorRange,
      uLabel,
      vLabel,
      sectionBasis,
      resolution,
      hoveredPoint,
      recipVecs,
      intensityLabel
    );
  }, [canvasSize, sectionData, offset, scale, isDark, showIntensity, colorRange, uLabel, vLabel, sectionBasis, resolution, hoveredPoint, recipVecs, intensityLabel]);

  const handleMouseDown = useCallback(
    (e: React.MouseEvent) => {
      setDragging(true);
      setDragStart({ x: e.clientX - offset.x, y: e.clientY - offset.y });
      // Clear hover on drag start
      if (hoveredRef.current !== null) {
        hoveredRef.current = null;
        setHoveredPoint(null);
      }
    },
    [offset]
  );

  const handleMouseMove = useCallback(
    (e: React.MouseEvent) => {
      if (dragging) {
        setOffset({
          x: e.clientX - dragStart.x,
          y: e.clientY - dragStart.y,
        });
        return;
      }

      // Hover hit-test
      const rect = canvasRef.current?.getBoundingClientRect();
      if (!rect) return;

      const mouseX = e.clientX - rect.left;
      const mouseY = e.clientY - rect.top;
      const cX = canvasSize.width / 2 + offset.x;
      const cY = canvasSize.height / 2 + offset.y;
      const recipU = (mouseX - cX) / scale;
      const recipV = (cY - mouseY) / scale;

      // Hit radius: 10px in reciprocal space units
      const hitR = 10 / scale;
      const hitR2 = hitR * hitR;

      let nearest: SectionPoint | null = null;
      let nearestD2 = hitR2;

      for (const pt of sectionData) {
        const du = pt.u - recipU;
        const dv = pt.v - recipV;
        const d2 = du * du + dv * dv;
        if (d2 < nearestD2) {
          nearestD2 = d2;
          nearest = pt;
        }
      }

      // Only update state if the hovered point actually changed
      if (nearest !== hoveredRef.current) {
        hoveredRef.current = nearest;
        setHoveredPoint(nearest);
      }
    },
    [dragging, dragStart, canvasSize, offset, scale, sectionData]
  );

  const handleMouseUp = useCallback(() => setDragging(false), []);

  const handleMouseLeave = useCallback(() => {
    setDragging(false);
    if (hoveredRef.current !== null) {
      hoveredRef.current = null;
      setHoveredPoint(null);
    }
  }, []);

  const handleWheel = useCallback(
    (e: React.WheelEvent) => {
      const rect = canvasRef.current?.getBoundingClientRect();
      if (!rect) return;

      const mouseX = e.clientX - rect.left - canvasSize.width / 2;
      const mouseY = e.clientY - rect.top - canvasSize.height / 2;

      const zoomFactor = e.deltaY > 0 ? 0.9 : 1.1;
      const newScale = Math.max(10, Math.min(50000, scale * zoomFactor));
      const ratio = newScale / scale;

      setOffset((prev) => ({
        x: mouseX - (mouseX - prev.x) * ratio,
        y: mouseY - (mouseY - prev.y) * ratio,
      }));
      setScale(newScale);
    },
    [scale, canvasSize]
  );

  // Double-click to reset view
  const handleDoubleClick = useCallback(() => {
    if (sectionData.length === 0) return;

    let minU = Infinity,
      maxU = -Infinity,
      minV = Infinity,
      maxV = -Infinity;
    for (const pt of sectionData) {
      if (pt.u < minU) minU = pt.u;
      if (pt.u > maxU) maxU = pt.u;
      if (pt.v < minV) minV = pt.v;
      if (pt.v > maxV) maxV = pt.v;
    }

    const rangeU = maxU - minU || 0.1;
    const rangeV = maxV - minV || 0.1;
    const padding = 0.15;

    const scaleX = canvasSize.width / (rangeU * (1 + 2 * padding));
    const scaleY = canvasSize.height / (rangeV * (1 + 2 * padding));
    const fitScale = Math.min(scaleX, scaleY);

    setScale(fitScale);
    setOffset({
      x: -((minU + maxU) / 2) * fitScale,
      y: ((minV + maxV) / 2) * fitScale,
    });
  }, [sectionData, canvasSize]);

  const cursor = dragging ? "grabbing" : hoveredPoint ? "crosshair" : "grab";

  return (
    <div
      ref={containerRef}
      style={{ width: "100%", height: "100%", overflow: "hidden" }}
    >
      <canvas
        ref={canvasRef}
        style={{ display: "block", cursor }}
        onMouseDown={handleMouseDown}
        onMouseMove={handleMouseMove}
        onMouseUp={handleMouseUp}
        onMouseLeave={handleMouseLeave}
        onWheel={handleWheel}
        onDoubleClick={handleDoubleClick}
      />
    </div>
  );
};
