"use-client";
import { useRef, useEffect } from "react";

interface LatticeViewerProps {
  url: URL | null;
}
const LatticeViewer: React.FC<LatticeViewerProps> = ({ url }) => {
  const canvasRef = useRef<HTMLCanvasElement | null>(null);

  useEffect(() => {
    const resizeCanvas = () => {
      const canvas = canvasRef.current;
      if (!canvas) return;

      const parent = canvas.parentElement;
      if (!parent) return;

      canvas.width = parent.clientWidth;
      canvas.height = parent.clientHeight;

      const ctx = canvas.getContext("2d");
      if (!ctx) return;

      // Clear canvas
      ctx.clearRect(0, 0, canvas.width, canvas.height);

      // Example drawing: A responsive rectangle
      ctx.fillStyle = "#3498db";
      ctx.fillRect(
        canvas.width / 4,
        canvas.height / 4,
        canvas.width / 2,
        canvas.height / 2
      );
    };

    resizeCanvas();
    window.addEventListener("resize", resizeCanvas);
    return () => window.removeEventListener("resize", resizeCanvas);
  }, []);

  return (
    <div className="flex justify-center items-center h-screen">
      <canvas ref={canvasRef} className="border border-gray-300" />
    </div>
  );
};

export default LatticeViewer;

function degreesToRadians(degrees: number): number {
  return (degrees * Math.PI) / 180;
}

function crossProduct(v1: number[], v2: number[]): number[] {
  return [
    v1[1] * v2[2] - v1[2] * v2[1],
    v1[2] * v2[0] - v1[0] * v2[2],
    v1[0] * v2[1] - v1[1] * v2[0],
  ];
}

function vectorMagnitude(v: number[]): number {
  return Math.sqrt(v[0] ** 2 + v[1] ** 2 + v[2] ** 2);
}

function reciprocalLatticeVectors(
  a: number,
  b: number,
  c: number,
  alpha: number,
  beta: number,
  gamma: number
): { aStar: number[]; bStar: number[]; cStar: number[] } {
  // Convert angles to radians
  alpha = degreesToRadians(alpha);
  beta = degreesToRadians(beta);
  gamma = degreesToRadians(gamma);

  // Direct lattice basis vectors
  let aVec: number[] = [a, 0, 0];
  let bVec: number[] = [b * Math.cos(gamma), b * Math.sin(gamma), 0];
  let cVec: number[] = [
    c * Math.cos(beta),
    (c * (Math.cos(alpha) - Math.cos(beta) * Math.cos(gamma))) /
      Math.sin(gamma),
    (c *
      Math.sqrt(
        1 -
          Math.cos(alpha) ** 2 -
          Math.cos(beta) ** 2 -
          Math.cos(gamma) ** 2 +
          2 * Math.cos(alpha) * Math.cos(beta) * Math.cos(gamma)
      )) /
      Math.sin(gamma),
  ];

  // Unit cell volume
  let volume: number =
    aVec[0] * (bVec[1] * cVec[2] - bVec[2] * cVec[1]) -
    aVec[1] * (bVec[0] * cVec[2] - bVec[2] * cVec[0]) +
    aVec[2] * (bVec[0] * cVec[1] - bVec[1] * cVec[0]);

  // Reciprocal lattice vectors
  let aStar: number[] = crossProduct(bVec, cVec).map((v) => v / volume);
  let bStar: number[] = crossProduct(cVec, aVec).map((v) => v / volume);
  let cStar: number[] = crossProduct(aVec, bVec).map((v) => v / volume);

  return { aStar, bStar, cStar };
}

// Example Usage
let reciprocalVectors = reciprocalLatticeVectors(3, 4, 5, 90, 90, 90);
console.log("Reciprocal Lattice Vectors:", reciprocalVectors);
