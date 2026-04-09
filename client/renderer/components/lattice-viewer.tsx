"use client";
import { useRef, useEffect } from "react";
export {
  degreesToRadians,
  crossProduct,
  vectorMagnitude,
  reciprocalLatticeVectors,
} from "../lib/reciprocal-lattice";

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

      ctx.clearRect(0, 0, canvas.width, canvas.height);

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
