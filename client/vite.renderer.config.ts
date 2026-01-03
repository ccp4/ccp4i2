// vite.renderer.config.ts
import { defineConfig } from "vite";
import react from "@vitejs/plugin-react";
import path from "path";

export default defineConfig({
  root: "./src/renderer",
  plugins: [react()],
  resolve: {
    alias: {
      react: path.resolve(__dirname, "node_modules/react"),
      "react-dom": path.resolve(__dirname, "node_modules/react-dom"),
      ReactDOM: path.resolve(__dirname, "node_modules/react-dom"),
    },
  },
  build: {
    outDir: "../../dist/renderer",
    emptyOutDir: true,
  },
});
