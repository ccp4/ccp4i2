import { defineConfig } from "vite";
import { fileURLToPath } from "url";
import path from "path";
const __dirname = path.dirname(__filename);

export default defineConfig({
  build: {
    outDir: "dist",
    target: "node18", // target Node, not browser
    ssr: true, // <== CRUCIAL: treat this as a server-side build
    rollupOptions: {
      input: path.resolve(__dirname, "main/ccp4i2-master.ts"),
      output: {
        entryFileNames: "main.js",
        format: "esm",
      },
      external: [
        "electron",
        "fs",
        "path",
        "os",
        "stubborn-fs",
        "node:process",
        "node:util",
        "node:os",
        "node:net",
        "node:url",
        "node:path",
        "node:fs",
        "node:child_process",
        "node:module",
        "next",
        "electron-store",
        "detect-port",
        "child_process",
        "express",
        "http-proxy-middleware",
      ], // add more if needed
    },
  },
});
