import { defineConfig } from "vitest/config";
import react from "@vitejs/plugin-react";
import path from "path";
import { fileURLToPath } from "url";

const __dirname = path.dirname(fileURLToPath(import.meta.url));

export default defineConfig({
  plugins: [react()],
  resolve: {
    alias: {
      "@": path.resolve(__dirname, "renderer"),
    },
  },
  test: {
    environment: "jsdom",
    globals: true,
    setupFiles: ["./renderer/__tests__/setup.ts"],
    include: ["renderer/__tests__/**/*.test.{ts,tsx}"],
    css: false,
    testTimeout: 15000,
  },
});
