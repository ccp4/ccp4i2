import { defineWorkspace } from "vitest/config";
import react from "@vitejs/plugin-react";
import path from "path";
import { fileURLToPath } from "url";

const __dirname = path.dirname(fileURLToPath(import.meta.url));
const alias = { "@": path.resolve(__dirname, "renderer") };

// Two projects so pure-logic tests (schema/format) run on a fast `node`
// environment with no DOM setup, while component tests keep jsdom + the
// next/* mocks in setup.ts. Split is by file extension: *.test.ts = logic,
// *.test.tsx = component.
export default defineWorkspace([
  {
    extends: false,
    resolve: { alias },
    test: {
      name: "node",
      environment: "node",
      include: ["renderer/__tests__/**/*.test.ts"],
    },
  },
  {
    extends: false,
    plugins: [react()],
    resolve: { alias },
    test: {
      name: "dom",
      environment: "jsdom",
      globals: true,
      setupFiles: ["./renderer/__tests__/setup.ts"],
      include: ["renderer/__tests__/**/*.test.tsx"],
      css: false,
      testTimeout: 15000,
    },
  },
]);
