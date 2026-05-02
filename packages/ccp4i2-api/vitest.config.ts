import { defineConfig } from "vitest/config";

export default defineConfig({
  test: {
    // Browser-like globals (window, document) for tests that touch
    // window.ccp4i2LocalSession or the AUTH_ERROR_EVENT dispatch path.
    environment: "jsdom",
    include: ["tests/js/**/*.test.ts"],
    // tsconfig.build.json has rootDir restrictions; the default
    // tsconfig.json applies to test compilation and is more permissive.
  },
});
