import { describe, expect, it } from "vitest";
import { validateScene } from "../lib/scene";
import { isFetchable } from "../lib/moorhen-scene-resolver";
import { extractProjectPathFromUniqueId, projectPathUrl } from "../lib/scene-project-path";

const PATH = "CCP4_JOBS/job_3/pandda2_out/processed_datasets/xtal-0002/xtal-0002-z_map.native.ccp4";

const sceneWith = (ref: Record<string, unknown>) => ({
  scene: "event",
  version: 1,
  files: [{ name: "zmap", kind: "map", ...ref }],
});

describe("projectPath refs", () => {
  it("round-trip through the loader URL that becomes a uniqueId", () => {
    const url = projectPathUrl(7, PATH);
    expect(url).toBe(`/api/proxy/ccp4i2/projects/7/project_file?path=${encodeURIComponent(PATH)}`);
    expect(extractProjectPathFromUniqueId(url)).toBe(PATH);
    expect(extractProjectPathFromUniqueId("/api/proxy/ccp4i2/files/12/download/")).toBeNull();
  });

  it("validate with a project qualifier and are fetchable", () => {
    const { errors } = validateScene(sceneWith({ projectPath: PATH, projectName: "BAZ2B" }));
    expect(errors).toEqual([]);
    expect(isFetchable({ name: "zmap", projectPath: PATH, projectName: "BAZ2B" })).toBe(true);
  });

  it("need a project", () => {
    const { errors } = validateScene(sceneWith({ projectPath: PATH }));
    expect(errors.map((e) => e.message).join()).toMatch(/project qualifier/);
  });

  it("stay inside the project directory", () => {
    for (const bad of ["/etc/passwd", "CCP4_JOBS/../../x.pdb"]) {
      const { errors } = validateScene(sceneWith({ projectPath: bad, projectName: "p" }));
      expect(errors.map((e) => e.message).join()).toMatch(/stays inside it/);
    }
  });

  it("are one source, not one of two", () => {
    const { errors } = validateScene(sceneWith({ projectPath: PATH, fileId: 3, projectName: "p" }));
    expect(errors.map((e) => e.message).join()).toMatch(/exactly one/);
  });
});
