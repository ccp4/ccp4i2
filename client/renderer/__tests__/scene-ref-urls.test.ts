import { beforeEach, describe, expect, it, vi } from "vitest";

const apiGet = vi.fn();
vi.mock("../api-fetch", () => ({ apiGet: (url: string) => apiGet(url) }));

import { resolveProjectPk, sceneRefUrl } from "../lib/scene-ref-urls";

const page = { pk: 4, uuid: "uuid-4", name: "BAZ2B" };

beforeEach(() => apiGet.mockReset());

describe("resolveProjectPk", () => {
  it("is the page's project when the ref names none", async () => {
    expect(await resolveProjectPk({}, page)).toBe(4);
    expect(apiGet).not.toHaveBeenCalled();
  });

  it("is the page's project when the ref names it, without a lookup", async () => {
    expect(await resolveProjectPk({ uuid: "uuid-4" }, page)).toBe(4);
    expect(await resolveProjectPk({ name: "BAZ2B" }, page)).toBe(4);
    expect(apiGet).not.toHaveBeenCalled();
  });

  it("looks another project up by uuid, and never falls back to the page's", async () => {
    apiGet.mockResolvedValueOnce([{ id: 9 }]);
    expect(await resolveProjectPk({ uuid: "uuid-9" }, page)).toBe(9);
    apiGet.mockResolvedValueOnce([]);
    expect(await resolveProjectPk({ uuid: "uuid-none" }, page)).toBeNull();
  });
});

describe("sceneRefUrl", () => {
  it("needs no project for a pdb, fileId or url ref", async () => {
    expect(await sceneRefUrl({ pdb: "1OGU" }, page)).toBe("/api/proxy/pdbe/entry-files/download/1ogu.cif");
    expect(await sceneRefUrl({ fileId: 12 }, page)).toBe("/api/proxy/ccp4i2/files/12/download/");
    expect(await sceneRefUrl({ url: "https://x/y.pdb" }, page)).toBe("https://x/y.pdb");
    expect(apiGet).not.toHaveBeenCalled();
  });

  it("resolves job+param through the project's jobs and files", async () => {
    apiGet
      .mockResolvedValueOnce([{ id: 77, number: "3" }])
      .mockResolvedValueOnce([{ id: 501, job_param_name: "XYZOUT", directory: 1 }]);
    expect(await sceneRefUrl({ job: 3, param: "XYZOUT" }, page)).toBe("/api/proxy/ccp4i2/files/501/download/");
    expect(apiGet.mock.calls.map((c) => c[0])).toEqual(["jobs/?project=4", "files/?job=77"]);
  });

  it("resolves projectPath against the page's project", async () => {
    expect(await sceneRefUrl({ projectPath: "CCP4_JOBS/job_3/a.pdb" }, page)).toBe(
      "/api/proxy/ccp4i2/projects/4/project_file?path=CCP4_JOBS%2Fjob_3%2Fa.pdb",
    );
  });

  it("has no URL for a bundle or inline ref", async () => {
    expect(await sceneRefUrl({}, page)).toBeNull();
  });
});
