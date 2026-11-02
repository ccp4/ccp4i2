/**
 * The scene fetches wait for a server-side build.
 *
 * ``apiGet`` aborts after 30 s unless told otherwise. A campaign summary
 * of a few hundred datasets took longer than that to assemble, so the
 * browser abandoned a response that was on its way and the viewer sat
 * empty. Each scene fetch must pass its own, much longer, ceiling.
 */
import { describe, it, expect, vi, beforeEach } from "vitest";

const apiGet = vi.fn();

vi.mock("../api-fetch", () => ({
  apiGet: (...args: any[]) => apiGet(...args),
  apiPost: vi.fn(),
  apiDelete: vi.fn(),
  apiPatch: vi.fn(),
  apiPut: vi.fn(),
  apiFetch: vi.fn(),
}));

vi.mock("../api", () => ({
  useApi: () => ({ get: vi.fn() }),
}));

import { useCampaignsApi, SCENE_BUILD_TIMEOUT_MS } from "@/lib/campaigns-api";

beforeEach(() => {
  apiGet.mockReset();
  apiGet.mockResolvedValue({ scene: {}, stats: {} });
});

describe("campaign scene fetches", () => {
  it("give the summary scene the scene-build ceiling, not the 30 s default", async () => {
    await useCampaignsApi().fetchSummaryScene(68);
    expect(apiGet).toHaveBeenCalledWith("projectgroups/68/summary_scene", {
      timeout: SCENE_BUILD_TIMEOUT_MS,
    });
    expect(SCENE_BUILD_TIMEOUT_MS).toBeGreaterThan(30_000);
  });

  it("give a site scene the same ceiling", async () => {
    await useCampaignsApi().fetchSiteScene(68, 3, { includeUnclear: true });
    expect(apiGet).toHaveBeenCalledWith(
      "projectgroups/68/sites/3/scene/?include=unclear",
      { timeout: SCENE_BUILD_TIMEOUT_MS },
    );
  });

  it("give a several-sites scene the same ceiling", async () => {
    await useCampaignsApi().fetchSitesScene(68, [3, 5]);
    expect(apiGet).toHaveBeenCalledWith(
      "projectgroups/68/sites/scene/?ids=3%2C5",
      { timeout: SCENE_BUILD_TIMEOUT_MS },
    );
  });
});
