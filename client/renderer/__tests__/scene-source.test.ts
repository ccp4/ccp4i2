import { describe, expect, it } from "vitest";
import { sceneSourceFromParams, sceneSourceUrl } from "../lib/scene-source";

describe("scene sources", () => {
  it("round-trip through the Moorhen page URL", () => {
    const url = sceneSourceUrl(4, { jobId: 48, method: "eventScene", kwargs: { dtag: "x1", event_idx: 2 } });
    expect(url).toMatch(/^\/ccp4i2\/moorhen-page\/project\/4\?/);
    const params = new URL(url as string, "http://h").searchParams;
    expect(sceneSourceFromParams(params)).toEqual({ jobId: 48, method: "eventScene", kwargs: { dtag: "x1", event_idx: 2 } });
  });

  it("open on the campaign viewer's member job when given one", () => {
    const url = sceneSourceUrl(4, { jobId: 48, method: "eventScene", kwargs: { dtag: "x1" } }, { campaignId: 1, jobId: 5 });
    expect(url).toMatch(/^\/ccp4i2\/moorhen-page\/campaign\/1\?job=5&sceneJob=48/);
    const params = new URL(url as string, "http://h").searchParams;
    expect(sceneSourceFromParams(params)).toEqual({ jobId: 48, method: "eventScene", kwargs: { dtag: "x1" } });
    // half a target is no target
    expect(sceneSourceUrl(4, { jobId: 48, method: "eventScene", kwargs: {} }, { campaignId: 1, jobId: NaN }))
      .toMatch(/^\/ccp4i2\/moorhen-page\/project\/4\?/);
  });

  it("call only scene methods", () => {
    expect(sceneSourceUrl(4, { jobId: 48, method: "fillDatasetsFromCampaign", kwargs: {} })).toBeNull();
    const params = new URLSearchParams({ sceneJob: "48", sceneMethod: "fillDatasetsFromCampaign" });
    expect(sceneSourceFromParams(params)).toBeNull();
  });

  it("refuse arguments that are not an object", () => {
    const params = new URLSearchParams({ sceneJob: "48", sceneMethod: "eventScene", sceneArgs: "[1]" });
    expect(sceneSourceFromParams(params)?.kwargs).toEqual({});
    params.set("sceneArgs", "{not json");
    expect(sceneSourceFromParams(params)).toBeNull();
  });
});
