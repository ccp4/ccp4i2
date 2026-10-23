/**
 * The site-view link: one site, or several viewed as one scene.
 */
import { describe, it, expect } from "vitest";
import { parseSiteParam, siteViewHref } from "../lib/site-view-link";

describe("siteViewHref", () => {
  it("names one site as before", () => {
    expect(siteViewHref(68, [12])).toBe("/ccp4i2/moorhen-page/campaign/68?summary=1&site=12");
  });

  it("names several in the order given", () => {
    expect(siteViewHref(68, [15, 12, 13])).toBe(
      "/ccp4i2/moorhen-page/campaign/68?summary=1&site=15,12,13"
    );
  });
});

describe("parseSiteParam", () => {
  it("reads what siteViewHref writes", () => {
    const href = siteViewHref(68, [15, 12, 13]);
    expect(parseSiteParam(new URL(href, "http://x").searchParams.get("site"))).toEqual([15, 12, 13]);
  });

  it("is empty for no parameter, an empty one, or nothing usable", () => {
    expect(parseSiteParam(null)).toEqual([]);
    expect(parseSiteParam(undefined)).toEqual([]);
    expect(parseSiteParam("")).toEqual([]);
    expect(parseSiteParam(",,")).toEqual([]);
    expect(parseSiteParam("x,-3,1.5,0")).toEqual([]);
  });

  it("drops junk and repeats but keeps the order of the rest", () => {
    expect(parseSiteParam("13, 12,x,13,15")).toEqual([13, 12, 15]);
  });

  it("does not take a prefix of a malformed id, as parseInt would", () => {
    expect(parseSiteParam("12abc")).toEqual([]);
  });
});
