import { describe, expect, it } from "vitest";
import { desktopLaunchCommand, isOpenableRoute, parseOpenRoute } from "../../main/ccp4i2-open-route";

describe("parseOpenRoute", () => {
  it("reads the route from the separate-argument form", () => {
    expect(parseOpenRoute(["app", "--open-route", "/ccp4i2/moorhen-page/session/12"]))
      .toBe("/ccp4i2/moorhen-page/session/12");
  });
  it("reads the route from the equals form", () => {
    expect(parseOpenRoute(["app", "--open-route=/ccp4i2/project/3"])).toBe("/ccp4i2/project/3");
  });
  it("returns null when absent", () => {
    expect(parseOpenRoute(["app", "--no-sandbox"])).toBeNull();
    expect(parseOpenRoute(["app", "--open-route"])).toBeNull();
  });
  it("refuses anything that is not one of the app's own routes", () => {
    for (const bad of [
      "https://example.com/ccp4i2/x",
      "//evil/ccp4i2/x",
      "/other/route",
      "/ccp4i2/with space",
      "ccp4i2/relative",
      "/ccp4i2/<script>",
    ]) {
      expect(isOpenableRoute(bad), bad).toBe(false);
      expect(parseOpenRoute(["app", "--open-route", bad]), bad).toBeNull();
    }
  });
  it("accepts a route with a query string", () => {
    expect(isOpenableRoute("/ccp4i2/moorhen-page/job-by-id/5?view=abc%3D")).toBe(true);
  });
});

describe("desktopLaunchCommand", () => {
  it("is the bare executable when packaged", () => {
    expect(desktopLaunchCommand("/Applications/X.app/Contents/MacOS/x", true, "/ignored"))
      .toEqual(["/Applications/X.app/Contents/MacOS/x"]);
  });
  it("adds the app directory in development", () => {
    expect(desktopLaunchCommand("/n/electron", false, "/repo/client"))
      .toEqual(["/n/electron", "/repo/client"]);
  });
});
