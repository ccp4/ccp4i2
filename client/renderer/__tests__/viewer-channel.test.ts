/**
 * Moorhen links go to a viewer that is already open, when one can take them.
 *
 * Every link used to open a new window, which started Moorhen from nothing
 * and loaded the data again. The channel (lib/viewer-channel) asks open
 * viewers first, and opens a window only when none answers. Exercised over
 * Node's own BroadcastChannel, which delivers between instances in one
 * process as browsers do between pages.
 */
import { afterEach, beforeEach, describe, expect, it, vi } from "vitest";
import { listenAsViewer, showInViewer, viewerKey } from "../lib/viewer-channel";

const ORIGIN = "http://localhost:3000";
const PROJECT_1 = "/ccp4i2/moorhen-page/project/1";

let opened: string[];
let stops: Array<() => void>;

beforeEach(() => {
  opened = [];
  stops = [];
  (globalThis as any).window = {
    location: { origin: ORIGIN, href: `${ORIGIN}${PROJECT_1}` },
    open: (href: string) => {
      opened.push(href);
      return null;
    },
    focus: () => {},
  };
});

afterEach(() => {
  stops.forEach((stop) => stop());
  delete (globalThis as any).window;
});

const viewer = (navigate = vi.fn()) => {
  stops.push(listenAsViewer(navigate));
  return navigate;
};

describe("viewerKey", () => {
  it("is the route for a project viewer, whatever scene it shows", () => {
    expect(viewerKey(`${PROJECT_1}?sceneJob=4&sceneMethod=eventScene`)).toBe(PROJECT_1);
  });

  it("includes the member job for a campaign viewer", () => {
    expect(viewerKey("/ccp4i2/moorhen-page/campaign/7?job=12&sceneJob=4")).toBe(
      "/ccp4i2/moorhen-page/campaign/7?job=12"
    );
  });

  it("is null for a page that is not a viewer", () => {
    expect(viewerKey("/ccp4i2/project/1")).toBeNull();
  });
});

describe("showInViewer", () => {
  const link = `${PROJECT_1}?sceneJob=4&sceneMethod=eventScene&sceneArgs=%7B%7D`;

  it("opens a window when no viewer is open", async () => {
    expect(await showInViewer(link)).toBe("opened");
    expect(opened).toEqual([link]);
  });

  it("hands the link to an open viewer of the same project, and opens nothing", async () => {
    const navigate = viewer();
    expect(await showInViewer(link)).toBe("reused");
    expect(navigate).toHaveBeenCalledWith(link);
    expect(opened).toEqual([]);
  });

  it("leaves a viewer of another project alone", async () => {
    const navigate = viewer();
    const other = "/ccp4i2/moorhen-page/project/2?sceneJob=4&sceneMethod=eventScene";
    expect(await showInViewer(other)).toBe("opened");
    expect(navigate).not.toHaveBeenCalled();
    expect(opened).toEqual([other]);
  });

  it("gives the link to one viewer when two could take it", async () => {
    const first = viewer();
    const second = viewer();
    expect(await showInViewer(link)).toBe("reused");
    expect(first.mock.calls.length + second.mock.calls.length).toBe(1);
  });
});
