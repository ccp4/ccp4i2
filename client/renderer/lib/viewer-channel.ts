/**
 * Send a Moorhen link to a viewer that is already open, rather than opening
 * another.
 *
 * Every "open in Moorhen" link used to open a new window, which then started
 * Moorhen from nothing and loaded the model and maps again: seconds per click,
 * for a click that usually only wants the camera somewhere else. Moorhen needs
 * its page cross-origin isolated (COOP/COEP, for SharedArrayBuffer), which
 * keeps the viewer in a page of its own on the web; it does not need a new page
 * per click.
 *
 * A BroadcastChannel reaches every same-origin page, however it was opened
 * (COOP and noopener cut the opener relationship, not this). The link asks
 * whether a viewer with the same key is open; the first to answer is told the
 * URL, and replaces its own with it. The viewer pages already apply a scene
 * named in their URL whenever it changes (moorhen-wrapper's sceneSource), and
 * reuse molecules they already hold, so the click costs a camera move and
 * whatever is new. No answer within the timeout: open a window as before.
 *
 * The key is the viewer's route, plus the member job for a campaign viewer:
 * that page reads ?job= only when it opens, so a link for another member opens
 * its own window rather than landing over the wrong model.
 */

const CHANNEL = "ccp4i2-moorhen-viewers";
const ANSWER_TIMEOUT_MS = 300;

type ViewerMessage =
  | { type: "probe"; requestId: string; key: string }
  | { type: "available"; requestId: string; viewerId: string }
  | { type: "show"; requestId: string; viewerId: string; href: string }
  | { type: "shown"; requestId: string; viewerId: string };

const randomId = () => Math.random().toString(36).slice(2) + Date.now().toString(36);

/** Which open viewer may take this link: null for a URL that is not a
 *  viewer page. */
export function viewerKey(href: string, origin = "http://localhost"): string | null {
  let url: URL;
  try {
    url = new URL(href, origin);
  } catch {
    return null;
  }
  if (!url.pathname.startsWith("/ccp4i2/moorhen-page/")) return null;
  if (url.pathname.startsWith("/ccp4i2/moorhen-page/campaign/")) {
    return `${url.pathname}?job=${url.searchParams.get("job") ?? ""}`;
  }
  return url.pathname;
}

const openWindow = (href: string) => window.open(href, "_blank", "noopener,noreferrer");

/** Show `href` in an open viewer that can take it, else open a window.
 *  Resolves to which of the two happened. */
export function showInViewer(href: string): Promise<"reused" | "opened"> {
  const key = viewerKey(href, window.location.origin);
  if (!key || typeof BroadcastChannel === "undefined") {
    openWindow(href);
    return Promise.resolve("opened");
  }
  const channel = new BroadcastChannel(CHANNEL);
  const requestId = randomId();
  return new Promise((resolve) => {
    let chosen: string | null = null;
    const finish = (how: "reused" | "opened") => {
      clearTimeout(timer);
      channel.close();
      resolve(how);
    };
    const timer = setTimeout(() => {
      openWindow(href);
      finish("opened");
    }, ANSWER_TIMEOUT_MS);
    channel.onmessage = (event: MessageEvent<ViewerMessage>) => {
      const message = event.data;
      if (!message || message.requestId !== requestId) return;
      if (message.type === "available" && chosen === null) {
        chosen = message.viewerId;
        channel.postMessage({ type: "show", requestId, viewerId: chosen, href } satisfies ViewerMessage);
      } else if (message.type === "shown" && message.viewerId === chosen) {
        finish("reused");
      }
    };
    channel.postMessage({ type: "probe", requestId, key } satisfies ViewerMessage);
  });
}

/** For a viewer page: answer links that this viewer can take, and show them
 *  by replacing the page's URL (`navigate`, e.g. router.replace). Returns the
 *  function that stops listening. */
export function listenAsViewer(navigate: (href: string) => void): () => void {
  if (typeof BroadcastChannel === "undefined") return () => {};
  const channel = new BroadcastChannel(CHANNEL);
  const viewerId = randomId();
  channel.onmessage = (event: MessageEvent<ViewerMessage>) => {
    const message = event.data;
    if (!message) return;
    if (message.type === "probe" && message.key === viewerKey(window.location.href)) {
      channel.postMessage({ type: "available", requestId: message.requestId, viewerId } satisfies ViewerMessage);
    } else if (message.type === "show" && message.viewerId === viewerId) {
      navigate(message.href);
      // Desktop: bring this window forward (main process, "focus-window").
      // The web cannot raise a tab from script; the link says where it went.
      window.electronAPI?.sendMessage?.("focus-window", {});
      window.focus();
      channel.postMessage({ type: "shown", requestId: message.requestId, viewerId } satisfies ViewerMessage);
    }
  };
  return () => channel.close();
}
