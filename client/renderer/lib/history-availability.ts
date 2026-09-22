import { usePathname } from "next/navigation";
import { useEffect, useState } from "react";

interface HistoryAvailability {
  canGoBack: boolean;
  canGoForward: boolean;
}

/**
 * Whether this window has anywhere to go back or forward to.
 *
 * Needed because the app opens plenty of pages with window.open, which
 * Electron turns into a fresh BrowserWindow with no session history at all —
 * a back arrow there is simply dead, and the page you came from is still on
 * screen behind it.
 *
 * Chromium's Navigation API answers both questions directly, which covers
 * Electron and the Chromium browsers. Elsewhere we can only tell that a window
 * has more than one entry, not which side of it we are on, so forward stays
 * enabled: a button that does nothing is a smaller fault than one that is
 * greyed out over a navigation the user could really have made.
 *
 * Both start false because at server-render time there is no history to ask
 * about, and a control that does nothing should not be offered in the meantime.
 */
export function useHistoryAvailability(): HistoryAvailability {
  const pathname = usePathname();
  const [availability, setAvailability] = useState<HistoryAvailability>({
    canGoBack: false,
    canGoForward: false,
  });

  useEffect(() => {
    const navigation = (window as any).navigation;

    const update = () =>
      setAvailability(
        navigation
          ? {
              canGoBack: !!navigation.canGoBack,
              canGoForward: !!navigation.canGoForward,
            }
          : { canGoBack: window.history.length > 1, canGoForward: true }
      );

    // currententrychange fires synchronously inside history.pushState, and the
    // app router calls pushState from a useInsertionEffect — so updating state
    // straight from the listener draws "useInsertionEffect must not schedule
    // updates" from React. A microtask puts it back on ordinary footing.
    const schedule = () => queueMicrotask(update);

    update();
    navigation?.addEventListener("currententrychange", schedule);
    window.addEventListener("popstate", schedule);
    return () => {
      navigation?.removeEventListener("currententrychange", schedule);
      window.removeEventListener("popstate", schedule);
    };
  }, [pathname]);

  return availability;
}
