import { usePathname } from "next/navigation";
import { useMemo } from "react";
import { useCCP4i2Window } from "../app-context";

/**
 * The project the current route is about, or null when it is about none.
 *
 * CCP4i2Context holds the last project visited and never clears it, which is
 * what job views and dialogs want. Chrome that is shared across every route
 * cannot use it: after opening a project, the preferences page would still
 * offer to export it. So the project id comes from the URL — except in the
 * stand-alone job window, whose URL names only the job, and where JobView
 * puts the owning project into the context for us.
 */
export function useProjectScope(): number | null {
  const pathname = usePathname() ?? "";
  const { projectId } = useCCP4i2Window();
  return useMemo(() => {
    const match = pathname.match(/^\/ccp4i2\/project\/(\d+)/);
    if (match) return Number(match[1]);
    return /^\/ccp4i2\/job\//.test(pathname) ? projectId ?? null : null;
  }, [pathname, projectId]);
}

/** Whether the current route is showing one particular job. */
export function useIsJobRoute(): boolean {
  const pathname = usePathname() ?? "";
  return (
    /^\/ccp4i2\/job\//.test(pathname) ||
    /^\/ccp4i2\/project\/\d+\/job\//.test(pathname)
  );
}
