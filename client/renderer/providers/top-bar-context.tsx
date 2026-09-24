"use client";
import {
  createContext,
  PropsWithChildren,
  useContext,
  useEffect,
  useMemo,
  useState,
} from "react";

interface TopBarState {
  /** What this page is called, shown in the bar when there is no project. */
  title?: string;
}

interface TopBarContextValue extends TopBarState {
  setTopBar: (state: TopBarState) => void;
}

const TopBarContext = createContext<TopBarContextValue>({ setTopBar: () => {} });

/**
 * Lets a page name itself in the one app bar mounted by the (authed) layout.
 * Pages no longer render their own bar, so this is how "Preferences" or a
 * campaign's name reaches it.
 */
export function TopBarProvider(props: PropsWithChildren) {
  const [state, setTopBar] = useState<TopBarState>({});
  const value = useMemo(() => ({ ...state, setTopBar }), [state]);
  return (
    <TopBarContext.Provider value={value}>
      {props.children}
    </TopBarContext.Provider>
  );
}

export const useTopBarState = () => useContext(TopBarContext);

export function useTopBar({ title }: TopBarState) {
  const { setTopBar } = useTopBarState();
  useEffect(() => {
    setTopBar({ title });
    return () => setTopBar({});
  }, [setTopBar, title]);
}
