"use client";
import React, {
  createContext,
  useContext,
  useState,
  useEffect,
  useCallback,
  useMemo,
  ReactNode,
} from "react";
import { ThemeProvider, createTheme } from "@mui/material/styles";
import CssBaseline from "@mui/material/CssBaseline";
import {
  lightPaletteOptions,
  darkPaletteOptions,
  lightCustomColors,
  darkCustomColors,
} from "./palette";
import { typographyOptions } from "./typography";
import { createComponentVariants } from "./variants";

export type ThemeMode = "light" | "dark";

/**
 * Check if running in an iframe (Teams context).
 */
function isRunningInIframe(): boolean {
  if (typeof window === "undefined") return false;
  try {
    return window.self !== window.top;
  } catch {
    return true;
  }
}

// CSS zoom level for Teams compact mode (0.85 = 85% size)
const TEAMS_ZOOM_LEVEL = 0.85;

interface ThemeContextType {
  mode: ThemeMode;
  toggleTheme: () => void;
  setTheme: (mode: ThemeMode) => void;
  customColors: typeof lightCustomColors;
  isCompact: boolean;
}

const ThemeContext = createContext<ThemeContextType | undefined>(undefined);

export const useTheme = () => {
  const context = useContext(ThemeContext);
  if (!context) {
    throw new Error("useTheme must be used within a ThemeProvider");
  }
  return context;
};

interface CCP4i2ThemeProviderProps {
  children: ReactNode;
}

export const CCP4i2ThemeProvider: React.FC<CCP4i2ThemeProviderProps> = ({
  children,
}) => {
  const [mode, setMode] = useState<ThemeMode>("light");
  const [mounted, setMounted] = useState(false);
  const [isCompact, setIsCompact] = useState(false);

  // Load theme preference from localStorage on mount
  // Also detect Teams iframe context for compact mode
  useEffect(() => {
    const savedTheme = localStorage.getItem("ccp4i2-theme") as ThemeMode;
    if (savedTheme && (savedTheme === "light" || savedTheme === "dark")) {
      setMode(savedTheme);
    } else {
      // Check system preference
      const prefersDark = window.matchMedia(
        "(prefers-color-scheme: dark)"
      ).matches;
      setMode(prefersDark ? "dark" : "light");
    }

    // Enable compact mode in Teams iframe context
    if (isRunningInIframe()) {
      setIsCompact(true);
      // Apply CSS zoom to make everything smaller in Teams
      document.documentElement.style.zoom = String(TEAMS_ZOOM_LEVEL);
    }

    setMounted(true);
  }, []);

  // Save theme preference to localStorage (only after initial mount)
  useEffect(() => {
    if (mounted) {
      localStorage.setItem("ccp4i2-theme", mode);
    }
  }, [mode, mounted]);

  // Listen for theme changes from Electron native menu
  useEffect(() => {
    if (typeof window !== "undefined" && window.electronAPI) {
      const handleMessage = (event: any, data: any) => {
        if (
          data.message === "theme-changed" &&
          (data.theme === "light" || data.theme === "dark")
        ) {
          setMode(data.theme);
        }
      };
      window.electronAPI.onMessage("message-from-main", handleMessage);
    }
  }, []);

  const toggleTheme = useCallback(() => {
    setMode((prevMode) => (prevMode === "light" ? "dark" : "light"));
  }, []);

  const setThemeMode = useCallback((newMode: ThemeMode) => {
    setMode(newMode);
  }, []);

  // Create theme based on current mode
  const theme = useMemo(() => {
    const paletteOptions =
      mode === "light" ? lightPaletteOptions : darkPaletteOptions;

    const baseThemeOptions = {
      palette: paletteOptions,
      typography: typographyOptions,
    };

    const baseTheme = createTheme(baseThemeOptions);

    return createTheme({
      ...baseTheme,
      components: createComponentVariants(baseTheme),
    });
  }, [mode]);

  const customColors = mode === "light" ? lightCustomColors : darkCustomColors;

  // Memoize context value to prevent unnecessary re-renders of consumers
  const contextValue = useMemo(() => ({
    mode,
    toggleTheme,
    setTheme: setThemeMode,
    customColors,
    isCompact,
  }), [mode, toggleTheme, setThemeMode, customColors, isCompact]);

  return (
    <ThemeContext.Provider value={contextValue}>
      <ThemeProvider theme={theme}>
        <CssBaseline />
        {children}
      </ThemeProvider>
    </ThemeContext.Provider>
  );
};
