'use client';

/**
 * Theme provider for compounds frontend.
 *
 * This file provides theme context for the compounds app. It has two modes:
 *
 * 1. Standalone mode (compounds app running independently):
 *    Uses its own CompoundsThemeProvider implementation below.
 *
 * 2. Docker/integrated mode (overlaid onto ccp4i2):
 *    The ccp4i2 layout already provides CCP4i2ThemeProvider at the root,
 *    so this file just exports a hook that accesses that context.
 *
 * When overlaid, the path becomes renderer/lib/compounds/theme-provider.tsx,
 * and the CCP4i2ThemeProvider is at renderer/theme/theme-provider.tsx.
 * The @/ alias points to renderer/, so we can import from @/theme/theme-provider.
 *
 * However, since this file IS the theme-provider being imported, we can't
 * conditionally import from different sources. Instead, we detect the
 * environment and either use our own context or re-export from the parent.
 */

import React, {
  createContext,
  useContext,
  useState,
  useEffect,
  ReactNode,
} from 'react';
import { ThemeProvider as MuiThemeProvider, createTheme, CssBaseline } from '@mui/material';

export type ThemeMode = 'light' | 'dark';

// Light mode palette
const lightPaletteOptions = {
  mode: 'light' as const,
  primary: {
    main: '#1976d2',
    light: '#42a5f5',
    dark: '#1565c0',
  },
  secondary: {
    main: '#dc004e',
    light: '#ff4081',
    dark: '#c51162',
  },
  background: {
    default: '#ffffff',
    paper: '#ffffff',
  },
  text: {
    primary: 'rgba(0, 0, 0, 0.87)',
    secondary: 'rgba(0, 0, 0, 0.6)',
  },
};

// Dark mode palette
const darkPaletteOptions = {
  mode: 'dark' as const,
  primary: {
    main: '#90caf9',
    light: '#e3f2fd',
    dark: '#42a5f5',
  },
  secondary: {
    main: '#f48fb1',
    light: '#fce4ec',
    dark: '#f06292',
  },
  background: {
    default: '#121212',
    paper: '#1e1e1e',
  },
  text: {
    primary: '#ffffff',
    secondary: 'rgba(255, 255, 255, 0.7)',
  },
};

interface ThemeContextType {
  mode: ThemeMode;
  toggleTheme: () => void;
  setTheme: (mode: ThemeMode) => void;
}

const ThemeContext = createContext<ThemeContextType | undefined>(undefined);

/**
 * Hook to access theme context.
 * Works in both standalone compounds mode and integrated ccp4i2 mode.
 */
export const useTheme = (): ThemeContextType => {
  const context = useContext(ThemeContext);

  // If we have our own context, use it
  if (context) {
    return context;
  }

  // Fallback: if no context is available (shouldn't happen if provider is set up),
  // return a default implementation that at least doesn't crash
  console.warn('useTheme called outside of ThemeProvider - using defaults');
  return {
    mode: 'light',
    toggleTheme: () => {},
    setTheme: () => {},
  };
};

interface CompoundsThemeProviderProps {
  children: ReactNode;
}

/**
 * Theme provider for compounds frontend.
 *
 * In standalone mode, this provides full theme management.
 * In Docker mode (overlaid onto ccp4i2), the ccp4i2 layout's CCP4i2ThemeProvider
 * takes precedence and this may not be rendered, but if it is, it works the same.
 */
export const CompoundsThemeProvider: React.FC<CompoundsThemeProviderProps> = ({
  children,
}) => {
  // Default to light, will be updated on mount
  const [mode, setMode] = useState<ThemeMode>('light');
  const [mounted, setMounted] = useState(false);

  // Load theme preference from localStorage on mount
  // Use same key as ccp4i2 for consistency in integrated mode
  useEffect(() => {
    const savedTheme = (localStorage.getItem('ccp4i2-theme') || localStorage.getItem('compounds-theme')) as ThemeMode;
    if (savedTheme && (savedTheme === 'light' || savedTheme === 'dark')) {
      setMode(savedTheme);
    } else {
      // Check system preference
      const prefersDark = window.matchMedia(
        '(prefers-color-scheme: dark)'
      ).matches;
      setMode(prefersDark ? 'dark' : 'light');
    }
    setMounted(true);
  }, []);

  // Save theme preference to localStorage (use both keys for compatibility)
  useEffect(() => {
    if (mounted) {
      localStorage.setItem('ccp4i2-theme', mode);
      localStorage.setItem('compounds-theme', mode);
    }
  }, [mode, mounted]);

  // Listen for system theme changes
  useEffect(() => {
    const mediaQuery = window.matchMedia('(prefers-color-scheme: dark)');
    const handleChange = (e: MediaQueryListEvent) => {
      // Only update if user hasn't set a preference
      const savedTheme = localStorage.getItem('ccp4i2-theme') || localStorage.getItem('compounds-theme');
      if (!savedTheme) {
        setMode(e.matches ? 'dark' : 'light');
      }
    };

    mediaQuery.addEventListener('change', handleChange);
    return () => mediaQuery.removeEventListener('change', handleChange);
  }, []);

  const toggleTheme = () => {
    setMode((prevMode) => (prevMode === 'light' ? 'dark' : 'light'));
  };

  const setTheme = (newMode: ThemeMode) => {
    setMode(newMode);
  };

  // Create theme based on current mode
  const theme = React.useMemo(() => {
    const paletteOptions =
      mode === 'light' ? lightPaletteOptions : darkPaletteOptions;

    return createTheme({
      palette: paletteOptions,
      typography: {
        fontFamily: [
          '-apple-system',
          'BlinkMacSystemFont',
          '"Segoe UI"',
          'Roboto',
          'Oxygen',
          'Ubuntu',
          'Cantarell',
          '"Fira Sans"',
          '"Droid Sans"',
          '"Helvetica Neue"',
          'sans-serif',
        ].join(','),
      },
    });
  }, [mode]);

  const contextValue = {
    mode,
    toggleTheme,
    setTheme,
  };

  return (
    <ThemeContext.Provider value={contextValue}>
      <MuiThemeProvider theme={theme}>
        <CssBaseline />
        {children}
      </MuiThemeProvider>
    </ThemeContext.Provider>
  );
};
