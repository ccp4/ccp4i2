/**
 * Re-export theme provider from the core theme directory.
 *
 * This file ensures compounds components that import from '@/lib/compounds/theme-provider'
 * use the same React context as the root layout's CCP4i2ThemeProvider.
 *
 * The Docker build copies compounds components but NOT the compounds theme-provider.tsx,
 * so this file is used instead, maintaining a single shared theme context.
 */
export { useTheme, CCP4i2ThemeProvider as CompoundsThemeProvider } from '../../theme/theme-provider';
export type { ThemeMode } from '../../theme/theme-provider';
