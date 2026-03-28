/*
 * Copyright (C) 2026 Newcastle University
 *
 * This file is part of CCP4i2.
 *
 * CCP4i2 is free software: you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License version 3,
 * modified in accordance with the provisions of the license to address
 * the requirements of UK law.
 *
 * See https://www.ccp4.ac.uk/ccp4license.php for details.
 */
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
