/**
 * next/font Roboto loader for the compounds testbed.
 *
 * Mirrors the production shell's `client/renderer/theme/typography.ts` so
 * the standalone compounds dev server (`npm run dev` against
 * apps/compounds/frontend) renders text with the same font-family
 * resolution chain as the deployed CCP4i2 web image.
 *
 * Without this, MUI typography in the testbed inherits a system-font
 * stack and html2canvas captures cleanly — so a card export bug that
 * only manifests under next/font's CSS-variable family
 * (`__Roboto_<hash>__` plus a metric-tuned fallback) is invisible
 * locally and only shows up after build → push → deploy. Aligning the
 * shells closes that loop.
 */

import { Roboto } from 'next/font/google';
import type { TypographyVariantsOptions } from '@mui/material/styles';

export const roboto = Roboto({
  weight: ['300', '400', '500', '700'],
  subsets: ['latin'],
  display: 'swap',
});

export const typographyOptions: TypographyVariantsOptions = {
  fontFamily: roboto.style.fontFamily,
};
