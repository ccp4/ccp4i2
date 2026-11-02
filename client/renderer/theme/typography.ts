import localFont from "next/font/local";
// import { TypographyOptions } from '@mui/material/styles';
import { TypographyVariantsOptions } from "@mui/material/styles";

// Roboto ships with the app (fonts/, latin subset, from @fontsource/roboto
// 5.2.6; OFL, see fonts/OFL.txt). next/font/google fetched it at build time,
// and a malformed reply from Google failed release builds at random.
export const roboto = localFont({
  src: [
    { path: "./fonts/roboto-latin-300-normal.woff2", weight: "300", style: "normal" },
    { path: "./fonts/roboto-latin-400-normal.woff2", weight: "400", style: "normal" },
    { path: "./fonts/roboto-latin-500-normal.woff2", weight: "500", style: "normal" },
    { path: "./fonts/roboto-latin-700-normal.woff2", weight: "700", style: "normal" },
  ],
  display: "swap",
});
export const typographyOptions: TypographyVariantsOptions = {
  fontFamily: roboto.style.fontFamily,
};
