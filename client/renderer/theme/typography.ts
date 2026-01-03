import { Roboto } from "next/font/google";
// import { TypographyOptions } from '@mui/material/styles';
import { TypographyVariantsOptions } from "@mui/material/styles";

export const roboto = Roboto({
  weight: ["300", "400", "500", "700"],
  subsets: ["latin"],
  display: "swap",
});
export const typographyOptions: TypographyVariantsOptions = {
  fontFamily: roboto.style.fontFamily,
};
