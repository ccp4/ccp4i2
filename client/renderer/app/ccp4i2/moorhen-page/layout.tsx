"use client";
import { PropsWithChildren } from "react";
import { CootProvider } from "../../../providers/coot-provider";
import {
  AppBar,
  Breadcrumbs,
  Link as MuiLink,
  Toolbar,
  Typography,
} from "@mui/material";
import NavigateNextIcon from "@mui/icons-material/NavigateNext";
import { useRouter } from "next/navigation";
import { NavigationShortcutsProvider } from "../../../providers/navigation-shortcuts-provider";
import { ThemeToggle } from "../../../components/theme-toggle";
import {
  MoorhenBreadcrumbProvider,
  useMoorhenBreadcrumbs,
} from "../../../providers/moorhen-breadcrumb-context";

function MoorhenAppBar() {
  const router = useRouter();
  const { breadcrumbs } = useMoorhenBreadcrumbs();

  return (
    <AppBar position="static">
      <Toolbar sx={{ gap: 2 }}>
        <ThemeToggle />
        {breadcrumbs.length > 0 ? (
          <Breadcrumbs
            separator={<NavigateNextIcon fontSize="small" sx={{ color: "rgba(255,255,255,0.7)" }} />}
            sx={{ color: "white", margin: "10px" }}
          >
            {breadcrumbs.map((crumb, index) =>
              index < breadcrumbs.length - 1 ? (
                <MuiLink
                  key={index}
                  color="inherit"
                  href={crumb.href}
                  onClick={(e) => {
                    e.preventDefault();
                    router.push(crumb.href);
                  }}
                  underline="hover"
                  sx={{ cursor: "pointer" }}
                >
                  {crumb.label}
                </MuiLink>
              ) : (
                <Typography key={index} variant="h6" color="inherit">
                  {crumb.label}
                </Typography>
              )
            )}
          </Breadcrumbs>
        ) : (
          <Typography
            variant="h6"
            style={{ textAlign: "center", margin: "10px" }}
          >
            Moorhen Viewer
          </Typography>
        )}
      </Toolbar>
    </AppBar>
  );
}

export default function MoorhenPageLayout(props: PropsWithChildren) {
  return (
    <CootProvider>
      <NavigationShortcutsProvider>
        <MoorhenBreadcrumbProvider>
          <MoorhenAppBar />
          {props.children}
        </MoorhenBreadcrumbProvider>
      </NavigationShortcutsProvider>
    </CootProvider>
  );
}
