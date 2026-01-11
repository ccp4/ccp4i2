"use client";
import {
  AppBar,
  Toolbar,
  Typography,
  Button,
  Box,
  IconButton,
  useTheme,
  useMediaQuery,
} from "@mui/material";
import { Home, ArrowBack, Help } from "@mui/icons-material";
import { useRouter, usePathname } from "next/navigation";
import HelpMenu from "./help-menu";

interface CCP4i2TopBarProps {
  title?: string;
  showBackButton?: boolean;
  backPath?: string;
}

/**
 * Simple top bar for CCP4i2 pages that don't have a project context
 * (e.g., project list, new project, import project)
 */
export default function CCP4i2TopBar({
  title = "CCP4i2",
  showBackButton = false,
  backPath,
}: CCP4i2TopBarProps) {
  const router = useRouter();
  const pathname = usePathname();
  const theme = useTheme();
  const isMobile = useMediaQuery(theme.breakpoints.down("sm"));

  // Determine if we're in the root of ccp4i2 (project list)
  const isProjectList = pathname === "/ccp4i2";

  const handleBack = () => {
    if (backPath) {
      router.push(backPath);
    } else {
      router.back();
    }
  };

  const handleHome = () => {
    // Go to app selector if available, otherwise projects list
    router.push("/");
  };

  return (
    <AppBar position="static" sx={{ mb: 2 }}>
      <Toolbar variant={isMobile ? "dense" : "regular"}>
        {/* Home / Back button */}
        {showBackButton ? (
          <IconButton
            edge="start"
            color="inherit"
            onClick={handleBack}
            aria-label="Go back"
            sx={{ mr: 1 }}
          >
            <ArrowBack />
          </IconButton>
        ) : (
          <IconButton
            edge="start"
            color="inherit"
            onClick={handleHome}
            aria-label="Home"
            sx={{ mr: 1 }}
          >
            <Home />
          </IconButton>
        )}

        {/* Title */}
        <Typography
          variant={isMobile ? "h6" : "h5"}
          component="h1"
          sx={{ flexGrow: 1 }}
        >
          {title}
        </Typography>

        {/* Help menu */}
        <HelpMenu />
      </Toolbar>
    </AppBar>
  );
}
