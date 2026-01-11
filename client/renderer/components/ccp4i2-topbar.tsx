"use client";
import { useState } from "react";
import {
  AppBar,
  Toolbar,
  Typography,
  Box,
  IconButton,
  useTheme,
  useMediaQuery,
  Menu,
  MenuItem,
  ListItemIcon,
  ListItemText,
  Divider,
  Chip,
  Tooltip,
} from "@mui/material";
import { Home, ArrowBack, Person, Logout } from "@mui/icons-material";
import { useRouter, usePathname } from "next/navigation";
import { useMsal } from "@azure/msal-react";
import HelpMenu from "./help-menu";
import { isElectron } from "../utils/platform";

interface CCP4i2TopBarProps {
  title?: string;
  showBackButton?: boolean;
  backPath?: string;
}

const REQUIRE_AUTH = process.env.NEXT_PUBLIC_REQUIRE_AUTH === "true";

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
  const [userMenuAnchor, setUserMenuAnchor] = useState<null | HTMLElement>(
    null
  );

  // Only use MSAL when auth is required and not in Electron
  const showUserMenu = REQUIRE_AUTH && !isElectron();
  const msalContext = showUserMenu ? useMsalSafe() : null;
  const accounts = msalContext?.accounts || [];
  const instance = msalContext?.instance;
  const currentUser = accounts[0];

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

  const handleUserMenuOpen = (event: React.MouseEvent<HTMLElement>) => {
    setUserMenuAnchor(event.currentTarget);
  };

  const handleUserMenuClose = () => {
    setUserMenuAnchor(null);
  };

  const handleLogout = () => {
    handleUserMenuClose();
    if (instance) {
      instance.logoutRedirect();
    }
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

        {/* User menu - only show when auth is required and not in Electron */}
        {showUserMenu && currentUser && (
          <>
            <Divider
              orientation="vertical"
              flexItem
              sx={{ mx: 1, borderColor: "rgba(255,255,255,0.3)" }}
            />
            <Tooltip title={currentUser.username || "User"}>
              <Chip
                icon={<Person fontSize="small" sx={{ color: "inherit" }} />}
                label={currentUser.name?.split(" ")[0] || "User"}
                size="small"
                onClick={handleUserMenuOpen}
                sx={{
                  color: "inherit",
                  borderColor: "rgba(255,255,255,0.5)",
                  cursor: "pointer",
                  "& .MuiChip-icon": { color: "inherit" },
                  "&:hover": {
                    bgcolor: "rgba(255,255,255,0.1)",
                  },
                }}
                variant="outlined"
              />
            </Tooltip>
            <Menu
              anchorEl={userMenuAnchor}
              open={Boolean(userMenuAnchor)}
              onClose={handleUserMenuClose}
              anchorOrigin={{ vertical: "bottom", horizontal: "right" }}
              transformOrigin={{ vertical: "top", horizontal: "right" }}
            >
              <MenuItem disabled>
                <ListItemText
                  primary={currentUser.name}
                  secondary={currentUser.username}
                  primaryTypographyProps={{ fontWeight: 500 }}
                />
              </MenuItem>
              <Divider />
              <MenuItem onClick={handleLogout}>
                <ListItemIcon>
                  <Logout fontSize="small" />
                </ListItemIcon>
                <ListItemText>Sign Out</ListItemText>
              </MenuItem>
            </Menu>
          </>
        )}
      </Toolbar>
    </AppBar>
  );
}

/**
 * Safe wrapper for useMsal that handles the case when MsalProvider is not present.
 * This prevents crashes in non-auth environments.
 */
function useMsalSafe() {
  try {
    return useMsal();
  } catch {
    return null;
  }
}
