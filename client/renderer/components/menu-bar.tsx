import {
  AppBar,
  Typography,
  IconButton,
  Menu,
  MenuItem,
  useMediaQuery,
  useTheme,
  Box,
} from "@mui/material";
import { MoreHoriz } from "@mui/icons-material";
import EditMenu from "./edit-menu";
import FileMenu from "./file-menu";
import HelpMenu from "./help-menu";
import UtilMenu from "./util-menu";
import ViewMenu from "./view-menu";
import { useEffect, useState, useRef } from "react";
import { useCCP4i2Window } from "../app-context";
import { useApi } from "../api";
import { Job, Project } from "../types/models";
import EditableTypography from "./editable-typography";
import HistoryToolbar from "./history-toolbar";
import { useRouter } from "next/navigation";
import { TagsOfProject } from "./tags-of-project";

export default function MenuBar() {
  const { projectId, jobId, devMode, setDevMode } = useCCP4i2Window();
  const api = useApi();
  const { data: project, mutate: mutateProject } = api.get<Project>(
    `projects/${projectId}`
  );
  const { data: job } = api.get<Job>(`jobs/${jobId}`);
  const router = useRouter();
  const theme = useTheme();

  // Responsive breakpoints - we'll use different strategies at different sizes
  const isXSmall = useMediaQuery(theme.breakpoints.down("sm")); // < 600px
  const isSmall = useMediaQuery(theme.breakpoints.down("md")); // < 900px
  const isMedium = useMediaQuery(theme.breakpoints.down("lg")); // < 1200px

  // More menu state
  const [moreAnchorEl, setMoreAnchorEl] = useState<null | HTMLElement>(null);
  const isMoreMenuOpen = Boolean(moreAnchorEl);

  const handleMoreMenuOpen = (event: React.MouseEvent<HTMLElement>) => {
    setMoreAnchorEl(event.currentTarget);
  };

  const handleMoreMenuClose = () => {
    setMoreAnchorEl(null);
  };

  // Determine which menu items to show based on screen size
  const getVisibleMenus = () => {
    if (isXSmall) {
      // Extra small: Only File menu visible, rest in overflow
      return {
        visible: ["file"],
        overflow: ["edit", "view", "util", "help", "tags"],
      };
    } else if (isSmall) {
      // Small: File, Edit, and View visible
      return {
        visible: ["file", "edit", "view"],
        overflow: ["util", "help", "tags"],
      };
    } else if (isMedium) {
      // Medium: All menus visible, but tags might be hidden
      return {
        visible: ["file", "edit", "view", "util", "help"],
        overflow: ["tags"],
      };
    } else {
      // Large: Everything visible
      return {
        visible: ["file", "edit", "view", "util", "help", "tags"],
        overflow: [],
      };
    }
  };

  const { visible, overflow } = getVisibleMenus();

  useEffect(() => {
    // Send a message to the main process to get the config
    if (window.electronAPI) {
      window.electronAPI.sendMessage("get-config");
      // Listen for messages from the main process
      window.electronAPI.onMessage(
        "message-from-main",
        (event: any, data: any) => {
          if (data.message === "get-config") {
            setDevMode(data.config.devMode);
          }
        }
      );
    } else console.log("window.electron is not available");
  }, []);

  return (
    <AppBar position="static">
      <HistoryToolbar>
        {/* Always visible menus based on screen size */}
        {visible.includes("file") && <FileMenu />}
        {visible.includes("edit") && <EditMenu />}
        {visible.includes("view") && <ViewMenu />}
        {visible.includes("util") && <UtilMenu />}
        {visible.includes("help") && <HelpMenu />}
        {visible.includes("tags") && project && (
          <TagsOfProject projectId={project.id} />
        )}

        {/* Overflow menu for hidden items */}
        {overflow.length > 0 && (
          <>
            <IconButton
              onClick={handleMoreMenuOpen}
              color="inherit"
              aria-label="More menu options"
            >
              <MoreHoriz />
            </IconButton>
            <Menu
              anchorEl={moreAnchorEl}
              open={isMoreMenuOpen}
              onClose={handleMoreMenuClose}
              transformOrigin={{ horizontal: "left", vertical: "top" }}
              anchorOrigin={{ horizontal: "left", vertical: "bottom" }}
            >
              {overflow.includes("edit") && (
                <Box sx={{ px: 2, py: 1 }}>
                  <EditMenu />
                </Box>
              )}
              {overflow.includes("view") && (
                <Box sx={{ px: 2, py: 1 }}>
                  <ViewMenu />
                </Box>
              )}
              {overflow.includes("util") && (
                <Box sx={{ px: 2, py: 1 }}>
                  <UtilMenu />
                </Box>
              )}
              {overflow.includes("help") && (
                <Box sx={{ px: 2, py: 1 }}>
                  <HelpMenu />
                </Box>
              )}
              {overflow.includes("tags") && project && (
                <Box sx={{ px: 2, py: 1 }}>
                  <TagsOfProject projectId={project.id} />
                </Box>
              )}
            </Menu>
          </>
        )}

        <Typography sx={{ flexGrow: 1 }} />

        {/* Project/Job info - always visible but might be truncated on mobile */}
        <Box
          sx={{
            display: "flex",
            alignItems: "center",
            minWidth: 0, // Allow shrinking
            maxWidth: isXSmall ? "200px" : "none", // Limit width on mobile
          }}
        >
          {job?.number && (
            <EditableTypography
              variant={isXSmall ? "body1" : "h5"}
              text={`Job ${job.number}: `}
            />
          )}
          {project && (
            <EditableTypography
              variant={isXSmall ? "body1" : "h5"}
              text={project.name}
              onDelay={(name) =>
                api
                  .patch(`projects/${project.id}`, { name: name })
                  .then((_) => {
                    mutateProject();
                  })
              }
            />
          )}
        </Box>
      </HistoryToolbar>
    </AppBar>
  );
}
