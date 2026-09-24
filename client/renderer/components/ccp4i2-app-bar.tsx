"use client";
import { useEffect, useMemo, useState } from "react";
import { useRouter } from "next/navigation";
import {
  AppBar,
  Box,
  Chip,
  Divider,
  IconButton,
  Menu,
  MenuItem,
  ListItemIcon,
  ListItemText,
  Stack,
  Toolbar,
  Tooltip,
  Typography,
  useMediaQuery,
  useTheme,
} from "@mui/material";
import {
  ArrowBack,
  ArrowForward,
  Home,
  LocalOffer as LocalOfferIcon,
  Logout,
  MoreHoriz,
  Person,
  Science as ScienceIcon,
} from "@mui/icons-material";
import { useMsal } from "@azure/msal-react";
import { CCP4i2MoorhenIcon } from "./General/CCP4i2Icons";
import EditMenu from "./edit-menu";
import FileMenu from "./file-menu";
import HelpMenu from "./help-menu";
import UtilMenu from "./util-menu";
import ViewMenu from "./view-menu";
import { useCCP4i2Window } from "../app-context";
import { useApi } from "../api";
import { Job, Project } from "../types/models";
import { useTopBarState } from "../providers/top-bar-context";
import { useIsJobRoute, useProjectScope } from "../lib/project-scope";
import { useHistoryAvailability } from "../lib/history-availability";
import { isElectron } from "../utils/platform";

interface CampaignInfo {
  campaign_id: number;
  campaign_name: string;
  membership_type: "parent" | "member";
  member_count?: number;
}

const REQUIRE_AUTH = process.env.NEXT_PUBLIC_REQUIRE_AUTH === "true";

// MUI's disabled grey all but disappears against a coloured AppBar, so fade
// the inherited colour instead.
const disabledArrowSx = {
  "&.Mui-disabled": { color: "inherit", opacity: 0.35 },
};

/**
 * The application's one top bar, mounted by the (authed) route-group layout.
 *
 * Everything that does not need a project — history, Home, the File/Edit/View/
 * Utilities/Help menus, the signed-in user — is present on every route. The
 * project-scoped extras (tags, campaign, Moorhen, the project/job label) appear
 * only on routes that actually have a project — see useProjectScope for why
 * that is not simply CCP4i2Context's projectId.
 */
export default function CCP4i2AppBar() {
  const router = useRouter();
  const theme = useTheme();
  const api = useApi();
  const { jobId, setDevMode } = useCCP4i2Window();
  const { title } = useTopBarState();
  const { canGoBack, canGoForward } = useHistoryAvailability();
  const projectId = useProjectScope();
  const showJob = useIsJobRoute();

  const { data: projectData } = api.get<Project>(
    projectId ? `projects/${projectId}` : null
  );
  const { data: jobData } = api.get<Job>(
    showJob && jobId ? `jobs/${jobId}` : null
  );
  // api.get sets keepPreviousData, so a response outlives its key going null.
  // Without re-checking the scope here, the bar went on naming the last
  // project you visited — and showing its tags — on the preferences page.
  const project = projectId ? projectData : null;
  const job = showJob ? jobData : null;
  const { data: campaignInfo } = api.get<Record<string, CampaignInfo>>(
    projectId
      ? `projectgroups/project_campaigns/?project_ids=${projectId}&include_members=true`
      : null
  );
  const projectCampaign =
    projectId && campaignInfo ? campaignInfo[String(projectId)] : null;

  // The project's tags, as objects. Older payloads carried bare ids, which
  // cannot be labelled, so those are simply not shown.
  const projectTags = Array.isArray(project?.tags)
    ? (project.tags as any[]).filter(
        (tag) => typeof tag === "object" && tag !== null
      )
    : [];

  const isXSmall = useMediaQuery(theme.breakpoints.down("sm")); // < 600px
  const isSmall = useMediaQuery(theme.breakpoints.down("md")); // < 900px
  const isMedium = useMediaQuery(theme.breakpoints.down("lg")); // < 1200px

  const [moreAnchorEl, setMoreAnchorEl] = useState<null | HTMLElement>(null);
  const [userMenuAnchor, setUserMenuAnchor] = useState<null | HTMLElement>(null);

  // Home only makes sense in web deployments, and MSAL only runs there too.
  // Both are deferred to an effect so the first client render matches the
  // server-rendered markup.
  const [isWeb, setIsWeb] = useState(false);
  useEffect(() => setIsWeb(!isElectron()), []);

  const showUserMenu = REQUIRE_AUTH && isWeb;
  const msalContext = useMsalSafe();
  const currentUser = showUserMenu ? msalContext?.accounts?.[0] : undefined;

  useEffect(() => {
    if (!window.electronAPI) return;
    window.electronAPI.sendMessage("get-config");
    window.electronAPI.onMessage(
      "message-from-main",
      (event: any, data: any) => {
        if (data.message === "get-config") {
          setDevMode(data.config.devMode);
        }
      }
    );
  }, [setDevMode]);

  const hasCampaign = !!projectCampaign;
  const hasTags = projectTags.length > 0;

  // Which menus stay on the bar, and which fall into the overflow, as the
  // window narrows. The project extras are the first to go because they are
  // shortcuts to things reachable elsewhere.
  const { visible, overflow } = useMemo(() => {
    const extras = [
      ...(hasTags ? ["tags"] : []),
      ...(hasCampaign ? ["campaign"] : []),
    ];
    if (isXSmall)
      return {
        visible: ["file"],
        overflow: ["edit", "view", "util", "help", ...extras],
      };
    if (isSmall)
      return {
        visible: ["file", "edit", "view"],
        overflow: ["util", "help", ...extras],
      };
    if (isMedium)
      return {
        visible: ["file", "edit", "view", "util", "help", ...(hasTags ? ["tags"] : [])],
        overflow: hasCampaign ? ["campaign"] : [],
      };
    return {
      visible: ["file", "edit", "view", "util", "help", ...extras],
      overflow: [],
    };
  }, [isXSmall, isSmall, isMedium, hasTags, hasCampaign]);

  const handleNavigateToCampaign = () => {
    if (projectCampaign) {
      router.push(`/ccp4i2/campaigns/${projectCampaign.campaign_id}`);
    }
  };

  // Open in a new window (Electron's window-open handler turns this into a
  // fresh BrowserWindow) rather than navigating this one — Moorhen has no way
  // back into i2, so replacing the project view would strand the user.
  const handleOpenMoorhen = () => {
    if (projectId) {
      window.open(
        `/ccp4i2/moorhen-page/project/${projectId}`,
        "_blank",
        "noopener,noreferrer"
      );
    }
  };

  const tagChips = projectTags.map((tag) => (
    <Tooltip
      key={tag.id}
      title={
        tag.display_path && tag.display_path !== tag.text
          ? `${tag.display_path} — show projects with this tag`
          : "Show projects with this tag"
      }
    >
      <Chip
        icon={<LocalOfferIcon />}
        label={tag.text}
        size="small"
        onClick={() => router.push(`/ccp4i2?tag=${tag.id}`)}
        sx={{
          bgcolor: "rgba(255, 255, 255, 0.15)",
          color: "inherit",
          "&:hover": { bgcolor: "rgba(255, 255, 255, 0.25)" },
          "& .MuiChip-icon": { color: "inherit" },
        }}
      />
    </Tooltip>
  ));

  return (
    <AppBar position="static">
      <Toolbar variant={isXSmall ? "dense" : "regular"} sx={{ gap: 0.5 }}>
        <Stack direction="row" spacing={0.5} sx={{ mr: 1 }}>
          {/* Greyed out rather than removed when there is nowhere to go, so
              the rest of the bar does not shift sideways on every navigation.
              The span is what lets the tooltip still work over a disabled
              button. */}
          <Tooltip title="Back">
            <span>
              <IconButton
                color="inherit"
                aria-label="Back"
                disabled={!canGoBack}
                onClick={() => router.back()}
                sx={disabledArrowSx}
              >
                <ArrowBack />
              </IconButton>
            </span>
          </Tooltip>
          <Tooltip title="Forward">
            <span>
              <IconButton
                color="inherit"
                aria-label="Forward"
                disabled={!canGoForward}
                onClick={() => router.forward()}
                sx={disabledArrowSx}
              >
                <ArrowForward />
              </IconButton>
            </span>
          </Tooltip>
          {isWeb && (
            <Tooltip title="Home">
              <IconButton color="inherit" aria-label="Home" onClick={() => router.push("/")}>
                <Home />
              </IconButton>
            </Tooltip>
          )}
        </Stack>

        {visible.includes("file") && <FileMenu />}
        {visible.includes("edit") && <EditMenu />}
        {visible.includes("view") && <ViewMenu />}
        {visible.includes("util") && <UtilMenu />}
        {visible.includes("help") && <HelpMenu />}

        {/* A project's classifications, shown while you are working inside
            it. Read-only on purpose: an editor in the menu bar was too much,
            but which tags a project carries is worth seeing. Clicking one
            opens the project list filtered to that tag. */}
        {visible.includes("tags") && (
          <Box sx={{ display: "flex", alignItems: "center", gap: 0.5, ml: 1 }}>
            {tagChips}
          </Box>
        )}

        {visible.includes("campaign") && projectCampaign && (
          <Tooltip title={`Go to ${projectCampaign.campaign_name} campaign`}>
            <Chip
              icon={<ScienceIcon />}
              label={projectCampaign.campaign_name}
              onClick={handleNavigateToCampaign}
              size="small"
              sx={{
                ml: 1,
                bgcolor: "rgba(255, 255, 255, 0.15)",
                color: "inherit",
                "&:hover": { bgcolor: "rgba(255, 255, 255, 0.25)" },
                "& .MuiChip-icon": { color: "inherit" },
              }}
            />
          </Tooltip>
        )}

        {/* Moorhen (project-scoped): opens the 3D viewer with this project's
            context for scene authoring + job/param file resolution. */}
        {projectId && (
          <Tooltip title="Open this project in Moorhen (new window)">
            <IconButton
              color="inherit"
              onClick={handleOpenMoorhen}
              aria-label="Open Moorhen for this project"
              size="small"
              sx={{ ml: 0.5 }}
            >
              <CCP4i2MoorhenIcon />
            </IconButton>
          </Tooltip>
        )}

        {overflow.length > 0 && (
          <>
            <IconButton
              onClick={(event) => setMoreAnchorEl(event.currentTarget)}
              color="inherit"
              aria-label="More menu options"
            >
              <MoreHoriz />
            </IconButton>
            <Menu
              anchorEl={moreAnchorEl}
              open={Boolean(moreAnchorEl)}
              onClose={() => setMoreAnchorEl(null)}
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
              {overflow.includes("tags") &&
                projectTags.map((tag) => (
                  <MenuItem
                    key={tag.id}
                    onClick={() => {
                      router.push(`/ccp4i2?tag=${tag.id}`);
                      setMoreAnchorEl(null);
                    }}
                  >
                    <LocalOfferIcon sx={{ mr: 1 }} fontSize="small" />
                    {tag.display_path ?? tag.text}
                  </MenuItem>
                ))}
              {overflow.includes("campaign") && projectCampaign && (
                <MenuItem
                  onClick={() => {
                    handleNavigateToCampaign();
                    setMoreAnchorEl(null);
                  }}
                >
                  <ScienceIcon sx={{ mr: 1 }} fontSize="small" />
                  {projectCampaign.campaign_name}
                </MenuItem>
              )}
            </Menu>
          </>
        )}

        {/* Where you are: the project (and job) on a project route, otherwise
            whatever the page called itself through useTopBar. */}
        <Typography
          component="h1"
          variant="body1"
          sx={{
            flex: 1,
            overflow: "hidden",
            textAlign: "right",
            textOverflow: "ellipsis",
            whiteSpace: "nowrap",
          }}
        >
          {project ? (
            <>
              {job && `Job ${job.number}: `}
              <Tooltip title={project.description || ""}>
                <span style={{ fontWeight: 500 }}>{project.name}</span>
              </Tooltip>
            </>
          ) : (
            title ?? "CCP4i2"
          )}
        </Typography>

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
                onClick={(event) => setUserMenuAnchor(event.currentTarget)}
                sx={{
                  color: "inherit",
                  borderColor: "rgba(255,255,255,0.5)",
                  cursor: "pointer",
                  "& .MuiChip-icon": { color: "inherit" },
                  "&:hover": { bgcolor: "rgba(255,255,255,0.1)" },
                }}
                variant="outlined"
              />
            </Tooltip>
            <Menu
              anchorEl={userMenuAnchor}
              open={Boolean(userMenuAnchor)}
              onClose={() => setUserMenuAnchor(null)}
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
              <MenuItem
                onClick={() => {
                  setUserMenuAnchor(null);
                  msalContext?.instance?.logoutRedirect();
                }}
              >
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
 * useMsal throws when no MsalProvider is above it, which is the normal case in
 * the desktop app (auth-provider renders its children bare in LocalSession
 * mode).
 */
function useMsalSafe() {
  try {
    return useMsal();
  } catch {
    return null;
  }
}
