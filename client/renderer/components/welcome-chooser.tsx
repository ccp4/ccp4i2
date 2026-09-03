"use client";
import React from "react";
import { useRouter } from "next/navigation";
import {
  Box,
  Card,
  CardActionArea,
  Container,
  Stack,
  Typography,
  alpha,
} from "@mui/material";
import {
  Add,
  Upload,
  ChevronRight,
} from "@mui/icons-material";
import { useAuth } from "@/lib/compounds/auth-context";

/**
 * First-project welcome. Shown whenever the database has no projects — which is
 * the genuine "nothing here yet" state, whether that's a brand-new install or a
 * migration that failed / a projects directory that was rolled back. It is not a
 * one-shot: if the user ends up back at zero projects, they see it again.
 *
 * Names the ways in rather than assuming the user knows where to look. Uses the
 * same routes as the projects toolbar, so there's one source of truth for each
 * path.
 */

interface Path {
  key: string;
  title: string;
  detail: string;
  icon: React.ReactNode;
  route: string;
  primary?: boolean;
  adminOnly?: boolean;
}

const PATHS: Path[] = [
  {
    key: "new",
    title: "Start a new project",
    detail: "An empty project to set up and run jobs in.",
    icon: <Add />,
    route: "/ccp4i2/new-project",
    primary: true,
  },
  {
    key: "import",
    title: "Import a project",
    detail:
      "From a CCP4i2 project zip, or a project folder already on this " +
      "machine. Imported by copy, so the original is left untouched.",
    icon: <Upload />,
    route: "/ccp4i2/import-project",
  },
  // "Migrate a legacy CCP4i2" is deliberately absent during the alpha: adopting
  // a legacy installation's projects where they lie risks damaging the
  // installation the user still depends on. "Import a project" above is the
  // supported route, and it copies rather than adopts.
];

export const WelcomeChooser: React.FC = () => {
  const router = useRouter();
  const { canAdminister } = useAuth();

  const paths = PATHS.filter((p) => !p.adminOnly || canAdminister);

  return (
    <Container maxWidth="sm" sx={{ py: { xs: 4, sm: 8 } }}>
      <Stack spacing={1} sx={{ mb: 4, textAlign: "center" }}>
        <Typography variant="h4" fontWeight={600}>
          Welcome to CCP4i2
        </Typography>
        <Typography variant="body1" color="text.secondary">
          There are no projects here yet. How would you like to begin?
        </Typography>
      </Stack>

      <Stack spacing={1.5}>
        {paths.map((p) => (
          <Card
            key={p.key}
            variant="outlined"
            sx={{
              borderColor: p.primary ? "primary.main" : "divider",
              bgcolor: p.primary ? "action.hover" : undefined,
            }}
          >
            <CardActionArea onClick={() => router.push(p.route)} sx={{ p: 2 }}>
              <Stack direction="row" spacing={2} alignItems="center">
                <Box
                  sx={{
                    width: 40,
                    height: 40,
                    borderRadius: 1.5,
                    flex: "none",
                    display: "grid",
                    placeItems: "center",
                    color: "primary.main",
                    bgcolor: (theme) => alpha(theme.palette.primary.main, 0.12),
                  }}
                >
                  {p.icon}
                </Box>
                <Box sx={{ flex: 1 }}>
                  <Typography fontWeight={600}>{p.title}</Typography>
                  <Typography variant="body2" color="text.secondary">
                    {p.detail}
                  </Typography>
                </Box>
                <ChevronRight color="action" />
              </Stack>
            </CardActionArea>
          </Card>
        ))}
      </Stack>
    </Container>
  );
};
