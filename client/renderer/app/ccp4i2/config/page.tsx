"use client";
import React, { useState } from "react";
import { ConfigContent } from "../../../components/config-content";
import { NavigationShortcutsProvider } from "../../../providers/navigation-shortcuts-provider";
import { TeamsDiagnostic } from "../../../components/teams-diagnostic";
import {
  Box,
  Paper,
  Typography,
  Accordion,
  AccordionSummary,
  AccordionDetails,
  Container,
  Grid2,
  Chip,
} from "@mui/material";
import { ExpandMore, Settings, Security, BugReport } from "@mui/icons-material";

export default function ConfigPage() {
  const [teamsExpanded, setTeamsExpanded] = useState(false);

  const handleTeamsToggle = (
    event: React.SyntheticEvent,
    isExpanded: boolean
  ) => {
    setTeamsExpanded(isExpanded);
  };

  return (
    <NavigationShortcutsProvider>
      <Container maxWidth="xl" sx={{ py: 4 }}>
        <Grid2 container spacing={4}>
          {/* Header Section */}
          <Grid2 size={{ xs: 12 }}>
            <Paper
              elevation={2}
              sx={{
                p: 3,
                background: "linear-gradient(135deg, #667eea 0%, #764ba2 100%)",
                color: "white",
                borderRadius: 2,
              }}
            >
              <Box display="flex" alignItems="center" gap={2}>
                <Settings sx={{ fontSize: 40 }} />
                <Box>
                  <Typography variant="h3" component="h1" fontWeight="bold">
                    Configuration
                  </Typography>
                  <Typography variant="h6" sx={{ opacity: 0.9 }}>
                    System settings and diagnostic tools
                  </Typography>
                </Box>
              </Box>
            </Paper>
          </Grid2>

          {/* Main Configuration Section */}
          <Grid2 size={{ xs: 12, lg: 8 }}>
            <Paper
              elevation={1}
              sx={{
                borderRadius: 2,
                overflow: "hidden",
                border: "1px solid",
                borderColor: "divider",
              }}
            >
              <Box
                sx={{
                  p: 2,
                  bgcolor: "grey.50",
                  borderBottom: "1px solid",
                  borderColor: "divider",
                  display: "flex",
                  alignItems: "center",
                  gap: 1,
                }}
              >
                <Settings color="primary" />
                <Typography variant="h6" fontWeight="600">
                  System Configuration
                </Typography>
              </Box>
              <Box sx={{ p: 3 }}>
                <ConfigContent />
              </Box>
            </Paper>
          </Grid2>

          {/* Info Panel */}
          <Grid2 size={{ xs: 12, lg: 4 }}>
            <Paper
              elevation={1}
              sx={{
                p: 3,
                borderRadius: 2,
                border: "1px solid",
                borderColor: "divider",
                height: "fit-content",
              }}
            >
              <Box display="flex" alignItems="center" gap={1} mb={2}>
                <Security color="primary" />
                <Typography variant="h6" fontWeight="600">
                  Configuration Status
                </Typography>
              </Box>

              <Box sx={{ display: "flex", flexDirection: "column", gap: 2 }}>
                <Box>
                  <Typography variant="subtitle2" color="text.secondary" mb={1}>
                    Environment
                  </Typography>
                  <Chip
                    label={
                      process.env.NODE_ENV === "development"
                        ? "Development"
                        : "Production"
                    }
                    color={
                      process.env.NODE_ENV === "development"
                        ? "warning"
                        : "success"
                    }
                    variant="outlined"
                    size="small"
                  />
                </Box>

                <Box>
                  <Typography variant="subtitle2" color="text.secondary" mb={1}>
                    Authentication
                  </Typography>
                  <Chip
                    label={
                      process.env.NEXT_PUBLIC_REQUIRE_AUTH === "true"
                        ? "Enabled"
                        : "Disabled"
                    }
                    color={
                      process.env.NEXT_PUBLIC_REQUIRE_AUTH === "true"
                        ? "success"
                        : "default"
                    }
                    variant="outlined"
                    size="small"
                  />
                </Box>

                <Box>
                  <Typography
                    variant="body2"
                    color="text.secondary"
                    sx={{ mt: 2 }}
                  >
                    Use the diagnostic tools below to troubleshoot
                    authentication and Teams integration issues.
                  </Typography>
                </Box>
              </Box>
            </Paper>
          </Grid2>

          {/* Teams Diagnostic Section - Collapsible */}
          <Grid2 size={{ xs: 12 }}>
            <Accordion
              expanded={teamsExpanded}
              onChange={handleTeamsToggle}
              sx={{
                border: "1px solid",
                borderColor: "divider",
                borderRadius: 2,
                "&:before": {
                  display: "none",
                },
                "&.Mui-expanded": {
                  margin: 0,
                },
              }}
              elevation={1}
            >
              <AccordionSummary
                expandIcon={<ExpandMore />}
                sx={{
                  bgcolor: "grey.50",
                  borderBottom: teamsExpanded ? "1px solid" : "none",
                  borderColor: "divider",
                  minHeight: 64,
                  "& .MuiAccordionSummary-content": {
                    alignItems: "center",
                  },
                  "& .MuiAccordionSummary-content.Mui-expanded": {
                    margin: "12px 0",
                  },
                }}
              >
                <Box display="flex" alignItems="center" gap={2}>
                  <BugReport color="primary" />
                  <Box>
                    <Typography variant="h6" fontWeight="600">
                      Teams Integration Diagnostic
                    </Typography>
                    <Typography variant="body2" color="text.secondary">
                      Test and troubleshoot Microsoft Teams authentication
                    </Typography>
                  </Box>
                  <Chip
                    label="Developer Tools"
                    color="info"
                    variant="outlined"
                    size="small"
                    sx={{ ml: "auto", mr: 2 }}
                  />
                </Box>
              </AccordionSummary>
              <AccordionDetails sx={{ p: 0 }}>
                <Box sx={{ p: 3, pt: 2 }}>
                  <TeamsDiagnostic />
                </Box>
              </AccordionDetails>
            </Accordion>
          </Grid2>
        </Grid2>
      </Container>
    </NavigationShortcutsProvider>
  );
}
