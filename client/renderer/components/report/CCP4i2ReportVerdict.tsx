import React, { JSX, useCallback, useEffect, useState } from "react";
import { apiText } from "../../api-fetch";
import {
  Box,
  Typography,
  Chip,
  Divider,
  List,
  ListItem,
  ListItemText,
  Paper,
  Grid,
  Skeleton,
  Alert,
  Accordion,
  AccordionSummary,
  AccordionDetails,
  Collapse,
  IconButton,
  Button,
  CircularProgress,
} from "@mui/material";
import {
  CheckCircle,
  Warning,
  Error,
  Info,
  ExpandMore,
  KeyboardArrowDown,
  KeyboardArrowRight,
  ContentCopy,
} from "@mui/icons-material";
import { CCP4i2ReportElementProps } from "./CCP4i2ReportElement";
import { useApi } from "../../api";
import { useProject } from "../../utils";
import { Job } from "../../types/models";
import { usePopcorn } from "../../providers/popcorn-provider";
import { useRouter } from "next/navigation";
import { useTheme } from "../../theme/theme-provider";

interface SuggestedParameters {
  [key: string]: string;
}

interface AdditionalNode {
  tag: string;
  value: number;
  displayName: string;
}

interface VerdictData {
  verdictScore: number;
  verdictMessage: string;
  bottomline: string;
  suggestedParameters: SuggestedParameters;
  additionalNodes: AdditionalNode[];
}

// SVG Meter Component
const VerdictMeter: React.FC<{ score: number; size?: number; isDark?: boolean }> = ({
  score,
  size = 120,
  isDark = false,
}) => {
  const radius = size / 2 - 10;
  const circumference = Math.PI * radius; // Half circle
  const strokeWidth = 8;
  const center = size / 2;

  // Calculate the stroke offset based on score (0-100)
  const offset = circumference - (score / 100) * circumference;

  // Determine color based on score ranges
  const getScoreColor = (score: number): string => {
    if (score >= 80) return "#4caf50"; // Green
    if (score >= 60) return "#2196f3"; // Blue
    if (score >= 40) return "#ff9800"; // Orange
    return "#f44336"; // Red
  };

  const scoreColor = getScoreColor(score);

  // Theme-aware colors
  const tickColor = isDark ? "#999" : "#666";
  const needleColor = isDark ? "#ccc" : "#333";
  const bgArcColor = isDark ? "#444" : "#e0e0e0";

  // Create tick marks for the scale
  const createTickMarks = () => {
    const ticks: JSX.Element[] = [];
    for (let i = 0; i <= 10; i++) {
      const angle = -Math.PI + (i * Math.PI) / 10;
      const x1 = center + (radius - 5) * Math.cos(angle);
      const y1 = center + (radius - 5) * Math.sin(angle);
      const x2 = center + (radius + 2) * Math.cos(angle);
      const y2 = center + (radius + 2) * Math.sin(angle);

      const isMainTick = i % 2 === 0;

      ticks.push(
        <line
          key={i}
          x1={x1}
          y1={y1}
          x2={x2}
          y2={y2}
          stroke={tickColor}
          strokeWidth={isMainTick ? 2 : 1}
        />
      );

      // Add labels for main ticks
      if (isMainTick) {
        const labelX = center + (radius + 12) * Math.cos(angle);
        const labelY = center + (radius + 12) * Math.sin(angle);
        ticks.push(
          <text
            key={`label-${i}`}
            x={labelX}
            y={labelY}
            textAnchor="middle"
            dominantBaseline="middle"
            fontSize="10"
            fill={tickColor}
          >
            {i * 10}
          </text>
        );
      }
    }
    return ticks;
  };

  // Calculate needle position
  const needleAngle = -Math.PI + (score / 100) * Math.PI;
  const needleLength = radius - 15;
  const needleX = center + needleLength * Math.cos(needleAngle);
  const needleY = center + needleLength * Math.sin(needleAngle);

  return (
    <Box display="flex" flexDirection="column" alignItems="center">
      <svg width={size} height={size * 0.7} style={{ overflow: "visible" }}>
        {/* Background arc */}
        <path
          d={`M ${center - radius} ${center} A ${radius} ${radius} 0 0 1 ${
            center + radius
          } ${center}`}
          fill="none"
          stroke={bgArcColor}
          strokeWidth={strokeWidth}
          strokeLinecap="round"
        />

        {/* Score arc */}
        <path
          d={`M ${center - radius} ${center} A ${radius} ${radius} 0 0 1 ${
            center + radius
          } ${center}`}
          fill="none"
          stroke={scoreColor}
          strokeWidth={strokeWidth}
          strokeLinecap="round"
          strokeDasharray={circumference}
          strokeDashoffset={offset}
          style={{
            transition: "stroke-dashoffset 1s ease-in-out",
          }}
        />

        {/* Tick marks and labels */}
        {createTickMarks()}

        {/* Needle */}
        <line
          x1={center}
          y1={center}
          x2={needleX}
          y2={needleY}
          stroke={needleColor}
          strokeWidth={3}
          strokeLinecap="round"
        />

        {/* Center dot */}
        <circle cx={center} cy={center} r={4} fill={needleColor} />

        {/* Score text */}
        <text
          x={center}
          y={center + 25}
          textAnchor="middle"
          dominantBaseline="middle"
          fontSize="16"
          fontWeight="bold"
          fill={scoreColor}
        >
          {score.toFixed(1)}
        </text>
      </svg>
    </Box>
  );
};

export const CCP4i2ReportVerdict: React.FC<CCP4i2ReportElementProps> = ({
  job,
}) => {
  const [verdictData, setVerdictData] = useState<VerdictData | null>(null);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);
  const [expanded, setExpanded] = useState(true);
  const [parametersExpanded, setParametersExpanded] = useState(false);
  const [creatingTask, setCreatingTask] = useState(false);
  const api = useApi();
  const { project, mutateJobs } = useProject(job.project);
  const { setMessage } = usePopcorn();
  const router = useRouter();
  const { mode } = useTheme();
  const isDark = mode === "dark";

  // Define the core nodes that should not be treated as statistics
  const coreNodes = new Set([
    "verdict_score",
    "verdict_message",
    "bottomline",
    "suggestedParameters",
  ]);

  //Function to create task with suggested parameters
  const createPatchedTask = useCallback(
    async (updatedParams: SuggestedParameters) => {
      try {
        setCreatingTask(true);
        const cloneResult: any = await api.post(`jobs/${job.id}/clone/`);
        if (cloneResult?.success === false) {
          setMessage(`Failed to clone job: ${cloneResult?.error || "Unknown error"}`, "error");
          return;
        }
        if (!cloneResult?.id) {
          setMessage("Failed to clone job: No job ID returned", "error");
          return;
        }
        await mutateJobs();

        // Apply all suggested parameters in series
        for (const [key, value] of Object.entries(updatedParams)) {
          if (key) {
            try {
              // Await each call to ensure serial execution
              const paramResult: any = await api.post(`jobs/${cloneResult.id}/set_parameter/`, {
                object_path: key,
                value,
              });
              if (paramResult?.success === false) {
                console.warn(`Failed to set parameter ${key}:`, paramResult?.error);
              }
            } catch (paramError) {
              console.warn(`Failed to set parameter ${key}:`, paramError);
            }
          }
        }
        const runResult: any = await api.post(`jobs/${cloneResult.id}/run/`);
        if (runResult?.success === false) {
          setMessage(`Failed to run job: ${runResult?.error || "Unknown error"}`, "error");
          return;
        }
        setMessage(
          `Submitted job ${runResult?.number}: ${runResult?.task_name}`,
          "success"
        );
        if (runResult?.id) {
          mutateJobs();
          router.push(`/project/${job.project}/job/${runResult.id}`);
        }

        console.log(
          `Created new job ${cloneResult.id} with suggested parameters`
        );
      } catch (error) {
        console.error("Error creating patched task:", error);
        setMessage(`Error creating task: ${error instanceof Error ? error.message : String(error)}`, "error");
      } finally {
        setCreatingTask(false);
      }
    },
    [job.id, api, mutateJobs, setMessage, router, job.project]
  );

  const handleApplySuggestedParameters = useCallback(async () => {
    if (verdictData?.suggestedParameters) {
      await createPatchedTask(verdictData.suggestedParameters);
    }
  }, [verdictData?.suggestedParameters, createPatchedTask]);

  useEffect(() => {
    const fetchVerdictData = async () => {
      if (!project) return;
      if (!job?.id) {
        setError("No job ID provided");
        setLoading(false);
        return;
      }

      try {
        setLoading(true);
        setError(null);

        // Fetch the program.xml file for this job
        const composite_path = `/api/proxy/projects/${
          project.id
        }/project_file?path=${encodeURIComponent(
          `CCP4_JOBS/job_${job.number}/program.xml`
        )}`;

        const xmlText = await apiText(composite_path, {
          headers: {
            "Content-Type": "application/xml",
          },
        });

        // Parse the XML
        const parser = new DOMParser();
        const xmlDoc = parser.parseFromString(xmlText, "text/xml");

        // Check for parsing errors
        const parserError = xmlDoc.querySelector("parsererror");
        if (parserError) {
          console.log("Failed to parse program.xml");
          return;
        }

        // Find the Verdict node
        const verdictNode = xmlDoc.querySelector("Verdict");
        if (!verdictNode) {
          console.log("No Verdict node found in program.xml");
          return;
        }

        // Extract verdict data
        const getTextContent = (selector: string): string => {
          const element = verdictNode.querySelector(selector);
          return element?.textContent || "";
        };

        const getNumberContent = (selector: string): number | undefined => {
          const text = getTextContent(selector);
          const num = parseFloat(text);
          return isNaN(num) ? undefined : num;
        };

        // Extract suggested parameters
        const suggestedParameters: SuggestedParameters = {};
        const suggestedParamsNode = verdictNode.querySelector(
          "suggestedParameters"
        );
        if (suggestedParamsNode) {
          Array.from(suggestedParamsNode.children).forEach((child) => {
            if (child.tagName && child.textContent) {
              suggestedParameters[child.tagName] = child.textContent;
            }
          });
        }

        // Extract additional nodes (anything that's not in coreNodes and is numeric)
        const additionalNodes: AdditionalNode[] = [];
        Array.from(verdictNode.children).forEach((child) => {
          const tagName = child.tagName;
          const textContent = child.textContent;

          if (!coreNodes.has(tagName) && textContent) {
            const numValue = parseFloat(textContent);
            if (!isNaN(numValue)) {
              additionalNodes.push({
                tag: tagName,
                value: numValue,
                displayName: formatStatisticName(tagName),
              });
            }
          }
        });

        const data: VerdictData = {
          verdictScore: getNumberContent("verdict_score") || 0,
          verdictMessage: getTextContent("verdict_message"),
          bottomline: getTextContent("bottomline"),
          suggestedParameters,
          additionalNodes,
        };

        setVerdictData(data);
      } catch (err) {
        console.error("Error fetching verdict data:", err);
        if (
          err &&
          typeof err === "object" &&
          "message" in err &&
          typeof (err as any).message === "string"
        ) {
          setError((err as any).message);
        } else {
          setError("Unknown error occurred");
        }
      } finally {
        setLoading(false);
      }
    };

    fetchVerdictData();
  }, [job?.id, project]);

  // Format statistic names for display
  const formatStatisticName = (tagName: string): string => {
    // Handle specific known cases
    const nameMap: { [key: string]: string } = {
      meanRfree: "Mean R-free",
      medianClash: "Median Clash Score",
      ramaOutliers: "Ramachandran Outliers (%)",
      rWork: "R-work",
      rFree: "R-free",
      bondRMSD: "Bond RMSD",
      angleRMSD: "Angle RMSD",
      clashScore: "Clash Score",
      molProbityScore: "MolProbity Score",
    };

    if (nameMap[tagName]) {
      return nameMap[tagName];
    }

    // Default formatting: convert camelCase to title case
    return tagName
      .replace(/([A-Z])/g, " $1")
      .replace(/^./, (str) => str.toUpperCase())
      .trim();
  };

  // Format values based on the statistic type
  const formatStatisticValue = (tag: string, value: number): string => {
    // Percentages
    if (
      tag.toLowerCase().includes("outlier") ||
      tag.toLowerCase().includes("percent")
    ) {
      return `${value.toFixed(1)}%`;
    }

    // R-factors and RMSD values (typically 3 decimal places)
    if (
      tag.toLowerCase().includes("rfree") ||
      tag.toLowerCase().includes("rwork") ||
      tag.toLowerCase().includes("rmsd")
    ) {
      return value.toFixed(3);
    }

    // Clash scores (1 decimal place)
    if (tag.toLowerCase().includes("clash")) {
      return value.toFixed(1);
    }

    // Default: 2 decimal places
    return value.toFixed(2);
  };

  // Determine verdict quality based on score
  const getVerdictQuality = (score: number) => {
    if (score >= 80)
      return {
        label: "Excellent",
        color: "success" as const,
        icon: <CheckCircle />,
      };
    if (score >= 60)
      return { label: "Good", color: "primary" as const, icon: <Info /> };
    if (score >= 40)
      return { label: "Fair", color: "warning" as const, icon: <Warning /> };
    return { label: "Poor", color: "error" as const, icon: <Error /> };
  };

  // Decode HTML entities in messages
  const decodeHtml = (html: string): string => {
    const textarea = document.createElement("textarea");
    textarea.innerHTML = html;
    return textarea.value;
  };

  // Clean and format HTML content for display
  const formatMessage = (htmlString: string): string => {
    if (!htmlString) return "";
    const decoded = decodeHtml(htmlString);
    // Remove HTML tags for clean text display
    return decoded.replace(/<[^>]*>/g, "").replace(/&nbsp;/g, " ");
  };

  // Format parameter names for display
  const formatParameterName = (paramName: string): string => {
    return paramName
      .replace(/_/g, " ")
      .toLowerCase()
      .replace(/\b\w/g, (l) => l.toUpperCase());
  };

  const handleAccordionChange = (
    _event: React.SyntheticEvent,
    isExpanded: boolean
  ) => {
    setExpanded(isExpanded);
  };

  const handleParametersToggle = () => {
    setParametersExpanded(!parametersExpanded);
  };

  // Render accordion header content
  const renderAccordionHeader = () => {
    if (loading) {
      return (
        <Box display="flex" alignItems="center">
          <Skeleton variant="circular" width={24} height={24} sx={{ mr: 1 }} />
          <Skeleton variant="text" width={200} />
        </Box>
      );
    }

    if (error || !verdictData) {
      return (
        <Box display="flex" alignItems="center">
          <Error sx={{ mr: 1, color: "error.main" }} />
          <Typography variant="subtitle1">
            Structure Quality Verdict - Error
          </Typography>
        </Box>
      );
    }

    const quality = getVerdictQuality(verdictData.verdictScore);

    return (
      <Box display="flex" alignItems="center" width="100%">
        {quality.icon}
        <Box ml={1} flexGrow={1}>
          <Typography variant="subtitle1">Structure Quality Verdict</Typography>
        </Box>
        <Chip
          label={`${quality.label} (${verdictData.verdictScore.toFixed(1)}%)`}
          color={quality.color}
          variant="filled"
          size="small"
          sx={{ mr: 2 }}
        />
      </Box>
    );
  };

  // Render accordion content
  const renderAccordionContent = () => {
    // Loading state
    if (loading) {
      return (
        <Box>
          <Skeleton variant="rectangular" height={300} />
          <Box mt={2}>
            <Skeleton variant="text" height={40} />
            <Skeleton variant="text" height={20} />
            <Skeleton variant="text" height={20} />
          </Box>
        </Box>
      );
    }

    // Error state
    if (error) {
      return (
        <Alert severity="error">
          <Typography variant="h6">Failed to load verdict data</Typography>
          <Typography variant="body2">{error}</Typography>
        </Alert>
      );
    }

    // No data state
    if (!verdictData) {
      return (
        <Alert severity="info">
          <Typography variant="body1">No verdict data available</Typography>
        </Alert>
      );
    }

    return (
      <>
        {/* Score Meter and Assessment */}
        <Grid container spacing={3} mb={3}>
          <Grid item xs={12} md={4}>
            <Paper elevation={1} sx={{ p: 2, textAlign: "center" }}>
              <Typography variant="h6" gutterBottom>
                Quality Score
              </Typography>
              <VerdictMeter score={verdictData.verdictScore} size={140} isDark={isDark} />
            </Paper>
          </Grid>
          <Grid item xs={12} md={8}>
            {verdictData.verdictMessage && (
              <Paper
                elevation={1}
                sx={{ p: 2, backgroundColor: isDark ? "grey.900" : "grey.50", height: "100%" }}
              >
                <Typography variant="body1" gutterBottom>
                  <strong>Assessment:</strong>
                </Typography>
                <Typography variant="body2">
                  {formatMessage(verdictData.verdictMessage)}
                </Typography>
              </Paper>
            )}
          </Grid>
        </Grid>

        {/* Statistics - render all additional nodes dynamically */}
        {verdictData.additionalNodes.length > 0 && (
          <>
            <Typography variant="h6" gutterBottom>
              Structure Statistics
            </Typography>
            <Grid container spacing={2} mb={3}>
              {verdictData.additionalNodes.map((node) => (
                <Grid item xs={12} sm={4} md={3} key={node.tag}>
                  <Paper elevation={1} sx={{ p: 2, textAlign: "center" }}>
                    <Typography variant="body2" color="textSecondary">
                      {node.displayName}
                    </Typography>
                    <Typography variant="h6">
                      {formatStatisticValue(node.tag, node.value)}
                    </Typography>
                  </Paper>
                </Grid>
              ))}
            </Grid>
          </>
        )}

        {/* Detailed recommendations */}
        {verdictData.bottomline && (
          <>
            <Typography variant="h6" gutterBottom>
              Detailed Recommendations
            </Typography>
            <Paper elevation={1} sx={{ p: 2, mb: 3 }}>
              <Typography variant="body2" sx={{ whiteSpace: "pre-line" }}>
                {formatMessage(verdictData.bottomline)}
              </Typography>
            </Paper>
          </>
        )}

        {/* Suggested parameters - now collapsible with action button */}
        {Object.keys(verdictData.suggestedParameters).length > 0 && (
          <>
            <Divider sx={{ my: 2 }} />
            <Box
              display="flex"
              alignItems="center"
              justifyContent="space-between"
              mb={1}
            >
              <Box
                display="flex"
                alignItems="center"
                sx={{ cursor: "pointer" }}
                onClick={handleParametersToggle}
              >
                <IconButton size="small">
                  {parametersExpanded ? (
                    <KeyboardArrowDown />
                  ) : (
                    <KeyboardArrowRight />
                  )}
                </IconButton>
                <Typography variant="h6">
                  Suggested Parameter Changes
                </Typography>
              </Box>

              <Button
                variant="contained"
                color="primary"
                size="small"
                startIcon={
                  creatingTask ? (
                    <CircularProgress size={16} color="inherit" />
                  ) : (
                    <ContentCopy />
                  )
                }
                onClick={handleApplySuggestedParameters}
                disabled={
                  creatingTask ||
                  Object.keys(verdictData.suggestedParameters).length === 0
                }
                sx={{ ml: 2 }}
              >
                {creatingTask ? "Creating..." : "Apply & Run"}
              </Button>
            </Box>

            <Collapse in={parametersExpanded} timeout="auto" unmountOnExit>
              <List dense sx={{ mt: 1 }}>
                {Object.entries(verdictData.suggestedParameters).map(
                  ([key, value]) => (
                    <ListItem key={key} sx={{ py: 0.5 }}>
                      <ListItemText
                        primary={
                          <Box
                            display="flex"
                            justifyContent="space-between"
                            alignItems="center"
                          >
                            <Typography variant="body2" fontWeight="medium">
                              {formatParameterName(key)}
                            </Typography>
                            <Chip
                              label={value}
                              size="small"
                              variant="outlined"
                              color="primary"
                            />
                          </Box>
                        }
                      />
                    </ListItem>
                  )
                )}
              </List>
            </Collapse>
          </>
        )}
      </>
    );
  };

  return (
    <Accordion
      expanded={expanded}
      onChange={handleAccordionChange}
      disableGutters
      elevation={1}
    >
      <AccordionSummary
        expandIcon={<ExpandMore />}
        aria-controls="verdict-content"
        id="verdict-header"
        sx={{
          backgroundColor: isDark ? "grey.800" : "grey.100",
          "&:hover": {
            backgroundColor: isDark ? "grey.700" : "grey.200",
          },
          minHeight: 48,
          "& .MuiAccordionSummary-content": {
            margin: "8px 0",
          },
          "& .MuiAccordionSummary-content.Mui-expanded": {
            margin: "8px 0",
          },
        }}
      >
        {renderAccordionHeader()}
      </AccordionSummary>

      <AccordionDetails sx={{ p: 2 }}>
        {renderAccordionContent()}
      </AccordionDetails>
    </Accordion>
  );
};
