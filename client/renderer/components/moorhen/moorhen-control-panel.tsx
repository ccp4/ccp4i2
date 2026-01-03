import { moorhen } from "moorhen/types/moorhen";
import { useSelector } from "react-redux";
import { Box, Stack, Typography } from "@mui/material";
import { CCP4i2HierarchyBrowser } from "./ccp4i2-hierarchy-browser";
import { MoorhenLoadedContent } from "./moorhen-loaded-panel";
import { useTheme } from "../../theme/theme-provider";

interface MoorhenControlPanelProps {
  onFileSelect: (fileId: number) => Promise<void>;
}

export const MoorhenControlPanel: React.FC<MoorhenControlPanelProps> = ({
  onFileSelect,
}) => {
  const { customColors } = useTheme();
  const cootInitialized = useSelector(
    (state: moorhen.State) => state.generalStates.cootInitialized
  );

  if (!cootInitialized) {
    return (
      <Box
        sx={{
          height: "calc(100vh - 75px)",
          overflowY: "auto",
          display: "flex",
          alignItems: "center",
          justifyContent: "center",
          maxWidth: "100%",
        }}
      >
        Loading...
      </Box>
    );
  }

  return (
    <Stack
      direction="column"
      sx={{
        height: "calc(100vh - 100px)",
        width: "100%",
      }}
    >
      {/* Upper section - CCP4i2HierarchyBrowser */}
      <Box
        sx={{
          flex: 1,
          minHeight: 0, // Allows flex item to shrink
          overflow: "hidden",
          display: "flex",
          flexDirection: "column",
        }}
      >
        <CCP4i2HierarchyBrowser onFileSelect={onFileSelect} />
      </Box>

      {/* Lower section - Molecules  */}
      <Box
        sx={{
          flex: 1,
          minHeight: 0, // Allows flex item to shrink
          overflow: "hidden",
          backgroundColor: customColors.ui.lightGray,
          border: `1px solid ${customColors.ui.mediumGray}`,
          borderTop: `2px solid ${customColors.ui.lightBlue}`,
          position: "relative",
          display: "flex",
          flexDirection: "column",
        }}
      >
        <Typography
          variant="caption"
          sx={{
            position: "absolute",
            top: 2,
            left: 6,
            fontSize: "0.7rem",
            fontWeight: 500,
            color: customColors.ui.lightBlue,
            backgroundColor: "rgba(255, 255, 255, 0.8)",
            px: 0.5,
            py: 0.25,
            borderRadius: "2px",
            zIndex: 1,
          }}
        >
          Molecules
        </Typography>
        <Box
          sx={{
            flex: 1,
            overflow: "auto",
            display: "flex",
            alignItems: "center",
            justifyContent: "center",
            padding: 2,
            pt: 3, // Extra top padding to account for label
          }}
        >
          <MoorhenLoadedContent onFileSelect={onFileSelect} type="Molecule" />
        </Box>
      </Box>

      {/* Lower section - Maps  */}
      <Box
        sx={{
          flex: 1,
          minHeight: 0, // Allows flex item to shrink
          overflow: "hidden",
          backgroundColor: customColors.ui.lightGray,
          border: `1px solid ${customColors.ui.mediumGray}`,
          borderTop: `2px solid ${customColors.ui.lightBlue}`,
          position: "relative",
          display: "flex",
          flexDirection: "column",
        }}
      >
        <Typography
          variant="caption"
          sx={{
            position: "absolute",
            top: 2,
            left: 6,
            fontSize: "0.7rem",
            fontWeight: 500,
            color: customColors.ui.lightBlue,
            backgroundColor: "rgba(255, 255, 255, 0.8)",
            px: 0.5,
            py: 0.25,
            borderRadius: "2px",
            zIndex: 1,
          }}
        >
          Maps
        </Typography>
        <Box
          sx={{
            flex: 1,
            overflow: "auto",
            display: "flex",
            alignItems: "center",
            justifyContent: "center",
            padding: 2,
            pt: 3, // Extra top padding to account for label
          }}
        >
          <MoorhenLoadedContent onFileSelect={onFileSelect} type="Map" />
        </Box>
      </Box>
    </Stack>
  );
};
