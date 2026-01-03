import { useEffect, useMemo, useState } from "react";
import $ from "jquery";
import { CCP4i2ReportElementProps } from "./CCP4i2ReportElement";
import {
  Accordion,
  AccordionSummary,
  AccordionDetails,
  Typography,
  Box,
} from "@mui/material";
import { useTheme } from "@mui/material/styles";
import ExpandMoreIcon from "@mui/icons-material/ExpandMore";
import { JobMenu } from "../../providers/job-context-menu";
import { CCP4i2ReportFile } from "./CCP4i2ReportFile";

export const CCP4i2ReportInputOutputData: React.FC<CCP4i2ReportElementProps> = (
  props
) => {
  const theme = useTheme();
  const [fileUUIDs, setFileUUIDs] = useState<string[]>([]);
  const [expanded, setExpanded] = useState(true);

  useEffect(() => {
    if (!props.item) return;
    const fileUUIDs: string[] = [];
    $(props.item)
      .find("div")
      .each((iDiv, div) => {
        const divId = $(div).attr("id");
        if (divId) {
          const matches = /input_file_(.*)/.exec(divId);
          if (matches) {
            fileUUIDs.push(matches[1]);
          }
        }
      });
    setFileUUIDs(fileUUIDs);
  }, [props.item]);

  const title = useMemo<string>(() => {
    const h5Nodes = $(props.item).find("h5");
    const h5s = h5Nodes
      .map((iItem, item) => {
        return $(item).text();
      })
      .toArray();
    return h5s.length > 0 ? h5s.join(", ") : "Input or Output data";
  }, [props.item]);

  const handleAccordionChange = (
    _event: React.SyntheticEvent,
    isExpanded: boolean
  ) => {
    setExpanded(isExpanded);
  };

  return (
    <>
      <Accordion
        expanded={expanded}
        onChange={handleAccordionChange}
        disableGutters
        elevation={1}
      >
        <AccordionSummary
          expandIcon={<ExpandMoreIcon />}
          aria-controls="input-output-content"
          id="input-output-header"
          sx={{
            backgroundColor: theme.palette.action.hover,
            "&:hover": {
              backgroundColor: theme.palette.action.selected,
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
          <Box display="flex" alignItems="center" width="100%">
            <Typography variant="subtitle1" sx={{ fontWeight: 500 }}>
              {title}
            </Typography>
            {fileUUIDs.length > 0 && (
              <Typography
                variant="body2"
                color="textSecondary"
                sx={{ ml: "auto", mr: 2 }}
              >
                {fileUUIDs.length} file{fileUUIDs.length !== 1 ? "s" : ""}
              </Typography>
            )}
          </Box>
        </AccordionSummary>

        <AccordionDetails sx={{ p: 2 }}>
          {fileUUIDs.length > 0 ? (
            fileUUIDs.map((fileUUID: string, iFile: number) => (
              <CCP4i2ReportFile {...props} uuid={fileUUID} key={iFile} />
            ))
          ) : (
            <Typography variant="body2" color="textSecondary">
              No files found
            </Typography>
          )}
        </AccordionDetails>
      </Accordion>
      <JobMenu />
    </>
  );
};
