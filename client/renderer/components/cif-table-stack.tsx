import React from "react";
import {
  Accordion,
  AccordionSummary,
  AccordionDetails,
  Stack,
  Typography,
  Table,
  TableBody,
  TableCell,
  TableRow,
} from "@mui/material";
import ExpandMoreIcon from "@mui/icons-material/ExpandMore";

// Minimal CIF parser for blocks and loops
function parseCIFBlocks(cifText: string): Record<string, any> {
  const blocks: Record<string, any> = {};
  const blockRegex = /(?:^|\n)data_(\S+)[\s\S]*?(?=(?:\n\s*data_\S+)|$)/g;
  let match: RegExpExecArray | null;

  while ((match = blockRegex.exec(cifText)) !== null) {
    const blockName = match[1];
    const blockContent = match[0]
      .replace(new RegExp(`^data_${blockName}`), "")
      .trim();
    blocks[blockName] = parseCIFLoops(blockContent);
  }
  return blocks;
}

function parseCIFLoops(
  blockText: string
): Array<{ keys: string[]; values: string[][] }> {
  const lines = blockText.split(/\r?\n/);
  const loops: Array<{ keys: string[]; values: string[][] }> = [];
  let i = 0;

  while (i < lines.length) {
    let line = lines[i].trim();
    if (line.startsWith("loop_")) {
      i++;
      const keys: string[] = [];
      while (i < lines.length && lines[i].trim().startsWith("_")) {
        keys.push(lines[i].trim());
        i++;
      }
      const values: string[][] = [];
      while (
        i < lines.length &&
        lines[i].trim() &&
        !lines[i].trim().startsWith("_") &&
        !lines[i].trim().startsWith("loop_")
      ) {
        values.push(lines[i].trim().split(/\s+/));
        i++;
      }
      loops.push({ keys, values });
      continue;
    }
    i++;
  }
  return loops;
}

export interface CifTableStackProps {
  cifText: string;
}

export const CifTableStack: React.FC<CifTableStackProps> = ({ cifText }) => {
  const blocks = parseCIFBlocks(cifText);

  return (
    <Stack spacing={2}>
      {Object.entries(blocks).map(([blockName, loops]) => (
        <Accordion key={blockName} defaultExpanded>
          <AccordionSummary expandIcon={<ExpandMoreIcon />}>
            <Typography variant="subtitle1">{`data_${blockName}`}</Typography>
          </AccordionSummary>
          <AccordionDetails>
            <Stack spacing={2}>
              {(loops as Array<{ keys: string[]; values: string[][] }>).map(
                (loop, idx) => (
                  <Table
                    key={idx}
                    size="small"
                    sx={{
                      mb: 2,
                      width: "100%",
                      tableLayout: "fixed",
                    }}
                  >
                    <TableBody>
                      <TableRow>
                        {loop.keys.map((key) => (
                          <TableCell
                            key={key}
                            variant="head"
                            sx={{
                              fontWeight: "bold",
                              wordBreak: "break-word",
                              maxWidth: 150,
                            }}
                          >
                            {key}
                          </TableCell>
                        ))}
                      </TableRow>
                      {loop.values.map((row, ridx) => (
                        <TableRow key={ridx}>
                          {row.map((val, cid) => (
                            <TableCell
                              key={cid}
                              sx={{
                                wordBreak: "break-word",
                                maxWidth: 150,
                              }}
                            >
                              {val}
                            </TableCell>
                          ))}
                        </TableRow>
                      ))}
                    </TableBody>
                  </Table>
                )
              )}
              {(!loops || (loops as Array<any>).length === 0) && (
                <Typography variant="body2" color="text.secondary">
                  No loops found in this block.
                </Typography>
              )}
            </Stack>
          </AccordionDetails>
        </Accordion>
      ))}
    </Stack>
  );
};
