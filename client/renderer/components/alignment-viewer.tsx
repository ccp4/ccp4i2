import React, { useMemo, useCallback, useState, useRef } from "react";
import {
  Box,
  Typography,
  Chip,
  IconButton,
  Stack,
  Switch,
  FormControlLabel,
  TextField,
  Accordion,
  AccordionSummary,
  AccordionDetails,
  Alert,
  Button,
} from "@mui/material";
import {
  ExpandLess,
  ExpandMore,
  CheckCircle,
  Palette,
  PaletteOutlined,
  HighlightAlt,
  Clear,
  KeyboardArrowUp,
} from "@mui/icons-material";
import { Editor } from "@monaco-editor/react";
import { useRunCheck } from "../providers/run-check-provider";
import { useJob } from "../utils";
import { Job } from "../types/models";
import { useCCP4i2Window } from "../app-context";
import { useTheme } from "../theme/theme-provider";

interface AlignmentViewerProps {
  alignment: string;
}

// Clustal color scheme for amino acids
const AMINO_ACID_COLORS: Record<string, string> = {
  A: "#80a0f0",
  R: "#f01505",
  N: "#00ff00",
  D: "#c048c0",
  C: "#f08080",
  Q: "#00ff00",
  E: "#c048c0",
  G: "#f09048",
  H: "#15a4a4",
  I: "#80a0f0",
  L: "#80a0f0",
  K: "#f01505",
  M: "#80a0f0",
  F: "#80a0f0",
  P: "#ffff00",
  S: "#00ff00",
  T: "#00ff00",
  W: "#80a0f0",
  Y: "#15a4a4",
  V: "#80a0f0",
  "-": "#ffffff",
  ".": "#ffffff",
  " ": "#ffffff",
};

interface SequenceBlock {
  sequences: Array<{
    name: string;
    sequence: string;
    startPos: number;
    endPos: number;
  }>;
  conservation: string;
  blockStart: number;
  blockEnd: number;
}

interface HighlightRegion {
  sequenceIndex: number;
  startPos: number;
  endPos: number;
  valid: boolean;
  error?: string;
}

export const AlignmentViewer: React.FC<AlignmentViewerProps> = ({
  alignment,
}) => {
  const { customColors } = useTheme();
  const AMINO_ACID_COLORS = customColors.aminoAcids;
  const [showColors, setShowColors] = useState(true);
  const [blockSize, setBlockSize] = useState(60);
  const [highlightInput, setHighlightInput] = useState("");
  const [showHighlightPanel, setShowHighlightPanel] = useState(false);
  const alignmentContainerRef = useRef<HTMLDivElement>(null);

  const { sequences, sequenceNames, clustalConservation, sequenceBlocks } =
    useMemo(() => {
      if (!alignment)
        return {
          sequences: [],
          sequenceNames: [],
          clustalConservation: "",
          sequenceBlocks: [],
        };

      const lines = alignment.trim().split("\n");
      const seqMap: Record<string, string> = {};
      const names: string[] = [];
      let conservation = "";

      // Skip the ClustalW header line (usually starts with "CLUSTAL")
      let startIndex = 0;
      for (let i = 0; i < lines.length; i++) {
        if (lines[i].trim().startsWith("CLUSTAL") || lines[i].trim() === "") {
          startIndex = i + 1;
        } else {
          break;
        }
      }

      // Parse sequence blocks
      let seqStartPos, seqBlockLength;
      for (let i = startIndex; i < lines.length; i++) {
        const line = lines[i];

        // Skip empty lines
        if (!line.trim()) continue;

        // Check if this is a conservation line (starts with spaces/symbols, no sequence name)
        if (line.trim().match(/^[\s\*\:\.\s]+$/)) {
          //const seqStartPos = line.search(/\S/); // Find position of first non-whitespace character
          conservation += line.substring(
            seqStartPos,
            seqStartPos + seqBlockLength
          ); // Take from sequence start position
          continue;
        }

        // Parse sequence line: "sequence_name    SEQUENCE_DATA"
        const trimmedLine = line.trim();
        const parts = trimmedLine.split(/\s+/);
        if (parts.length >= 2) {
          const seqName = parts[0];
          // Preserve trailing spaces by taking everything after the sequence name
          const seqDataStart = line.indexOf(parts[1]);
          seqStartPos = seqDataStart;
          const seqData = line.substring(seqDataStart);
          seqBlockLength = seqData.length;

          // Initialize sequence if first occurrence
          if (!seqMap[seqName]) {
            seqMap[seqName] = "";
            names.push(seqName);
          }

          seqMap[seqName] += seqData;
        }
      }

      const seqs = names.map((name) => seqMap[name]);
      const maxLength = Math.max(...seqs.map((s) => s.length), 0);

      // Create blocks using the current block size
      const blocks: SequenceBlock[] = [];

      for (let start = 0; start < maxLength; start += blockSize) {
        const end = Math.min(start + blockSize, maxLength);

        const blockSequences = seqs.map((seq, index) => {
          const blockSeq = seq.slice(start, end);

          // Calculate actual sequence positions (ignoring gaps)
          const beforeBlock = seq.slice(0, start);
          const startPos = beforeBlock.replace(/[-.\s]/g, "").length + 1;

          const endPos = seq.slice(0, end).replace(/[-.\s]/g, "").length;

          return {
            name: names[index],
            sequence: blockSeq,
            startPos: startPos,
            endPos: endPos || startPos - 1, // Handle case where block is all gaps
          };
        });

        const blockConservation = conservation
          .slice(start, end)
          .padEnd(end - start, " ");

        blocks.push({
          sequences: blockSequences,
          conservation: blockConservation,
          blockStart: start + 1,
          blockEnd: end,
        });
      }

      return {
        sequences: seqs,
        sequenceNames: names,
        clustalConservation: conservation,
        sequenceBlocks: blocks,
      };
    }, [alignment, blockSize]);

  // Parse highlight regions from input
  const highlightRegions = useMemo((): HighlightRegion[] => {
    if (!highlightInput.trim()) return [];

    return highlightInput
      .split(",")
      .map((entry) => entry.trim())
      .filter((entry) => entry.length > 0)
      .map((entry) => {
        const match = entry.match(/^(\d+):(\d+)-(\d+)$/);
        if (!match) {
          return {
            sequenceIndex: -1,
            startPos: -1,
            endPos: -1,
            valid: false,
            error: `Invalid format: "${entry}". Use format "N:start-end"`,
          };
        }

        const sequenceIndex = parseInt(match[1]) - 1; // Convert to 0-based
        const startPos = parseInt(match[2]);
        const endPos = parseInt(match[3]);

        // Validate sequence index
        if (sequenceIndex < 0 || sequenceIndex >= sequences.length) {
          return {
            sequenceIndex,
            startPos,
            endPos,
            valid: false,
            error: `Sequence ${
              sequenceIndex + 1
            } does not exist. Available sequences: 1-${sequences.length}`,
          };
        }

        // Validate position range
        if (startPos > endPos) {
          return {
            sequenceIndex,
            startPos,
            endPos,
            valid: false,
            error: `Invalid range: start position (${startPos}) > end position (${endPos})`,
          };
        }

        const maxSeqPos = sequences[sequenceIndex].replace(
          /[-.\s]/g,
          ""
        ).length;
        if (startPos < 1 || endPos > maxSeqPos) {
          return {
            sequenceIndex,
            startPos,
            endPos,
            valid: false,
            error: `Position out of range for sequence ${
              sequenceIndex + 1
            }. Valid range: 1-${maxSeqPos}`,
          };
        }

        return {
          sequenceIndex,
          startPos,
          endPos,
          valid: true,
        };
      });
  }, [highlightInput, sequences]);

  // Convert sequence position to alignment position
  const getAlignmentPositions = useCallback(
    (
      sequenceIndex: number,
      seqStart: number,
      seqEnd: number
    ): { start: number; end: number } | null => {
      if (sequenceIndex < 0 || sequenceIndex >= sequences.length) return null;

      const sequence = sequences[sequenceIndex];
      let currentSeqPos = 0;
      let alignmentStart = -1;
      let alignmentEnd = -1;

      for (let i = 0; i < sequence.length; i++) {
        const residue = sequence[i];
        if (residue !== "-" && residue !== "." && residue !== " ") {
          currentSeqPos++;
          if (currentSeqPos === seqStart && alignmentStart === -1) {
            alignmentStart = i;
          }
          if (currentSeqPos === seqEnd) {
            alignmentEnd = i;
            break;
          }
        }
      }

      if (alignmentStart === -1 || alignmentEnd === -1) return null;
      return { start: alignmentStart, end: alignmentEnd };
    },
    [sequences]
  );

  // Get highlight regions for a specific sequence and block
  const getBlockHighlightRegions = useCallback(
    (sequenceIndex: number, blockStart: number, blockEnd: number) => {
      const validRegions = highlightRegions.filter(
        (region) => region.valid && region.sequenceIndex === sequenceIndex
      );

      return validRegions
        .map((region) => {
          const alignmentPositions = getAlignmentPositions(
            region.sequenceIndex,
            region.startPos,
            region.endPos
          );

          if (!alignmentPositions) return null;

          // Check if region overlaps with this block
          // blockStart and blockEnd are 1-based, alignmentPositions are 0-based
          const blockStartIndex = blockStart - 1; // Convert to 0-based
          const blockEndIndex = blockEnd - 1; // Convert to 0-based

          if (
            alignmentPositions.end < blockStartIndex ||
            alignmentPositions.start > blockEndIndex
          ) {
            return null;
          }

          // Calculate positions relative to the block
          const blockRelativeStart = Math.max(
            0,
            alignmentPositions.start - blockStartIndex
          );
          const blockRelativeEnd = Math.min(
            blockEndIndex - blockStartIndex,
            alignmentPositions.end - blockStartIndex
          );

          return {
            start: blockRelativeStart,
            end: blockRelativeEnd,
            region,
          };
        })
        .filter(Boolean);
    },
    [highlightRegions, getAlignmentPositions]
  );

  const calculateConservation = useCallback(
    (position: number): number => {
      if (sequences.length === 0) return 0;

      const residues: Record<string, number> = {};
      let validResidues = 0;

      sequences.forEach((seq) => {
        const residue = seq[position];
        if (residue && residue !== "-" && residue !== "." && residue !== " ") {
          residues[residue] = (residues[residue] || 0) + 1;
          validResidues++;
        }
      });

      if (validResidues === 0) return 0;

      const maxCount = Math.max(...Object.values(residues));
      return maxCount / validResidues;
    },
    [sequences]
  );

  const getConservationColor = useCallback(
    (symbol: string): string => {
      // Return white background when colors are disabled
      if (!showColors) return "#ffffff";

      switch (symbol) {
        case "*":
          return "#ff0000"; // Fully conserved (red)
        case ":":
          return "#0000ff"; // Strongly similar (blue)
        case ".":
          return "#00aa00"; // Weakly similar (green)
        default:
          return "#ffffff"; // No conservation (white)
      }
    },
    [showColors]
  );

  const formatSequenceName = useCallback(
    (name: string, maxWidth: number = 20): string => {
      // Reduced from 25 to 20
      if (name.length <= maxWidth) return name.padEnd(maxWidth);
      return name.slice(0, maxWidth - 3) + "...";
    },
    []
  );

  const getResidueBackgroundColor = useCallback(
    (residue: string): string => {
      if (!showColors) return "#ffffff";
      return AMINO_ACID_COLORS[residue.toUpperCase()] || "#ffffff";
    },
    [showColors]
  );

  const getResidueTextColor = useCallback(
    (residue: string): string => {
      if (!showColors) {
        return residue === "-" || residue === "." ? "#999" : "#000";
      }
      return residue === "-" || residue === "." ? "#999" : "#000";
    },
    [showColors]
  );

  const scrollToTop = useCallback(() => {
    if (alignmentContainerRef.current) {
      alignmentContainerRef.current.scrollTop = 0;
    }
  }, []);

  const validHighlights = highlightRegions.filter((r) => r.valid);
  const invalidHighlights = highlightRegions.filter((r) => !r.valid);

  if (!alignment) {
    return (
      <Box sx={{ p: 2, textAlign: "center", color: "text.secondary" }}>
        <Typography>No alignment data available</Typography>
      </Box>
    );
  }

  return (
    <Box sx={{ width: "100%", p: 0.5 }}>
      {" "}
      {/* Reduced padding from 1 to 0.5 */}
      <Box
        sx={{
          display: "flex",
          justifyContent: "space-between",
          alignItems: "center",
          mb: 1, // Reduced from 2 to 1
        }}
      >
        <Typography variant="h6">
          Protein Sequence Alignment (ClustalW)
        </Typography>

        <Box sx={{ display: "flex", gap: 1, alignItems: "center" }}>
          {" "}
          {/* Reduced gap from 2 to 1 */}
          {/* Block size selector */}
          <Box sx={{ display: "flex", gap: 0.5 }}>
            {" "}
            {/* Reduced gap from 1 to 0.5 */}
            {[40, 60, 80, 120].map((size) => (
              <Chip
                key={size}
                label={`${size}`}
                size="small"
                color={blockSize === size ? "primary" : "default"}
                variant={blockSize === size ? "filled" : "outlined"}
                onClick={() => setBlockSize(size)}
                sx={{ cursor: "pointer" }}
              />
            ))}
          </Box>
          {/* Highlight toggle */}
          <Chip
            icon={<HighlightAlt />}
            label={`Highlight (${validHighlights.length})`}
            size="small"
            color={showHighlightPanel ? "primary" : "default"}
            variant={showHighlightPanel ? "filled" : "outlined"}
            onClick={() => setShowHighlightPanel(!showHighlightPanel)}
            sx={{ cursor: "pointer" }}
          />
          {/* Color toggle */}
          <FormControlLabel
            control={
              <Switch
                checked={showColors}
                onChange={(e) => setShowColors(e.target.checked)}
                size="small"
              />
            }
            label={
              <Box sx={{ display: "flex", alignItems: "center", gap: 0.5 }}>
                {showColors ? (
                  <Palette fontSize="small" />
                ) : (
                  <PaletteOutlined fontSize="small" />
                )}
                <Typography variant="caption">Colors</Typography>
              </Box>
            }
          />
        </Box>
      </Box>
      {/* Highlight Panel */}
      {showHighlightPanel && (
        <Accordion expanded sx={{ mb: 1 }}>
          {" "}
          {/* Reduced from mb: 2 to mb: 1 */}
          <AccordionSummary>
            <Typography variant="subtitle2">
              Sequence Region Highlighting
            </Typography>
          </AccordionSummary>
          <AccordionDetails>
            <Box sx={{ display: "flex", flexDirection: "column", gap: 2 }}>
              <Box>
                <Typography variant="body2" sx={{ mb: 1 }}>
                  Enter regions to highlight (comma-separated):
                </Typography>
                <TextField
                  fullWidth
                  multiline
                  rows={3}
                  value={highlightInput}
                  onChange={(e) => setHighlightInput(e.target.value)}
                  placeholder="1:23-45, 2:10-20, 1:100-150"
                  size="small"
                  helperText="Format: sequenceNumber:startPosition-endPosition (e.g., 1:23-45 highlights residues 23-45 of sequence 1)"
                  InputProps={{
                    endAdornment: highlightInput && (
                      <IconButton
                        size="small"
                        onClick={() => setHighlightInput("")}
                        sx={{ mt: -3 }}
                      >
                        <Clear fontSize="small" />
                      </IconButton>
                    ),
                  }}
                />
              </Box>

              {/* Quick test buttons */}
              <Box>
                <Typography variant="body2" sx={{ mb: 1 }}>
                  Quick tests:
                </Typography>
                <Box sx={{ display: "flex", gap: 1, flexWrap: "wrap" }}>
                  <Chip
                    label="Test 1:1-5"
                    size="small"
                    color="primary"
                    variant="outlined"
                    onClick={() => {
                      setHighlightInput("1:1-5");
                      setTimeout(scrollToTop, 100);
                    }}
                    sx={{ cursor: "pointer" }}
                  />
                  <Button
                    size="small"
                    startIcon={<KeyboardArrowUp />}
                    onClick={scrollToTop}
                    variant="outlined"
                  >
                    Jump to Top
                  </Button>
                </Box>
              </Box>

              {/* Sequence reference */}
              <Box>
                <Typography variant="body2" sx={{ mb: 1 }}>
                  Available sequences:
                </Typography>
                <Box sx={{ display: "flex", flexWrap: "wrap", gap: 0.5 }}>
                  {sequenceNames.map((name, index) => (
                    <Chip
                      key={index}
                      label={`${index + 1}: ${name}`}
                      size="small"
                      variant="outlined"
                    />
                  ))}
                </Box>
              </Box>

              {/* Debug info */}
              {process.env.NODE_ENV === "development" &&
                validHighlights.length > 0 && (
                  <Box
                    sx={{ p: 1, backgroundColor: "#f5f5f5", borderRadius: 1 }}
                  >
                    <Typography variant="caption">Debug Info:</Typography>
                    {validHighlights.map((region, index) => {
                      const alignmentPos = getAlignmentPositions(
                        region.sequenceIndex,
                        region.startPos,
                        region.endPos
                      );
                      return (
                        <Box
                          key={index}
                          sx={{ fontSize: "10px", fontFamily: "monospace" }}
                        >
                          Seq {region.sequenceIndex + 1}: {region.startPos}-
                          {region.endPos} → Alignment: {alignmentPos?.start}-
                          {alignmentPos?.end}
                        </Box>
                      );
                    })}
                  </Box>
                )}

              {/* Validation messages */}
              {invalidHighlights.length > 0 && (
                <Box>
                  {invalidHighlights.map((region, index) => (
                    <Alert key={index} severity="error" sx={{ mb: 1 }}>
                      {region.error}
                    </Alert>
                  ))}
                </Box>
              )}

              {/* Valid highlights summary */}
              {validHighlights.length > 0 && (
                <Box>
                  <Typography variant="body2" sx={{ mb: 1 }}>
                    Active highlights:
                  </Typography>
                  <Box sx={{ display: "flex", flexWrap: "wrap", gap: 0.5 }}>
                    {validHighlights.map((region, index) => (
                      <Chip
                        key={index}
                        label={`${sequenceNames[region.sequenceIndex]}: ${
                          region.startPos
                        }-${region.endPos}`}
                        size="small"
                        color="warning"
                        variant="filled"
                      />
                    ))}
                  </Box>
                </Box>
              )}
            </Box>
          </AccordionDetails>
        </Accordion>
      )}
      <Box
        ref={alignmentContainerRef}
        sx={{
          border: "1px solid #ddd",
          borderRadius: 1,
          overflow: "auto",
          maxHeight: "700px", // Increased from 600px to 700px
          fontFamily: "Courier, monospace",
          fontSize: "13px", // Increased from 11px to 13px
          lineHeight: 1.2, // Reduced from 1.3 to 1.2
          backgroundColor: "#fafafa",
          p: 0.5, // Reduced from 1 to 0.5
        }}
      >
        {sequenceBlocks.map((block, blockIndex) => (
          <Box key={blockIndex} sx={{ mb: 2, pageBreakInside: "avoid" }}>
            {" "}
            {/* Reduced from mb: 3 to mb: 2 */}
            {/* Block header with position range */}
            <Box
              sx={{
                mb: 0.5, // Reduced from 1 to 0.5
                p: 0.3, // Reduced from 0.5 to 0.3
                backgroundColor: "#e8e8e8",
                borderRadius: 0.5,
                textAlign: "center",
                fontWeight: "bold",
                fontSize: "10px",
              }}
            >
              Alignment positions {block.blockStart} - {block.blockEnd}
            </Box>
            {/* Position ruler */}
            <Box sx={{ display: "flex", mb: 0.3 }}>
              {" "}
              {/* Reduced from mb: 0.5 to mb: 0.3 */}
              <Box sx={{ width: "160px" }} />{" "}
              {/* Reduced from 220px to 160px */}
              <Box sx={{ width: "40px" }} /> {/* Reduced from 50px to 40px */}
              <Box
                sx={{
                  fontFamily: "Courier, monospace",
                  fontSize: "9px", // Increased from 8px to 9px
                  display: "flex",
                }}
              >
                {Array.from(
                  { length: block.sequences[0]?.sequence.length || 0 },
                  (_, i) => (
                    <span
                      key={i}
                      style={{ width: "12px", textAlign: "center" }} // Increased from 10px to 12px
                    >
                      {(block.blockStart + i) % 10 === 0
                        ? (block.blockStart + i) % 100
                        : ""}
                    </span>
                  )
                )}
              </Box>
              <Box sx={{ width: "40px" }} /> {/* Reduced from 50px to 40px */}
            </Box>
            {/* Sequence rows */}
            {block.sequences.map((seqData, seqIndex) => {
              const highlightRegionsRaw = getBlockHighlightRegions(
                seqIndex,
                block.blockStart,
                block.blockEnd
              );
              const highlightRegions = highlightRegionsRaw.filter(
                (highlight): highlight is NonNullable<typeof highlight> =>
                  highlight !== null
              );

              return (
                <Box
                  key={seqIndex}
                  sx={{ display: "flex", alignItems: "center", mb: 0.1 }} // Reduced from mb: 0.2 to mb: 0.1
                >
                  {/* Sequence name */}
                  <Box
                    sx={{
                      width: "160px", // Reduced from 220px to 160px
                      backgroundColor: "#f8f8f8",
                      border: "1px solid #ddd",
                      p: 0.2, // Reduced from 0.3 to 0.2
                      fontSize: "10px", // Increased from 9px to 10px
                      textAlign: "right",
                      overflow: "hidden",
                      textOverflow: "ellipsis",
                      whiteSpace: "nowrap",
                    }}
                    title={seqData.name}
                  >
                    {formatSequenceName(seqData.name, 20)}{" "}
                    {/* Reduced max width to 20 */}
                  </Box>

                  {/* Start position */}
                  <Box
                    sx={{
                      width: "40px", // Reduced from 50px to 40px
                      textAlign: "right",
                      pr: 0.5, // Reduced from 1 to 0.5
                      fontSize: "10px", // Increased from 9px to 10px
                      color: "#666",
                    }}
                  >
                    {seqData.startPos}
                  </Box>

                  {/* Sequence data with highlight overlays */}
                  <Box
                    sx={{
                      position: "relative",
                      display: "flex",
                      fontFamily: "Courier, monospace",
                    }}
                  >
                    {/* Highlight overlays */}
                    {highlightRegions.map((highlight, highlightIndex) => (
                      <Box
                        key={highlightIndex}
                        sx={{
                          position: "absolute",
                          left: `${highlight.start * 12}px`, // Increased from 10px to 12px
                          width: `${
                            (highlight.end - highlight.start + 1) * 12 // Increased from 10px to 12px
                          }px`,
                          height: "16px", // Increased from 14px to 16px
                          backgroundColor: "rgba(255, 107, 53, 0.7)",
                          border: "5px solid #000000",
                          borderRadius: "3px",
                          pointerEvents: "none",
                          zIndex: 1,
                          top: 0,
                          boxShadow: "0 0 4px rgba(0, 0, 0, 1.0)",
                        }}
                        title={`Highlighted region: ${highlight.region.startPos}-${highlight.region.endPos}`}
                      />
                    ))}

                    {/* Sequence residues */}
                    {seqData.sequence.split("").map((residue, pos) => {
                      // Calculate the actual sequence position for this residue
                      const alignmentPos = block.blockStart + pos - 1; // 0-based alignment position
                      const sequencePos = sequences[seqIndex]
                        .slice(0, alignmentPos + 1)
                        .replace(/[-.\s]/g, "").length;

                      return (
                        <Box
                          key={pos}
                          sx={{
                            width: "12px", // Increased from 10px to 12px
                            height: "16px", // Increased from 14px to 16px
                            backgroundColor: getResidueBackgroundColor(residue),
                            color: getResidueTextColor(residue),
                            display: "flex",
                            alignItems: "center",
                            justifyContent: "center",
                            fontSize: "11px", // Increased from 9px to 11px
                            fontWeight: "bold",
                            border: showColors
                              ? "0.5px solid rgba(0,0,0,0.1)"
                              : "0.5px solid rgba(0,0,0,0.3)",
                            cursor: "pointer",
                            position: "relative",
                            zIndex: 2,
                            "&:hover": {
                              outline: "1px solid #1976d2",
                              zIndex: 3,
                            },
                          }}
                          title={`${seqData.name} - Alignment pos: ${
                            block.blockStart + pos
                          }, Residue: ${residue}, Seq pos: ${
                            residue === "-" ||
                            residue === "." ||
                            residue === " "
                              ? "gap"
                              : sequencePos
                          }`}
                        >
                          {residue}
                        </Box>
                      );
                    })}
                  </Box>

                  {/* End position */}
                  <Box
                    sx={{
                      width: "40px", // Reduced from 50px to 40px
                      textAlign: "left",
                      pl: 0.5, // Reduced from 1 to 0.5
                      fontSize: "10px", // Increased from 9px to 10px
                      color: "#666",
                    }}
                  >
                    {seqData.endPos}
                  </Box>
                </Box>
              );
            })}
            {/* Conservation line */}
            <Box sx={{ display: "flex", alignItems: "center", mt: 0.3 }}>
              {" "}
              {/* Reduced from mt: 0.5 to mt: 0.3 */}
              <Box
                sx={{
                  width: "160px", // Reduced from 220px to 160px
                  textAlign: "right",
                  pr: 0.5, // Reduced from 1 to 0.5
                  fontSize: "10px", // Increased from 9px to 10px
                  color: "#666",
                  fontStyle: "italic",
                }}
              >
                Conservation
              </Box>
              <Box sx={{ width: "40px" }} /> {/* Reduced from 50px to 40px */}
              <Box sx={{ display: "flex" }}>
                {block.conservation.split("").map((symbol, pos) => (
                  <Box
                    key={pos}
                    sx={{
                      width: "12px", // Increased from 10px to 12px
                      height: "16px", // Increased from 14px to 16px
                      backgroundColor: getConservationColor(symbol),
                      border: "0.5px solid rgba(0,0,0,0.1)",
                      display: "flex",
                      alignItems: "center",
                      justifyContent: "center",
                      fontSize: "9px", // Increased from 8px to 9px
                      fontWeight: "bold",
                      color: !showColors ? "#000" : "#000",
                    }}
                    title={`Position ${block.blockStart + pos}: ${
                      symbol === "*"
                        ? "Fully conserved"
                        : symbol === ":"
                          ? "Strongly similar"
                          : symbol === "."
                            ? "Weakly similar"
                            : "Variable"
                    }`}
                  >
                    {symbol}
                  </Box>
                ))}
              </Box>
              <Box sx={{ width: "40px" }} /> {/* Reduced from 50px to 40px */}
            </Box>
          </Box>
        ))}
      </Box>
      {/* Summary and Legend - Keep existing but with reduced margin */}
      <Box sx={{ mt: 0.5 }}>
        {" "}
        {/* Reduced from mt: 1 to mt: 0.5 */}
        <Box sx={{ display: "flex", gap: 1, flexWrap: "wrap", mb: 1 }}>
          <Chip
            size="small"
            label={`${sequences.length} sequences`}
            icon={<CheckCircle />}
          />
          <Chip
            size="small"
            label={`${sequenceBlocks.length} blocks`}
            color="secondary"
            variant="outlined"
          />
          <Chip
            size="small"
            label={`${blockSize} chars/block`}
            color="primary"
            variant="outlined"
          />
          {validHighlights.length > 0 && (
            <Chip
              size="small"
              label={`${validHighlights.length} highlighted regions`}
              color="warning"
              variant="filled"
            />
          )}
        </Box>
        {/* Conservation Legend */}
        <Box sx={{ fontSize: "11px", color: "text.secondary" }}>
          <Typography variant="caption" display="block" gutterBottom>
            Conservation symbols:
          </Typography>
          <Box sx={{ display: "flex", gap: 2, flexWrap: "wrap" }}>
            <Box sx={{ display: "flex", alignItems: "center", gap: 0.5 }}>
              <Box
                sx={{
                  width: "12px",
                  height: "12px",
                  backgroundColor: showColors ? "#ff0000" : "#ffffff",
                  border: "1px solid #ccc",
                }}
              />
              <Typography variant="caption">* Fully conserved</Typography>
            </Box>
            <Box sx={{ display: "flex", alignItems: "center", gap: 0.5 }}>
              <Box
                sx={{
                  width: "12px",
                  height: "12px",
                  backgroundColor: showColors ? "#0000ff" : "#ffffff",
                  border: "1px solid #ccc",
                }}
              />
              <Typography variant="caption">: Strongly similar</Typography>
            </Box>
            <Box sx={{ display: "flex", alignItems: "center", gap: 0.5 }}>
              <Box
                sx={{
                  width: "12px",
                  height: "12px",
                  backgroundColor: showColors ? "#00aa00" : "#ffffff",
                  border: "1px solid #ccc",
                }}
              />
              <Typography variant="caption">. Weakly similar</Typography>
            </Box>
            {validHighlights.length > 0 && (
              <Box sx={{ display: "flex", alignItems: "center", gap: 0.5 }}>
                <Box
                  sx={{
                    width: "20px",
                    height: "12px",
                    backgroundColor: "rgba(255, 107, 53, 0.7)",
                    border: "5px solid #000000",
                    borderRadius: "3px",
                    boxShadow: "0 0 4px rgba(0, 0, 0, 1.0)",
                  }}
                />
                <Typography variant="caption">Highlighted regions</Typography>
              </Box>
            )}
          </Box>
        </Box>
      </Box>
    </Box>
  );
};
