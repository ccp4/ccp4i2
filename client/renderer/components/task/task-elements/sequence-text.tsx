import { Box, Tooltip } from "@mui/material";
import {
  PolymerType,
  describeOddLetters,
  nonCanonicalLetters,
  previewWindow,
  sequenceSegments,
} from "../../../lib/sequence-letters";

/**
 * A sequence with the letters the server will reject marked: red, bold and
 * underlined, each with its position in a tooltip. `maxLength` makes a
 * one-line preview, centred on the first rejected letter when there is
 * one past the start, so the preview shows the problem.
 */
export const SequenceText: React.FC<{
  sequence: string;
  polymerType: PolymerType;
  maxLength?: number;
  wrap?: boolean;
  /** Called with the 1-based position (in the whitespace-free sequence) of a marked letter that was clicked. */
  onLetterClick?: (position: number) => void;
}> = ({ sequence, polymerType, maxLength, wrap = false, onLetterClick }) => {
  const seq = (sequence || "").replace(/\s/g, "");
  if (!seq) return null;
  const window = maxLength ? previewWindow(seq, polymerType, maxLength) : null;
  const text = window ? window.text : seq;
  const offset = window ? window.start - 1 : 0;
  const segments = sequenceSegments(text, polymerType);
  return (
    <Box
      component="span"
      sx={{
        fontFamily: "monospace",
        whiteSpace: wrap ? "pre-wrap" : "nowrap",
        wordBreak: wrap ? "break-all" : undefined,
      }}
    >
      {window?.truncatedBefore ? "\u2026" : ""}
      {segments.map((s) =>
        s.canonical ? (
          <span key={s.start}>{s.text}</span>
        ) : (
          <Tooltip
            key={s.start}
            title={`Not a ${(polymerType || "").toUpperCase() === "PROTEIN" ? "standard amino acid" : "nucleotide"} letter: ${describeOddLetters(
              nonCanonicalLetters(s.text, polymerType).map((o) => ({ ...o, position: o.position + offset + s.start - 1 })),
              5,
            )}${onLetterClick ? ". Click to select it in the editor" : ""}`}
          >
            <Box
              component="span"
              onClick={onLetterClick ? () => onLetterClick(offset + s.start) : undefined}
              sx={{
                color: "error.main",
                fontWeight: "bold",
                textDecoration: "underline",
                textDecorationStyle: "wavy",
                bgcolor: "rgba(211, 47, 47, 0.12)",
                borderRadius: 0.5,
                cursor: onLetterClick ? "pointer" : undefined,
              }}
            >
              {s.text}
            </Box>
          </Tooltip>
        ),
      )}
      {window?.truncatedAfter ? "\u2026" : ""}
    </Box>
  );
};
