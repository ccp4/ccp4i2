import React, { useCallback, useEffect, useState } from "react";
import {
  Alert,
  Box,
  Button,
  Chip,
  IconButton,
  MenuItem,
  Paper,
  Select,
  TextField,
  ToggleButton,
  ToggleButtonGroup,
  Tooltip,
  Typography,
  useTheme,
} from "@mui/material";
import { Add, Delete } from "@mui/icons-material";

import {
  Assembly,
  ChainInfo,
  IMPLICIT_ROLE,
  MODES,
  MODE_HELP,
  Mode,
  ParsedBody,
  Segment,
  referenceBounds,
  roleLabel,
} from "./dm-spec";
import { PreviewBody, rmsdTone } from "./use-ncs-preview";
import { bodyColour } from "./coverage-strip";

/**
 * Rigid bodies, as picks rather than as a mini-language.
 *
 * A body is a set of residue ranges that move together. Typed as
 * "cyclin:10-95,CDK:45-60" it needs a paragraph to explain, and the role has
 * to match a word in the assembly exactly or nothing lines up. Here the role
 * is a dropdown of the assembly's own columns and the residue numbers are
 * bounded by the chain's real numbering, so neither mistake is available.
 *
 * Each body carries its own verdict: for every copy, how many CA atoms matched
 * and what RMSD the superposition reached. That is the answer to "do these
 * residues really move as one unit", and it arrives before the job runs.
 */

interface BodiesEditorProps {
  bodies: ParsedBody[];
  assembly: Assembly;
  chains: ChainInfo[];
  preview?: PreviewBody[];
  disabled?: boolean;
  onChange: (next: ParsedBody[]) => void;
  onHover?: (index: number | null) => void;
}

/** A residue number that only reaches the parameter when it is a number and
 *  inside the chain. Typing "2" on the way to "298" must not write a body
 *  that covers one residue. */
const ResidueField: React.FC<{
  value: number;
  bounds: { lo: number; hi: number } | null;
  disabled?: boolean;
  label: string;
  onCommit: (value: number) => void;
}> = ({ value, bounds, disabled, label, onCommit }) => {
  const [draft, setDraft] = useState(String(value));
  useEffect(() => setDraft(String(value)), [value]);

  const commit = useCallback(() => {
    const parsed = Number(draft);
    if (!Number.isFinite(parsed) || !Number.isInteger(parsed)) {
      setDraft(String(value));
      return;
    }
    const clamped = bounds
      ? Math.min(Math.max(parsed, bounds.lo), bounds.hi)
      : parsed;
    setDraft(String(clamped));
    if (clamped !== value) onCommit(clamped);
  }, [draft, value, bounds, onCommit]);

  return (
    <TextField
      size="small"
      variant="standard"
      type="number"
      value={draft}
      disabled={disabled}
      onChange={(e) => setDraft(e.target.value)}
      onBlur={commit}
      onKeyDown={(e) => {
        if (e.key === "Enter") (e.target as HTMLInputElement).blur();
      }}
      sx={{ width: "5rem" }}
      // aria-label has to reach the input itself; on the TextField it lands on
      // the wrapper, where nothing that reads the form will find it.
      slotProps={{
        htmlInput: { min: bounds?.lo, max: bounds?.hi, "aria-label": label },
      }}
    />
  );
};

const CopyVerdict: React.FC<{ body?: PreviewBody }> = ({ body }) => {
  if (!body) return null;
  if (body.error) {
    return (
      <Typography variant="caption" color="error.main">
        {body.error}
      </Typography>
    );
  }
  if (body.mode === "exclude") {
    return (
      <Typography variant="caption" color="text.secondary">
        excluded — no operators are fitted for this body
      </Typography>
    );
  }
  if (!body.copies.length) return null;
  return (
    <Box sx={{ display: "flex", gap: 0.75, flexWrap: "wrap", alignItems: "center" }}>
      {body.copies.map((copy) => {
        if (copy.skipped) {
          return (
            <Chip
              key={copy.label}
              size="small"
              variant="outlined"
              label={`${copy.label}: not in this copy`}
            />
          );
        }
        if (copy.error || copy.rmsd === undefined) {
          return (
            <Chip
              key={copy.label}
              size="small"
              color="error"
              variant="outlined"
              label={`${copy.label}: ${copy.error ?? "no fit"}`}
            />
          );
        }
        return (
          <Tooltip
            key={copy.label}
            title={`${copy.nCA} CA atoms matched between the reference and copy ${copy.label}`}
          >
            <Chip
              size="small"
              color={rmsdTone(copy.rmsd)}
              variant="outlined"
              label={`${copy.label}: ${copy.nCA} CA, ${copy.rmsd.toFixed(2)} Å`}
            />
          </Tooltip>
        );
      })}
    </Box>
  );
};

export const BodiesEditor: React.FC<BodiesEditorProps> = ({
  bodies,
  assembly,
  chains,
  preview,
  disabled,
  onChange,
  onHover,
}) => {
  const theme = useTheme();
  const roles = assembly.roles.length ? assembly.roles : [IMPLICIT_ROLE];
  const showRole = roles.length > 1;

  const patch = useCallback(
    (index: number, change: Partial<ParsedBody>) =>
      onChange(bodies.map((body, i) => (i === index ? { ...body, ...change } : body))),
    [bodies, onChange]
  );

  const patchSegment = useCallback(
    (bodyIndex: number, segmentIndex: number, change: Partial<Segment>) =>
      patch(bodyIndex, {
        segments: bodies[bodyIndex].segments.map((segment, i) =>
          i === segmentIndex ? { ...segment, ...change } : segment
        ),
      }),
    [bodies, patch]
  );

  const newSegment = useCallback(
    (role: string): Segment => {
      const bounds = referenceBounds(assembly, role, chains);
      return { role, lo: bounds?.lo ?? 1, hi: bounds?.hi ?? 1 };
    },
    [assembly, chains]
  );

  return (
    <Box>
      <Typography variant="body2" color="text.secondary" sx={{ mb: 1 }}>
        Each body is superposed across the copies on its own, so two bodies can
        follow different NCS operators — the reason to use this task rather than
        plain NCS averaging. Split the copy where it hinges.
      </Typography>

      {bodies.map((body, index) => {
        const verdict = preview?.[index];
        return (
          <Paper
            key={index}
            variant="outlined"
            sx={{
              p: 1,
              mb: 1,
              borderLeft: `4px solid ${bodyColour(theme, index)}`,
              opacity: body.mode === "exclude" ? 0.75 : 1,
            }}
            onMouseEnter={() => onHover?.(index)}
            onMouseLeave={() => onHover?.(null)}
          >
            <Box sx={{ display: "flex", alignItems: "center", gap: 1, flexWrap: "wrap" }}>
              <Typography variant="body2" sx={{ fontWeight: 600, minWidth: "4rem" }}>
                body {index + 1}
              </Typography>

              <ToggleButtonGroup
                size="small"
                exclusive
                value={body.mode}
                disabled={disabled}
                onChange={(_, mode: Mode | null) => mode && patch(index, { mode })}
              >
                {MODES.map((mode) => (
                  <Tooltip key={mode} title={MODE_HELP[mode]} describeChild>
                    <ToggleButton value={mode} sx={{ textTransform: "none", px: 1.25 }}>
                      {mode}
                    </ToggleButton>
                  </Tooltip>
                ))}
              </ToggleButtonGroup>

              <Box sx={{ flex: 1 }} />
              <Tooltip title="Remove this rigid body">
                <span>
                  <IconButton
                    size="small"
                    color="error"
                    disabled={disabled}
                    aria-label={`Remove body ${index + 1}`}
                    onClick={() => onChange(bodies.filter((_, i) => i !== index))}
                  >
                    <Delete fontSize="small" />
                  </IconButton>
                </span>
              </Tooltip>
            </Box>

            {body.error ? (
              <Alert severity="warning" sx={{ mt: 1 }}>
                This body is stored as <code>{body.raw}</code>, which cannot be
                read: {body.error}. Delete it and add the ranges again.
              </Alert>
            ) : (
              <Box sx={{ display: "flex", flexWrap: "wrap", gap: 1, alignItems: "center", mt: 0.5 }}>
                {body.segments.map((segment, segmentIndex) => {
                  const bounds = referenceBounds(assembly, segment.role, chains);
                  return (
                    <Box
                      key={segmentIndex}
                      sx={{
                        display: "flex",
                        alignItems: "center",
                        gap: 0.5,
                        border: 1,
                        borderColor: "divider",
                        borderRadius: 1,
                        px: 1,
                        py: 0.25,
                      }}
                    >
                      {showRole && (
                        <Select
                          size="small"
                          variant="standard"
                          value={roles.includes(segment.role) ? segment.role : ""}
                          disabled={disabled}
                          onChange={(e) =>
                            patchSegment(index, segmentIndex, {
                              role: String(e.target.value),
                            })
                          }
                          sx={{ minWidth: "6rem" }}
                          displayEmpty
                          renderValue={(role) =>
                            role ? (
                              roleLabel(String(role), assembly)
                            ) : (
                              <Typography component="span" variant="body2" color="error.main">
                                {segment.role}?
                              </Typography>
                            )
                          }
                        >
                          {roles.map((role) => (
                            <MenuItem key={role} value={role}>
                              {roleLabel(role, assembly)}
                            </MenuItem>
                          ))}
                        </Select>
                      )}
                      <ResidueField
                        value={segment.lo}
                        bounds={bounds}
                        disabled={disabled}
                        label={`body ${index + 1} segment ${segmentIndex + 1} first residue`}
                        onCommit={(lo) =>
                          patchSegment(index, segmentIndex, {
                            lo,
                            hi: Math.max(lo, segment.hi),
                          })
                        }
                      />
                      <Typography variant="body2">–</Typography>
                      <ResidueField
                        value={segment.hi}
                        bounds={bounds}
                        disabled={disabled}
                        label={`body ${index + 1} segment ${segmentIndex + 1} last residue`}
                        onCommit={(hi) =>
                          patchSegment(index, segmentIndex, {
                            hi,
                            lo: Math.min(hi, segment.lo),
                          })
                        }
                      />
                      {bounds && (
                        <Typography variant="caption" color="text.secondary">
                          of {bounds.lo}–{bounds.hi}
                        </Typography>
                      )}
                      <IconButton
                        size="small"
                        disabled={disabled || body.segments.length < 2}
                        aria-label="Remove segment"
                        onClick={() =>
                          patch(index, {
                            segments: body.segments.filter((_, i) => i !== segmentIndex),
                          })
                        }
                      >
                        <Delete fontSize="inherit" />
                      </IconButton>
                    </Box>
                  );
                })}

                <Tooltip
                  title={
                    showRole
                      ? "Add a range — from another entity, to make a body that crosses chains"
                      : "Add another range to this body"
                  }
                >
                  <span>
                    <Button
                      size="small"
                      startIcon={<Add />}
                      disabled={disabled}
                      onClick={() =>
                        patch(index, {
                          segments: [...body.segments, newSegment(roles[0])],
                        })
                      }
                    >
                      range
                    </Button>
                  </span>
                </Tooltip>
              </Box>
            )}

            <Box sx={{ mt: 0.75 }}>
              <CopyVerdict body={verdict} />
            </Box>
          </Paper>
        );
      })}

      <Button
        size="small"
        startIcon={<Add />}
        disabled={disabled}
        onClick={() =>
          onChange([...bodies, { segments: [newSegment(roles[0])], mode: "average" }])
        }
      >
        Add rigid body
      </Button>
    </Box>
  );
};

export default BodiesEditor;
