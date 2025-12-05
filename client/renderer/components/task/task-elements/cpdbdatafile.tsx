import {
  Box,
  Checkbox,
  Chip,
  Collapse,
  Divider,
  FormControlLabel,
  FormGroup,
  IconButton,
  Stack,
  TextField,
  Tooltip,
  Typography,
} from "@mui/material";
import ExpandMoreIcon from "@mui/icons-material/ExpandMore";
import ExpandLessIcon from "@mui/icons-material/ExpandLess";
import { CCP4i2TaskElement, CCP4i2TaskElementProps } from "./task-element";
import { useCallback, useEffect, useMemo, useState } from "react";
import { useJob } from "../../../utils";
import { CSimpleDataFileElement } from "./csimpledatafile";

/**
 * Chain info from CPdbDataComposition
 * [nres, first_resid, last_resid]
 */
type ChainInfo = [number, string, string];

/**
 * Digest structure returned by the backend for CPdbDataFile
 */
interface CPdbDataFileDigest {
  composition?: {
    chains?: string[];
    chainInfo?: ChainInfo[];
    peptides?: string[];
    nucleics?: string[];
    solventChains?: string[];
    saccharides?: string[];
    monomers?: string[];
    nChains?: number;
    nResidues?: number;
    nAtoms?: number;
    nresSolvent?: number;
    elements?: string[];
    containsHydrogen?: boolean;
  };
  status?: string;
  reason?: string;
}

/**
 * Category options for quick selection toggles
 */
const CATEGORY_OPTIONS = [
  { key: "protein", label: "Protein", keyword: "protein" },
  { key: "nucleic", label: "Nucleic", keyword: "nucleic" },
  { key: "ligand", label: "Ligand", keyword: "ligand" },
  { key: "sugar", label: "Sugar", keyword: "sugar" },
  { key: "solvent", label: "Solvent", keyword: "solvent" },
] as const;

/**
 * Atom-level selection options
 */
const ATOM_OPTIONS = [
  { key: "backbone", label: "Backbone only", keyword: "backbone" },
  { key: "sidechain", label: "Sidechain only", keyword: "sidechain" },
] as const;

/**
 * Build a selection string from the UI state
 */
function buildSelectionString(
  selectedChains: Set<string>,
  allChains: string[],
  includedCategories: Set<string>,
  excludedCategories: Set<string>,
  atomFilter: string | null,
  residueRanges: Record<string, [string, string]>
): string {
  const parts: string[] = [];

  // Chain selections with residue ranges
  const chainParts: string[] = [];
  for (const chain of selectedChains) {
    const range = residueRanges[chain];
    if (range && (range[0] || range[1])) {
      const [start, end] = range;
      if (start && end) {
        chainParts.push(`${chain}/${start}-${end}`);
      } else if (start) {
        chainParts.push(`${chain}/${start}`);
      } else {
        chainParts.push(chain);
      }
    } else {
      chainParts.push(chain);
    }
  }

  // If not all chains selected, add chain selection
  if (chainParts.length > 0 && chainParts.length < allChains.length) {
    if (chainParts.length === 1) {
      parts.push(chainParts[0]);
    } else {
      parts.push(`{${chainParts.join(" or ")}}`);
    }
  } else if (chainParts.length > 0) {
    // Check if any chains have residue ranges
    const chainsWithRanges = chainParts.filter((p) => p.includes("/"));
    if (chainsWithRanges.length > 0) {
      if (chainsWithRanges.length === 1) {
        parts.push(chainsWithRanges[0]);
      } else {
        parts.push(`{${chainsWithRanges.join(" or ")}}`);
      }
    }
  }

  // Add included categories
  if (includedCategories.size > 0) {
    const catParts = Array.from(includedCategories);
    if (catParts.length === 1) {
      parts.push(catParts[0]);
    } else {
      parts.push(`{${catParts.join(" or ")}}`);
    }
  }

  // Build the main expression
  let expr = parts.length > 0 ? parts.join(" and ") : "";

  // Add atom filter if set
  if (atomFilter && expr) {
    expr = `${expr} and ${atomFilter}`;
  } else if (atomFilter) {
    expr = atomFilter;
  }

  // Add excluded categories
  if (excludedCategories.size > 0) {
    const excludeParts = Array.from(excludedCategories);
    const excludeExpr =
      excludeParts.length === 1
        ? excludeParts[0]
        : `{${excludeParts.join(" or ")}}`;
    if (expr) {
      expr = `${expr} and not ${excludeExpr}`;
    } else {
      expr = `not ${excludeExpr}`;
    }
  }

  return expr;
}

/**
 * CPdbDataFileElement - Renders a CPdbDataFile selector with selection builder.
 *
 * When a file is selected, the digest is fetched and provides:
 * - Chain selection checkboxes
 * - Residue range inputs per chain
 * - Category toggles (protein, nucleic, ligand, solvent, etc.)
 * - Atom-level filters (backbone, sidechain)
 * - Live preview of the generated selection string
 */
export const CPdbDataFileElement: React.FC<CCP4i2TaskElementProps> = (
  props
) => {
  const { job, itemName, qualifiers } = props;
  const { useTaskItem, useFileDigest, setParameter, mutateContainer } = useJob(
    job.id
  );
  const { item } = useTaskItem(itemName);

  // Selection text path
  const selectionItemName = useMemo(() => {
    return `${item?._objectPath}.selection.text`;
  }, [item]);

  const { value: selectionString, update: updateSelectionString } =
    useTaskItem(selectionItemName);

  // Get the file digest
  const { data: fileDigest, mutate: mutateDigest } = useFileDigest(
    item?._objectPath
  ) as { data: CPdbDataFileDigest | undefined; mutate: () => void };

  // UI state
  const [selectedChains, setSelectedChains] = useState<Set<string>>(new Set());
  const [includedCategories, setIncludedCategories] = useState<Set<string>>(
    new Set()
  );
  const [excludedCategories, setExcludedCategories] = useState<Set<string>>(
    new Set()
  );
  const [atomFilter, setAtomFilter] = useState<string | null>(null);
  const [residueRanges, setResidueRanges] = useState<
    Record<string, [string, string]>
  >({});
  const [isBuilderExpanded, setIsBuilderExpanded] = useState(false);
  const [manualEdit, setManualEdit] = useState(false);

  // Extract composition data
  const composition = fileDigest?.composition;
  const chains = composition?.chains || [];
  const chainInfo = composition?.chainInfo || [];
  const peptideChains = new Set(composition?.peptides || []);
  const nucleicChains = new Set(composition?.nucleics || []);
  const solventChains = new Set(composition?.solventChains || []);
  const hasContent = chains.length > 0;

  // Initialize selections when digest loads
  useEffect(() => {
    if (chains.length > 0 && selectedChains.size === 0) {
      setSelectedChains(new Set(chains));
    }
  }, [chains]);

  // Generate selection string from UI state
  const generatedSelection = useMemo(() => {
    if (!hasContent) return "";
    return buildSelectionString(
      selectedChains,
      chains,
      includedCategories,
      excludedCategories,
      atomFilter,
      residueRanges
    );
  }, [
    selectedChains,
    chains,
    includedCategories,
    excludedCategories,
    atomFilter,
    residueRanges,
    hasContent,
  ]);

  // Update backend when selection changes (debounced effect)
  useEffect(() => {
    if (manualEdit || !updateSelectionString || job.status !== 1) return;

    const timeoutId = setTimeout(async () => {
      if (generatedSelection !== selectionString) {
        try {
          await updateSelectionString(generatedSelection);
          await mutateContainer();
        } catch (error) {
          console.error("Error updating selection:", error);
        }
      }
    }, 300);

    return () => clearTimeout(timeoutId);
  }, [
    generatedSelection,
    selectionString,
    updateSelectionString,
    mutateContainer,
    manualEdit,
    job.status,
  ]);

  // Handle chain checkbox change
  const handleChainToggle = useCallback((chain: string, checked: boolean) => {
    setSelectedChains((prev) => {
      const next = new Set(prev);
      if (checked) {
        next.add(chain);
      } else {
        next.delete(chain);
      }
      return next;
    });
    setManualEdit(false);
  }, []);

  // Handle select all/none chains
  const handleSelectAllChains = useCallback(
    (selectAll: boolean) => {
      if (selectAll) {
        setSelectedChains(new Set(chains));
      } else {
        setSelectedChains(new Set());
      }
      setManualEdit(false);
    },
    [chains]
  );

  // Handle category include/exclude toggle
  const handleCategoryToggle = useCallback(
    (category: string, mode: "include" | "exclude" | "none") => {
      setIncludedCategories((prev) => {
        const next = new Set(prev);
        if (mode === "include") {
          next.add(category);
        } else {
          next.delete(category);
        }
        return next;
      });
      setExcludedCategories((prev) => {
        const next = new Set(prev);
        if (mode === "exclude") {
          next.add(category);
        } else {
          next.delete(category);
        }
        return next;
      });
      setManualEdit(false);
    },
    []
  );

  // Handle residue range change
  const handleRangeChange = useCallback(
    (chain: string, field: "start" | "end", value: string) => {
      setResidueRanges((prev) => {
        const current = prev[chain] || ["", ""];
        return {
          ...prev,
          [chain]: field === "start" ? [value, current[1]] : [current[0], value],
        };
      });
      setManualEdit(false);
    },
    []
  );

  // Handle manual text edit
  const handleManualTextChange = useCallback(
    async (value: string) => {
      setManualEdit(true);
      if (updateSelectionString && job.status === 1) {
        try {
          await updateSelectionString(value);
          await mutateContainer();
        } catch (error) {
          console.error("Error updating selection:", error);
        }
      }
    },
    [updateSelectionString, mutateContainer, job.status]
  );

  // Get chain type label
  const getChainTypeLabel = (chain: string): string => {
    const types: string[] = [];
    if (peptideChains.has(chain)) types.push("protein");
    if (nucleicChains.has(chain)) types.push("nucleic");
    if (solventChains.has(chain)) types.push("solvent");
    return types.join(", ") || "other";
  };

  // Get chain info for display
  const getChainDisplayInfo = (
    chain: string,
    index: number
  ): { nRes: number; firstRes: string; lastRes: string } => {
    const info = chainInfo[index];
    if (info) {
      return { nRes: info[0], firstRes: info[1], lastRes: info[2] };
    }
    return { nRes: 0, firstRes: "", lastRes: "" };
  };

  // Determine visibility
  const inferredVisibility = useMemo(() => {
    if (!props.visibility) return true;
    if (typeof props.visibility === "function") {
      return props.visibility();
    }
    return props.visibility;
  }, [props.visibility]);

  const overriddenQualifiers = useMemo(() => {
    return { ...item?._qualifiers, ...qualifiers };
  }, [item, qualifiers]);

  // Should show selection builder?
  const showSelectionBuilder = overriddenQualifiers.ifAtomSelection && hasContent;

  // Force expanded if selection is set
  const forceExpanded = !!selectionString && selectionString.length > 0;

  if (!inferredVisibility) return null;

  return (
    <CSimpleDataFileElement {...props} forceExpanded={forceExpanded}>
      {showSelectionBuilder && (
        <Stack spacing={1} sx={{ mt: 1 }}>
          {/* Selection Builder Header */}
          <Stack
            direction="row"
            alignItems="center"
            justifyContent="space-between"
          >
            <Typography variant="subtitle2" color="text.secondary">
              Atom Selection
            </Typography>
            <IconButton
              size="small"
              onClick={() => setIsBuilderExpanded(!isBuilderExpanded)}
            >
              {isBuilderExpanded ? <ExpandLessIcon /> : <ExpandMoreIcon />}
            </IconButton>
          </Stack>

          {/* Collapsed view: just the text field */}
          {!isBuilderExpanded && (
            <TextField
              size="small"
              fullWidth
              value={selectionString || ""}
              onChange={(e) => handleManualTextChange(e.target.value)}
              placeholder="e.g., protein and not solvent"
              disabled={job.status !== 1}
              InputProps={{
                sx: { fontFamily: "monospace", fontSize: "0.85rem" },
              }}
            />
          )}

          {/* Expanded view: full builder */}
          <Collapse in={isBuilderExpanded}>
            <Stack spacing={2} sx={{ p: 1, bgcolor: "action.hover", borderRadius: 1 }}>
              {/* Chain Selection */}
              <Box>
                <Stack direction="row" alignItems="center" spacing={1} sx={{ mb: 0.5 }}>
                  <Typography variant="caption" fontWeight="medium">
                    Chains
                  </Typography>
                  <Typography
                    variant="caption"
                    color="primary"
                    sx={{ cursor: "pointer" }}
                    onClick={() => handleSelectAllChains(true)}
                  >
                    All
                  </Typography>
                  <Typography variant="caption" color="text.secondary">
                    |
                  </Typography>
                  <Typography
                    variant="caption"
                    color="primary"
                    sx={{ cursor: "pointer" }}
                    onClick={() => handleSelectAllChains(false)}
                  >
                    None
                  </Typography>
                </Stack>
                <FormGroup row>
                  {chains.map((chain, idx) => {
                    const info = getChainDisplayInfo(chain, idx);
                    const typeLabel = getChainTypeLabel(chain);
                    return (
                      <Tooltip
                        key={chain}
                        title={`${info.nRes} residues (${info.firstRes}-${info.lastRes}) - ${typeLabel}`}
                      >
                        <FormControlLabel
                          control={
                            <Checkbox
                              checked={selectedChains.has(chain)}
                              onChange={(e) =>
                                handleChainToggle(chain, e.target.checked)
                              }
                              disabled={job.status !== 1}
                              size="small"
                            />
                          }
                          label={
                            <Stack direction="row" spacing={0.5} alignItems="center">
                              <Typography variant="body2">{chain}</Typography>
                              <Chip
                                label={typeLabel}
                                size="small"
                                variant="outlined"
                                sx={{ height: 16, fontSize: "0.65rem" }}
                              />
                            </Stack>
                          }
                          sx={{ mr: 2 }}
                        />
                      </Tooltip>
                    );
                  })}
                </FormGroup>
              </Box>

              {/* Residue Ranges */}
              {Array.from(selectedChains).length > 0 &&
                Array.from(selectedChains).length <= 4 && (
                  <Box>
                    <Typography variant="caption" fontWeight="medium" sx={{ mb: 0.5 }}>
                      Residue Ranges (optional)
                    </Typography>
                    <Stack direction="row" spacing={2} flexWrap="wrap">
                      {Array.from(selectedChains).map((chain) => {
                        const idx = chains.indexOf(chain);
                        const info = getChainDisplayInfo(chain, idx);
                        return (
                          <Stack key={chain} direction="row" spacing={0.5} alignItems="center">
                            <Typography variant="caption" sx={{ minWidth: 20 }}>
                              {chain}:
                            </Typography>
                            <TextField
                              size="small"
                              placeholder={info.firstRes || "start"}
                              value={residueRanges[chain]?.[0] || ""}
                              onChange={(e) =>
                                handleRangeChange(chain, "start", e.target.value)
                              }
                              disabled={job.status !== 1}
                              sx={{ width: 60 }}
                              inputProps={{ style: { fontSize: "0.75rem" } }}
                            />
                            <Typography variant="caption">-</Typography>
                            <TextField
                              size="small"
                              placeholder={info.lastRes || "end"}
                              value={residueRanges[chain]?.[1] || ""}
                              onChange={(e) =>
                                handleRangeChange(chain, "end", e.target.value)
                              }
                              disabled={job.status !== 1}
                              sx={{ width: 60 }}
                              inputProps={{ style: { fontSize: "0.75rem" } }}
                            />
                          </Stack>
                        );
                      })}
                    </Stack>
                  </Box>
                )}

              <Divider />

              {/* Category Filters */}
              <Box>
                <Typography variant="caption" fontWeight="medium" sx={{ mb: 0.5 }}>
                  Include / Exclude Categories
                </Typography>
                <Stack direction="row" spacing={1} flexWrap="wrap">
                  {CATEGORY_OPTIONS.map((cat) => {
                    const isIncluded = includedCategories.has(cat.keyword);
                    const isExcluded = excludedCategories.has(cat.keyword);
                    return (
                      <Chip
                        key={cat.key}
                        label={cat.label}
                        size="small"
                        variant={isIncluded || isExcluded ? "filled" : "outlined"}
                        color={isIncluded ? "primary" : isExcluded ? "error" : "default"}
                        onClick={() => {
                          if (isIncluded) {
                            handleCategoryToggle(cat.keyword, "exclude");
                          } else if (isExcluded) {
                            handleCategoryToggle(cat.keyword, "none");
                          } else {
                            handleCategoryToggle(cat.keyword, "include");
                          }
                        }}
                        disabled={job.status !== 1}
                        sx={{ cursor: "pointer" }}
                      />
                    );
                  })}
                </Stack>
                <Typography variant="caption" color="text.secondary" sx={{ mt: 0.5, display: "block" }}>
                  Click to cycle: none → include (blue) → exclude (red) → none
                </Typography>
              </Box>

              {/* Atom Level Filters */}
              <Box>
                <Typography variant="caption" fontWeight="medium" sx={{ mb: 0.5 }}>
                  Atom Filter
                </Typography>
                <Stack direction="row" spacing={1}>
                  <Chip
                    label="All atoms"
                    size="small"
                    variant={atomFilter === null ? "filled" : "outlined"}
                    color={atomFilter === null ? "primary" : "default"}
                    onClick={() => setAtomFilter(null)}
                    disabled={job.status !== 1}
                  />
                  {ATOM_OPTIONS.map((opt) => (
                    <Chip
                      key={opt.key}
                      label={opt.label}
                      size="small"
                      variant={atomFilter === opt.keyword ? "filled" : "outlined"}
                      color={atomFilter === opt.keyword ? "primary" : "default"}
                      onClick={() => setAtomFilter(opt.keyword)}
                      disabled={job.status !== 1}
                    />
                  ))}
                </Stack>
              </Box>

              <Divider />

              {/* Generated Selection String */}
              <Box>
                <Typography variant="caption" fontWeight="medium" sx={{ mb: 0.5 }}>
                  Selection String
                </Typography>
                <TextField
                  size="small"
                  fullWidth
                  multiline
                  minRows={1}
                  maxRows={3}
                  value={manualEdit ? selectionString : generatedSelection}
                  onChange={(e) => handleManualTextChange(e.target.value)}
                  placeholder="e.g., protein and not solvent"
                  disabled={job.status !== 1}
                  InputProps={{
                    sx: { fontFamily: "monospace", fontSize: "0.85rem" },
                  }}
                  helperText={
                    manualEdit
                      ? "Manually edited - UI controls disabled"
                      : "Auto-generated from selections above"
                  }
                />
              </Box>
            </Stack>
          </Collapse>
        </Stack>
      )}

      {/* Fallback: just show text input if no digest */}
      {overriddenQualifiers.ifAtomSelection && !hasContent && (
        <CCP4i2TaskElement
          {...props}
          itemName={selectionItemName}
          qualifiers={{
            ...useTaskItem(selectionItemName).item?._qualifiers,
            guiLabel: "Selection string",
          }}
        />
      )}
    </CSimpleDataFileElement>
  );
};
