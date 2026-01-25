'use client';

import { useState, useEffect, useCallback, useRef } from 'react';
import {
  Box,
  Paper,
  Typography,
  TextField,
  Autocomplete,
  FormControl,
  FormControlLabel,
  Checkbox,
  CircularProgress,
  Chip,
  Select,
  MenuItem,
  InputLabel,
  ToggleButton,
  ToggleButtonGroup,
} from '@mui/material';
import {
  Science,
  Description,
  Search,
} from '@mui/icons-material';
import {
  Predicates,
  AggregationType,
  OutputFormat,
  ProtocolInfo,
} from '@/types/compounds/aggregation';
import { Target } from '@/types/compounds/models';
import { fetchTargets, fetchProtocols } from '@/lib/compounds/aggregation-api';

/** State exposed for URL sharing */
export interface PredicateBuilderState {
  targetNames: string[];
  protocolNames: string[];
  compoundSearch: string;
  outputFormat: OutputFormat;
}

interface PredicateBuilderProps {
  /** Initial target ID (for target detail page) - legacy support */
  initialTargetId?: string;
  /** Initial target names (for URL sharing) */
  initialTargetNames?: string[];
  /** Initial protocol names (for URL sharing) */
  initialProtocolNames?: string[];
  /** Initial compound search text (for compound detail navigation) */
  initialCompoundSearch?: string;
  /** Initial output format (for URL sharing) */
  initialOutputFormat?: OutputFormat;
  /** Called when predicates change (debounced for text fields) */
  onChange: (
    predicates: Predicates,
    outputFormat: OutputFormat,
    aggregations: AggregationType[]
  ) => void;
  /** Called when state changes (for URL sharing) */
  onStateChange?: (state: PredicateBuilderState) => void;
  /** Whether a query is currently running */
  loading?: boolean;
  /** Debounce delay in ms for text fields (default 500) */
  debounceMs?: number;
}

const AGGREGATION_OPTIONS: { value: AggregationType; label: string; description: string }[] = [
  { value: 'geomean', label: 'Geometric Mean', description: 'Average of log-transformed values' },
  { value: 'count', label: 'Count', description: 'Number of measurements' },
  { value: 'stdev', label: 'Std Dev', description: 'Sample standard deviation' },
  { value: 'list', label: 'List', description: 'All values as comma-separated list' },
];

export function PredicateBuilder({
  initialTargetId,
  initialTargetNames,
  initialProtocolNames,
  initialCompoundSearch,
  initialOutputFormat,
  onChange,
  onStateChange,
  loading = false,
  debounceMs = 500,
}: PredicateBuilderProps) {
  // Target selection
  const [selectedTargets, setSelectedTargets] = useState<Target[]>([]);
  const [targetOptions, setTargetOptions] = useState<Target[]>([]);
  const [targetLoading, setTargetLoading] = useState(false);
  const [targetSearch, setTargetSearch] = useState('');

  // Protocol selection
  const [selectedProtocols, setSelectedProtocols] = useState<ProtocolInfo[]>([]);
  const [protocolOptions, setProtocolOptions] = useState<ProtocolInfo[]>([]);
  const [protocolLoading, setProtocolLoading] = useState(false);
  const [protocolSearch, setProtocolSearch] = useState('');

  // Other predicates
  const [compoundSearch, setCompoundSearch] = useState(initialCompoundSearch || '');
  const [status, setStatus] = useState<string>('valid');

  // Output options
  const [outputFormat, setOutputFormat] = useState<OutputFormat>(initialOutputFormat || 'compact');
  const [aggregations, setAggregations] = useState<AggregationType[]>(['geomean', 'count']);

  // Debounce timer ref
  const debounceTimerRef = useRef<NodeJS.Timeout | null>(null);

  // Track if initial load is complete
  const [isInitialized, setIsInitialized] = useState(false);

  // Build predicates object
  const buildPredicates = useCallback((): Predicates => {
    const predicates: Predicates = {};
    if (selectedTargets.length > 0) {
      predicates.targets = selectedTargets.map((t) => t.id);
    }
    if (selectedProtocols.length > 0) {
      predicates.protocols = selectedProtocols.map((p) => p.id);
    }
    if (compoundSearch.trim()) {
      predicates.compound_search = compoundSearch.trim();
    }
    if (status) {
      predicates.status = status as any;
    }
    return predicates;
  }, [selectedTargets, selectedProtocols, compoundSearch, status]);

  // Check if we have any predicates
  const hasPredicates = useCallback(() => {
    return selectedTargets.length > 0 ||
      selectedProtocols.length > 0 ||
      compoundSearch.trim().length > 0;
  }, [selectedTargets, selectedProtocols, compoundSearch]);

  // Trigger onChange (debounced version)
  const triggerChangeDebounced = useCallback(() => {
    if (debounceTimerRef.current) {
      clearTimeout(debounceTimerRef.current);
    }
    debounceTimerRef.current = setTimeout(() => {
      if (hasPredicates() && aggregations.length > 0) {
        onChange(buildPredicates(), outputFormat, aggregations);
      }
    }, debounceMs);
  }, [buildPredicates, hasPredicates, onChange, outputFormat, aggregations, debounceMs]);

  // Trigger onChange immediately (for non-text fields)
  const triggerChangeImmediate = useCallback(() => {
    if (debounceTimerRef.current) {
      clearTimeout(debounceTimerRef.current);
    }
    // Small delay to allow state to update
    setTimeout(() => {
      if (hasPredicates() && aggregations.length > 0) {
        onChange(buildPredicates(), outputFormat, aggregations);
      }
    }, 50);
  }, [buildPredicates, hasPredicates, onChange, outputFormat, aggregations]);

  // Load initial targets and protocols if provided
  useEffect(() => {
    const initializeSelections = async () => {
      const promises: Promise<void>[] = [];

      // Load targets by ID (legacy) or by names
      if (initialTargetId || (initialTargetNames && initialTargetNames.length > 0)) {
        promises.push(
          fetchTargets().then((targets) => {
            if (initialTargetId) {
              const target = targets.find((t) => t.id === initialTargetId);
              if (target) {
                setSelectedTargets([target]);
              }
            } else if (initialTargetNames) {
              const matchedTargets = targets.filter((t) =>
                initialTargetNames.some((name) => t.name.toLowerCase() === name.toLowerCase())
              );
              if (matchedTargets.length > 0) {
                setSelectedTargets(matchedTargets);
              }
            }
          })
        );
      }

      // Load protocols by names
      if (initialProtocolNames && initialProtocolNames.length > 0) {
        promises.push(
          fetchProtocols({}).then((protocols) => {
            const matchedProtocols = protocols.filter((p) =>
              initialProtocolNames.some((name) => p.name.toLowerCase() === name.toLowerCase())
            );
            if (matchedProtocols.length > 0) {
              setSelectedProtocols(matchedProtocols);
            }
          })
        );
      }

      await Promise.all(promises);
      setIsInitialized(true);
    };

    initializeSelections();
  }, [initialTargetId, initialTargetNames, initialProtocolNames]);

  // Load target options on search
  useEffect(() => {
    setTargetLoading(true);
    fetchTargets({ search: targetSearch })
      .then(setTargetOptions)
      .finally(() => setTargetLoading(false));
  }, [targetSearch]);

  // Load protocol options on search (filtered by selected targets)
  useEffect(() => {
    setProtocolLoading(true);
    const targetId = selectedTargets.length === 1 ? selectedTargets[0].id : undefined;
    fetchProtocols({ target: targetId, search: protocolSearch })
      .then((protocols) => setProtocolOptions(protocols as any))
      .finally(() => setProtocolLoading(false));
  }, [protocolSearch, selectedTargets]);

  // Auto-submit when initial values are loaded
  useEffect(() => {
    if (isInitialized && hasPredicates() && aggregations.length > 0) {
      onChange(buildPredicates(), outputFormat, aggregations);
    }
    // Only run once after initialization
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [isInitialized]);

  // Notify parent of state changes for URL sharing
  useEffect(() => {
    if (isInitialized && onStateChange) {
      onStateChange({
        targetNames: selectedTargets.map((t) => t.name),
        protocolNames: selectedProtocols.map((p) => p.name),
        compoundSearch,
        outputFormat,
      });
    }
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [selectedTargets, selectedProtocols, compoundSearch, outputFormat, isInitialized]);

  // Trigger on output format or aggregation changes (immediate)
  useEffect(() => {
    if (isInitialized && hasPredicates() && aggregations.length > 0) {
      triggerChangeImmediate();
    }
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [outputFormat, aggregations]);

  // Trigger on compound search changes (debounced)
  useEffect(() => {
    if (isInitialized) {
      triggerChangeDebounced();
    }
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [compoundSearch]);

  // Trigger on selection changes (immediate)
  useEffect(() => {
    if (isInitialized) {
      triggerChangeImmediate();
    }
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [selectedTargets, selectedProtocols, status]);

  // Cleanup debounce timer
  useEffect(() => {
    return () => {
      if (debounceTimerRef.current) {
        clearTimeout(debounceTimerRef.current);
      }
    };
  }, []);

  const handleAggregationChange = (agg: AggregationType) => {
    setAggregations((prev) =>
      prev.includes(agg) ? prev.filter((a) => a !== agg) : [...prev, agg]
    );
  };

  return (
    <Paper sx={{ p: 2, mb: 2 }}>
      <Typography variant="subtitle1" sx={{ mb: 1.5, fontWeight: 600 }}>
        Query Builder
      </Typography>

      {/* Row 1: Filters */}
      <Box sx={{ display: 'flex', gap: 2, flexWrap: 'wrap', mb: 1.5 }}>
        {/* Target selection */}
        <Autocomplete
          multiple
          options={targetOptions}
          value={selectedTargets}
          onChange={(_, newValue) => setSelectedTargets(newValue)}
          onInputChange={(_, value) => setTargetSearch(value)}
          getOptionLabel={(option) => option.name}
          isOptionEqualToValue={(option, value) => option.id === value.id}
          loading={targetLoading}
          size="small"
          sx={{ minWidth: 220, flex: 1 }}
          renderInput={(params) => (
            <TextField
              {...params}
              label="Targets"
              placeholder="Select..."
              InputProps={{
                ...params.InputProps,
                startAdornment: (
                  <>
                    <Science sx={{ mr: 0.5, color: 'action.active', fontSize: 18 }} />
                    {params.InputProps.startAdornment}
                  </>
                ),
                endAdornment: (
                  <>
                    {targetLoading && <CircularProgress color="inherit" size={16} />}
                    {params.InputProps.endAdornment}
                  </>
                ),
              }}
            />
          )}
          renderTags={(value, getTagProps) =>
            value.map((option, index) => (
              <Chip
                {...getTagProps({ index })}
                key={option.id}
                label={option.name}
                size="small"
              />
            ))
          }
        />

        {/* Protocol selection */}
        <Autocomplete
          multiple
          options={protocolOptions}
          value={selectedProtocols}
          onChange={(_, newValue) => setSelectedProtocols(newValue)}
          onInputChange={(_, value) => setProtocolSearch(value)}
          getOptionLabel={(option) => option.name}
          isOptionEqualToValue={(option, value) => option.id === value.id}
          loading={protocolLoading}
          size="small"
          sx={{ minWidth: 220, flex: 1 }}
          renderInput={(params) => (
            <TextField
              {...params}
              label="Protocols"
              placeholder="Select..."
              InputProps={{
                ...params.InputProps,
                startAdornment: (
                  <>
                    <Description sx={{ mr: 0.5, color: 'action.active', fontSize: 18 }} />
                    {params.InputProps.startAdornment}
                  </>
                ),
                endAdornment: (
                  <>
                    {protocolLoading && <CircularProgress color="inherit" size={16} />}
                    {params.InputProps.endAdornment}
                  </>
                ),
              }}
            />
          )}
          renderTags={(value, getTagProps) =>
            value.map((option, index) => (
              <Chip
                {...getTagProps({ index })}
                key={option.id}
                label={option.name}
                size="small"
              />
            ))
          }
        />

        {/* Compound search */}
        <TextField
          label="Compound Search"
          placeholder="NCL-00026..."
          value={compoundSearch}
          onChange={(e) => setCompoundSearch(e.target.value)}
          size="small"
          sx={{ minWidth: 180, flex: 0.8 }}
          InputProps={{
            startAdornment: <Search sx={{ mr: 0.5, color: 'action.active', fontSize: 18 }} />,
          }}
        />

        {/* Status filter */}
        <FormControl size="small" sx={{ minWidth: 130 }}>
          <InputLabel>Status</InputLabel>
          <Select
            value={status}
            label="Status"
            onChange={(e) => setStatus(e.target.value)}
          >
            <MenuItem value="valid">Valid</MenuItem>
            <MenuItem value="invalid">Invalid</MenuItem>
            <MenuItem value="unassigned">Unassigned</MenuItem>
            <MenuItem value="">All</MenuItem>
          </Select>
        </FormControl>
      </Box>

      {/* Row 2: Output options */}
      <Box sx={{ display: 'flex', gap: 2, flexWrap: 'wrap', alignItems: 'center' }}>
        {/* Output format toggle */}
        <Box sx={{ display: 'flex', alignItems: 'center', gap: 1 }}>
          <Typography variant="body2" color="text.secondary">
            Format:
          </Typography>
          <ToggleButtonGroup
            value={outputFormat}
            exclusive
            onChange={(_, value) => value && setOutputFormat(value)}
            size="small"
          >
            <ToggleButton value="compact" sx={{ px: 1.5, py: 0.5 }}>
              Compact
            </ToggleButton>
            <ToggleButton value="medium" sx={{ px: 1.5, py: 0.5 }}>
              Medium
            </ToggleButton>
            <ToggleButton value="long" sx={{ px: 1.5, py: 0.5 }}>
              Long
            </ToggleButton>
          </ToggleButtonGroup>
        </Box>

        {/* Aggregation functions - shown for compact and medium formats */}
        {(outputFormat === 'compact' || outputFormat === 'medium') && (
          <Box sx={{ display: 'flex', alignItems: 'center', gap: 0.5, ml: 1 }}>
            <Typography variant="body2" color="text.secondary">
              Aggregations:
            </Typography>
            {AGGREGATION_OPTIONS.map((opt) => (
              <FormControlLabel
                key={opt.value}
                control={
                  <Checkbox
                    checked={aggregations.includes(opt.value)}
                    onChange={() => handleAggregationChange(opt.value)}
                    size="small"
                    sx={{ py: 0 }}
                  />
                }
                label={<Typography variant="body2">{opt.label}</Typography>}
                title={opt.description}
                sx={{ mr: 1, ml: 0 }}
              />
            ))}
          </Box>
        )}

        {/* Status indicator */}
        <Box sx={{ display: 'flex', gap: 1, alignItems: 'center', ml: 'auto' }}>
          {loading && (
            <Box sx={{ display: 'flex', alignItems: 'center', gap: 0.5 }}>
              <CircularProgress size={16} />
              <Typography variant="body2" color="text.secondary">
                Loading...
              </Typography>
            </Box>
          )}
          {!loading && !hasPredicates() && (
            <Typography variant="body2" color="text.secondary">
              Select a filter to see results
            </Typography>
          )}
          {!loading && hasPredicates() && (outputFormat === 'compact' || outputFormat === 'medium') && aggregations.length === 0 && (
            <Typography variant="body2" color="warning.main">
              Select an aggregation
            </Typography>
          )}
        </Box>
      </Box>
    </Paper>
  );
}
