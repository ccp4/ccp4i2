'use client';

import { useState, useEffect, useCallback, useRef, useMemo } from 'react';
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
  Button,
} from '@mui/material';
import {
  Science,
  Description,
  Search,
  PlayArrow,
} from '@mui/icons-material';
import {
  Predicates,
  AggregationType,
  OutputFormat,
  ProtocolInfo,
} from '@/types/compounds/aggregation';
import { Target } from '@/types/compounds/models';
import { fetchTargets, fetchProtocols } from '@/lib/compounds/aggregation-api';

/** State exposed for URL sharing and saving */
export interface PredicateBuilderState {
  targets: Target[];
  targetNames: string[];
  protocolNames: string[];
  compoundSearch: string;
  outputFormat: OutputFormat;
  aggregations: AggregationType[];
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
  /** Called when user clicks Run Query */
  onRunQuery: (
    predicates: Predicates,
    outputFormat: OutputFormat,
    aggregations: AggregationType[]
  ) => void;
  /** Called when state changes (for URL sharing) */
  onStateChange?: (state: PredicateBuilderState) => void;
  /** Whether a query is currently running */
  loading?: boolean;
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
  onRunQuery,
  onStateChange,
  loading = false,
}: PredicateBuilderProps) {
  // Target selection
  const [selectedTargets, setSelectedTargets] = useState<Target[]>([]);
  const [targetOptions, setTargetOptions] = useState<Target[]>([]);
  const [targetLoading, setTargetLoading] = useState(false);

  // Protocol selection
  const [selectedProtocols, setSelectedProtocols] = useState<ProtocolInfo[]>([]);
  const [protocolOptions, setProtocolOptions] = useState<ProtocolInfo[]>([]);
  const [protocolLoading, setProtocolLoading] = useState(false);

  // Other predicates
  const [compoundSearch, setCompoundSearch] = useState(initialCompoundSearch || '');
  const [status, setStatus] = useState<string>('valid');

  // Output options
  const [outputFormat, setOutputFormat] = useState<OutputFormat>(initialOutputFormat || 'compact');
  const [aggregations, setAggregations] = useState<AggregationType[]>(['geomean', 'count']);

  // Track initialization state
  const [isInitialized, setIsInitialized] = useState(false);
  const hasRunInitialQuery = useRef(false);
  const hasStartedInit = useRef(false);

  // Memoize initial values to prevent infinite loops from array reference changes
  const memoizedInitialTargetId = useMemo(() => initialTargetId, [initialTargetId]);
  const memoizedInitialTargetNames = useMemo(
    () => initialTargetNames,
    // eslint-disable-next-line react-hooks/exhaustive-deps
    [initialTargetNames?.join(',')]
  );
  const memoizedInitialProtocolNames = useMemo(
    () => initialProtocolNames,
    // eslint-disable-next-line react-hooks/exhaustive-deps
    [initialProtocolNames?.join(',')]
  );

  // Check if we have URL params that should trigger auto-query
  const hasUrlParams = !!(initialTargetId || initialTargetNames?.length || initialProtocolNames?.length || initialCompoundSearch);

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

  // Check if query can be run
  const canRunQuery = hasPredicates() &&
    (outputFormat === 'long' || aggregations.length > 0);

  // Handle Run Query button click
  const handleRunQuery = useCallback(() => {
    if (canRunQuery) {
      onRunQuery(buildPredicates(), outputFormat, aggregations);
    }
  }, [canRunQuery, onRunQuery, buildPredicates, outputFormat, aggregations]);

  // Load initial targets and protocols if URL params provided
  useEffect(() => {
    // Prevent re-running initialization
    if (hasStartedInit.current) return;
    hasStartedInit.current = true;

    const initializeSelections = async () => {
      const promises: Promise<void>[] = [];

      // Load targets by ID (legacy) or by names
      if (memoizedInitialTargetId || (memoizedInitialTargetNames && memoizedInitialTargetNames.length > 0)) {
        promises.push(
          fetchTargets().then((targets) => {
            setTargetOptions(targets); // Also populate options
            if (memoizedInitialTargetId) {
              const target = targets.find((t) => t.id === memoizedInitialTargetId);
              if (target) {
                setSelectedTargets([target]);
              }
            } else if (memoizedInitialTargetNames) {
              const matchedTargets = targets.filter((t) =>
                memoizedInitialTargetNames.some((name) => t.name.toLowerCase() === name.toLowerCase())
              );
              if (matchedTargets.length > 0) {
                setSelectedTargets(matchedTargets);
              }
            }
          })
        );
      }

      // Load protocols by names
      if (memoizedInitialProtocolNames && memoizedInitialProtocolNames.length > 0) {
        promises.push(
          fetchProtocols({}).then((protocols) => {
            setProtocolOptions(protocols); // Also populate options
            const matchedProtocols = protocols.filter((p) =>
              memoizedInitialProtocolNames.some((name) => p.name.toLowerCase() === name.toLowerCase())
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
  }, [memoizedInitialTargetId, memoizedInitialTargetNames, memoizedInitialProtocolNames]);

  // Auto-run query once if URL params were provided
  useEffect(() => {
    if (isInitialized && hasUrlParams && !hasRunInitialQuery.current && hasPredicates()) {
      hasRunInitialQuery.current = true;
      const predicates = buildPredicates();
      // Use aggregations for compact/medium, ignore for long
      onRunQuery(predicates, outputFormat, aggregations);
    }
  }, [isInitialized, hasUrlParams, hasPredicates, buildPredicates, outputFormat, aggregations, onRunQuery]);

  // Load target options when user types in autocomplete
  const handleTargetSearch = useCallback((search: string) => {
    setTargetLoading(true);
    fetchTargets({ search })
      .then(setTargetOptions)
      .finally(() => setTargetLoading(false));
  }, []);

  // Load protocol options when user types in autocomplete
  const handleProtocolSearch = useCallback((search: string) => {
    const targetId = selectedTargets.length === 1 ? selectedTargets[0].id : undefined;
    setProtocolLoading(true);
    fetchProtocols({ target: targetId, search })
      .then((protocols) => setProtocolOptions(protocols))
      .finally(() => setProtocolLoading(false));
  }, [selectedTargets]);

  // Load options when autocomplete opens (if not already loaded)
  const handleTargetOpen = useCallback(() => {
    if (targetOptions.length === 0) {
      handleTargetSearch('');
    }
  }, [targetOptions.length, handleTargetSearch]);

  const handleProtocolOpen = useCallback(() => {
    if (protocolOptions.length === 0) {
      handleProtocolSearch('');
    }
  }, [protocolOptions.length, handleProtocolSearch]);

  // Notify parent of state changes for URL sharing
  useEffect(() => {
    if (onStateChange) {
      onStateChange({
        targets: selectedTargets,
        targetNames: selectedTargets.map((t) => t.name),
        protocolNames: selectedProtocols.map((p) => p.name),
        compoundSearch,
        outputFormat,
        aggregations,
      });
    }
  }, [selectedTargets, selectedProtocols, compoundSearch, outputFormat, aggregations, onStateChange]);

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
          onInputChange={(_, value, reason) => {
            if (reason === 'input') handleTargetSearch(value);
          }}
          onOpen={handleTargetOpen}
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
          onInputChange={(_, value, reason) => {
            if (reason === 'input') handleProtocolSearch(value);
          }}
          onOpen={handleProtocolOpen}
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

      {/* Row 2: Output options and Run button */}
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

        {/* Run Query button */}
        <Box sx={{ display: 'flex', gap: 1, alignItems: 'center', ml: 'auto' }}>
          <Button
            variant="contained"
            size="small"
            startIcon={loading ? <CircularProgress size={16} color="inherit" /> : <PlayArrow />}
            onClick={handleRunQuery}
            disabled={!canRunQuery || loading}
          >
            {loading ? 'Running...' : 'Run Query'}
          </Button>
        </Box>
      </Box>
    </Paper>
  );
}
