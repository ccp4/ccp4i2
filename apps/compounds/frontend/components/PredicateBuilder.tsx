'use client';

import { useState, useEffect, useCallback, useRef } from 'react';
import {
  Box,
  Paper,
  Typography,
  TextField,
  Autocomplete,
  FormControl,
  FormLabel,
  FormGroup,
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
} from '@/types/aggregation';
import { Target, Protocol } from '@/types/models';
import { fetchTargets, fetchProtocols } from '@/lib/aggregation-api';

interface PredicateBuilderProps {
  /** Initial target ID (for target detail page) */
  initialTargetId?: string;
  /** Initial compound search text (for compound detail navigation) */
  initialCompoundSearch?: string;
  /** Called when predicates change (debounced for text fields) */
  onChange: (
    predicates: Predicates,
    outputFormat: OutputFormat,
    aggregations: AggregationType[]
  ) => void;
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
  initialCompoundSearch,
  onChange,
  loading = false,
  debounceMs = 500,
}: PredicateBuilderProps) {
  // Target selection
  const [selectedTargets, setSelectedTargets] = useState<Target[]>([]);
  const [targetOptions, setTargetOptions] = useState<Target[]>([]);
  const [targetLoading, setTargetLoading] = useState(false);
  const [targetSearch, setTargetSearch] = useState('');

  // Protocol selection
  const [selectedProtocols, setSelectedProtocols] = useState<Protocol[]>([]);
  const [protocolOptions, setProtocolOptions] = useState<Protocol[]>([]);
  const [protocolLoading, setProtocolLoading] = useState(false);
  const [protocolSearch, setProtocolSearch] = useState('');

  // Other predicates
  const [compoundSearch, setCompoundSearch] = useState(initialCompoundSearch || '');
  const [status, setStatus] = useState<string>('valid');

  // Output options
  const [outputFormat, setOutputFormat] = useState<OutputFormat>('compact');
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

  // Load initial target if provided
  useEffect(() => {
    if (initialTargetId) {
      fetchTargets().then((targets) => {
        const target = targets.find((t) => t.id === initialTargetId);
        if (target) {
          setSelectedTargets([target]);
        }
        setIsInitialized(true);
      });
    } else {
      setIsInitialized(true);
    }
  }, [initialTargetId]);

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
    <Paper sx={{ p: 3, mb: 3 }}>
      <Typography variant="h6" gutterBottom>
        Query Builder
      </Typography>

      <Box sx={{ display: 'flex', flexDirection: 'column', gap: 3 }}>
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
          renderInput={(params) => (
            <TextField
              {...params}
              label="Targets"
              placeholder="Select targets..."
              InputProps={{
                ...params.InputProps,
                startAdornment: (
                  <>
                    <Science sx={{ mr: 1, color: 'action.active' }} />
                    {params.InputProps.startAdornment}
                  </>
                ),
                endAdornment: (
                  <>
                    {targetLoading && <CircularProgress color="inherit" size={20} />}
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
                icon={<Science fontSize="small" />}
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
          renderInput={(params) => (
            <TextField
              {...params}
              label="Protocols"
              placeholder="Select protocols..."
              InputProps={{
                ...params.InputProps,
                startAdornment: (
                  <>
                    <Description sx={{ mr: 1, color: 'action.active' }} />
                    {params.InputProps.startAdornment}
                  </>
                ),
                endAdornment: (
                  <>
                    {protocolLoading && <CircularProgress color="inherit" size={20} />}
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
                icon={<Description fontSize="small" />}
              />
            ))
          }
        />

        {/* Compound search */}
        <TextField
          label="Compound Search"
          placeholder="NCL-00026... or compound name"
          value={compoundSearch}
          onChange={(e) => setCompoundSearch(e.target.value)}
          InputProps={{
            startAdornment: <Search sx={{ mr: 1, color: 'action.active' }} />,
          }}
        />

        {/* Status filter */}
        <FormControl size="small" sx={{ minWidth: 150 }}>
          <InputLabel>Status Filter</InputLabel>
          <Select
            value={status}
            label="Status Filter"
            onChange={(e) => setStatus(e.target.value)}
          >
            <MenuItem value="valid">Valid only</MenuItem>
            <MenuItem value="invalid">Invalid only</MenuItem>
            <MenuItem value="unassigned">Unassigned only</MenuItem>
            <MenuItem value="">All statuses</MenuItem>
          </Select>
        </FormControl>

        {/* Output format toggle */}
        <Box>
          <Typography variant="subtitle2" gutterBottom>
            Output Format
          </Typography>
          <ToggleButtonGroup
            value={outputFormat}
            exclusive
            onChange={(_, value) => value && setOutputFormat(value)}
            size="small"
          >
            <ToggleButton value="compact">
              Compact (one row per compound)
            </ToggleButton>
            <ToggleButton value="medium">
              Medium (one row per protocol)
            </ToggleButton>
            <ToggleButton value="long">
              Long (one row per measurement)
            </ToggleButton>
          </ToggleButtonGroup>
        </Box>

        {/* Aggregation functions - shown for compact and medium formats */}
        {(outputFormat === 'compact' || outputFormat === 'medium') && (
          <FormControl component="fieldset">
            <FormLabel component="legend">Aggregation Functions</FormLabel>
            <FormGroup row>
              {AGGREGATION_OPTIONS.map((opt) => (
                <FormControlLabel
                  key={opt.value}
                  control={
                    <Checkbox
                      checked={aggregations.includes(opt.value)}
                      onChange={() => handleAggregationChange(opt.value)}
                    />
                  }
                  label={opt.label}
                  title={opt.description}
                />
              ))}
            </FormGroup>
          </FormControl>
        )}

        {/* Status indicator */}
        <Box sx={{ display: 'flex', gap: 2, alignItems: 'center', minHeight: 36 }}>
          {loading && (
            <Box sx={{ display: 'flex', alignItems: 'center', gap: 1 }}>
              <CircularProgress size={20} />
              <Typography variant="body2" color="text.secondary">
                Running query...
              </Typography>
            </Box>
          )}
          {!loading && !hasPredicates() && (
            <Typography variant="body2" color="text.secondary">
              Select at least one target, protocol, or enter a compound search to see results
            </Typography>
          )}
          {!loading && hasPredicates() && (outputFormat === 'compact' || outputFormat === 'medium') && aggregations.length === 0 && (
            <Typography variant="body2" color="warning.main">
              Select at least one aggregation function
            </Typography>
          )}
        </Box>
      </Box>
    </Paper>
  );
}
