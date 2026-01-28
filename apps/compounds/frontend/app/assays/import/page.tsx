'use client';

import { useState, useMemo } from 'react';
import { useRouter, useSearchParams } from 'next/navigation';
import Link from 'next/link';
import {
  Container,
  Paper,
  Typography,
  Box,
  Button,
  Alert,
  List,
  ListItemButton,
  ListItemIcon,
  ListItemText,
  Chip,
  Divider,
  Skeleton,
  TextField,
  InputAdornment,
} from '@mui/material';
import {
  ArrowBack,
  Description,
  GridOn,
  Search,
  Science,
  Warning,
} from '@mui/icons-material';
import { PageHeader } from '@/components/compounds/PageHeader';
import { useCompoundsApi } from '@/lib/compounds/api';
import { routes } from '@/lib/compounds/routes';
import type { Protocol, ImportType } from '@/types/compounds/models';

// Helper to check if protocol is plate-based (not ADME or table of values)
function isPlateBasedProtocol(importType: ImportType): boolean {
  return importType === 'raw_data' || importType === 'ms_intact';
}

const IMPORT_TYPE_LABELS: Record<ImportType, string> = {
  raw_data: 'Raw Data (Dose-Response)',
  ms_intact: 'MS-Intact',
  table_of_values: 'Table of Values',
  pharmaron_adme: 'Pharmaron ADME',
};

export default function ImportAssayPage() {
  const router = useRouter();
  const searchParams = useSearchParams();
  const api = useCompoundsApi();

  // Get target from URL query params (passed from AddAssayMenu)
  const targetId = searchParams.get('target');

  const [searchQuery, setSearchQuery] = useState('');

  // Fetch all protocols
  const { data: protocolsData, isLoading } = api.get<Protocol[]>('protocols/');

  // Filter to plate-based protocols only
  const plateBasedProtocols = useMemo(() => {
    if (!protocolsData) return [];
    return protocolsData.filter(p => isPlateBasedProtocol(p.import_type));
  }, [protocolsData]);

  // Further filter by search query
  const filteredProtocols = useMemo(() => {
    if (!searchQuery.trim()) return plateBasedProtocols;
    const query = searchQuery.toLowerCase();
    return plateBasedProtocols.filter(p =>
      p.name.toLowerCase().includes(query)
    );
  }, [plateBasedProtocols, searchQuery]);

  // Handle protocol selection - navigate to protocol page with openUpload flag
  const handleProtocolSelect = (protocol: Protocol) => {
    const params = new URLSearchParams();
    params.set('openUpload', 'true');
    if (targetId) {
      params.set('target', targetId);
    }
    router.push(`${routes.assays.protocol(protocol.id)}?${params.toString()}`);
  };

  return (
    <Container maxWidth="md" sx={{ py: 4 }}>
      <PageHeader
        breadcrumbs={[
          { label: 'Assays', href: routes.assays.list() },
          { label: 'Import Plate Data' },
        ]}
      />

      <Box sx={{ display: 'flex', alignItems: 'center', gap: 2, mb: 3 }}>
        <Button
          component={Link}
          href={routes.assays.list()}
          startIcon={<ArrowBack />}
          size="small"
        >
          Back to Assays
        </Button>
        <Typography variant="h4" sx={{ flex: 1 }}>
          Import Plate Data
        </Typography>
      </Box>

      <Alert severity="info" sx={{ mb: 3 }}>
        <Typography variant="body2">
          Select a protocol to import plate reader data. Each protocol has a configured plate layout
          that defines how data is extracted from your Excel files.
        </Typography>
      </Alert>

      <Paper sx={{ p: 3 }}>
        <Box sx={{ display: 'flex', alignItems: 'center', gap: 2, mb: 2 }}>
          <Science color="primary" />
          <Typography variant="h6">
            Select Protocol
          </Typography>
          <Chip
            label={`${plateBasedProtocols.length} plate-based protocols`}
            size="small"
            variant="outlined"
          />
        </Box>

        <Typography variant="body2" color="text.secondary" sx={{ mb: 3 }}>
          Choose the assay protocol that matches your data. The protocol&apos;s plate layout
          configuration will be used to extract data series from your spreadsheet.
        </Typography>

        {/* Search field */}
        <TextField
          placeholder="Search protocols..."
          value={searchQuery}
          onChange={(e) => setSearchQuery(e.target.value)}
          size="small"
          fullWidth
          sx={{ mb: 2 }}
          slotProps={{
            input: {
              startAdornment: (
                <InputAdornment position="start">
                  <Search color="action" />
                </InputAdornment>
              ),
            },
          }}
        />

        <Divider sx={{ mb: 2 }} />

        {isLoading ? (
          <Box>
            {[1, 2, 3, 4].map((i) => (
              <Skeleton key={i} variant="rectangular" height={72} sx={{ mb: 1, borderRadius: 1 }} />
            ))}
          </Box>
        ) : filteredProtocols.length === 0 ? (
          <Box sx={{ py: 4, textAlign: 'center' }}>
            {searchQuery ? (
              <>
                <Search sx={{ fontSize: 48, color: 'grey.400', mb: 1 }} />
                <Typography color="text.secondary">
                  No protocols matching &quot;{searchQuery}&quot;
                </Typography>
              </>
            ) : (
              <>
                <Description sx={{ fontSize: 48, color: 'grey.400', mb: 1 }} />
                <Typography color="text.secondary">
                  No plate-based protocols found
                </Typography>
                <Typography variant="body2" color="text.secondary" sx={{ mt: 1 }}>
                  Create a protocol with &quot;Raw Data&quot; or &quot;MS-Intact&quot; import type first.
                </Typography>
              </>
            )}
          </Box>
        ) : (
          <List disablePadding>
            {filteredProtocols.map((protocol) => (
              <ListItemButton
                key={protocol.id}
                onClick={() => handleProtocolSelect(protocol)}
                sx={{
                  borderRadius: 1,
                  mb: 1,
                  border: '1px solid',
                  borderColor: 'divider',
                  '&:hover': {
                    borderColor: 'primary.main',
                    bgcolor: 'primary.50',
                  },
                }}
              >
                <ListItemIcon>
                  <Description color="primary" />
                </ListItemIcon>
                <ListItemText
                  primary={
                    <Box sx={{ display: 'flex', alignItems: 'center', gap: 1 }}>
                      <Typography fontWeight={500}>{protocol.name}</Typography>
                      <Chip
                        label={IMPORT_TYPE_LABELS[protocol.import_type]}
                        size="small"
                        variant="outlined"
                        color="primary"
                      />
                    </Box>
                  }
                  secondary={
                    <Box sx={{ display: 'flex', alignItems: 'center', gap: 2, mt: 0.5 }}>
                      {protocol.plate_layout_name ? (
                        <Chip
                          icon={<GridOn />}
                          label={protocol.plate_layout_name}
                          size="small"
                          variant="outlined"
                        />
                      ) : (
                        <Chip
                          icon={<Warning />}
                          label="No plate layout"
                          size="small"
                          color="warning"
                          variant="outlined"
                        />
                      )}
                      {protocol.assays_count !== undefined && (
                        <Typography variant="caption" color="text.secondary">
                          {protocol.assays_count} assay{protocol.assays_count !== 1 ? 's' : ''}
                        </Typography>
                      )}
                    </Box>
                  }
                />
              </ListItemButton>
            ))}
          </List>
        )}
      </Paper>

      {/* Alternative options */}
      <Paper sx={{ p: 3, mt: 3 }}>
        <Typography variant="subtitle2" gutterBottom>
          Other Import Options
        </Typography>
        <Typography variant="body2" color="text.secondary" sx={{ mb: 2 }}>
          For pre-analyzed data or vendor-specific formats, use these dedicated importers:
        </Typography>
        <Box sx={{ display: 'flex', gap: 2, flexWrap: 'wrap' }}>
          <Button
            component={Link}
            href={targetId ? `${routes.assays.importTableOfValues()}?target=${targetId}` : routes.assays.importTableOfValues()}
            variant="outlined"
            size="small"
          >
            Import Table of Values
          </Button>
          <Button
            component={Link}
            href={targetId ? `${routes.assays.importAdme()}?target=${targetId}` : routes.assays.importAdme()}
            variant="outlined"
            size="small"
          >
            Import Pharmaron ADME
          </Button>
        </Box>
      </Paper>
    </Container>
  );
}
