'use client';

import { useState } from 'react';
import { useRouter } from 'next/navigation';
import Link from 'next/link';
import {
  Container,
  Typography,
  Box,
  Chip,
  Button,
  ToggleButtonGroup,
  ToggleButton,
  Tooltip,
  IconButton,
  Snackbar,
} from '@mui/material';
import {
  Science,
  Add,
  LocalShipping,
  Apps,
  ViewList,
  ViewModule,
  Biotech,
  FiberNew,
  DeleteForever,
  Hub,
} from '@mui/icons-material';
import { useSWRConfig } from 'swr';
import { PageHeader } from '@/components/compounds/PageHeader';
import { DataTable, Column } from '@/components/data-table';
import { TargetCardGrid } from '@/components/compounds/TargetCardGrid';
import { TargetCreateDialog } from '@/components/compounds/TargetCreateDialog';
import { TargetDeletionDialog } from '@/components/compounds/TargetDeletionDialog';
import { useCompoundsApi } from '@/lib/compounds/api';
import { routes } from '@/lib/compounds/routes';
import { Target } from '@/types/compounds/models';
import { useAuth } from '@/lib/compounds/auth-context';

type ViewMode = 'table' | 'grid';

export default function TargetsPage() {
  const router = useRouter();
  const { mutate } = useSWRConfig();
  const api = useCompoundsApi();
  const { canContribute, canAdminister } = useAuth();
  const { data: targets, isLoading } = api.get<Target[]>('targets/');
  const [createDialogOpen, setCreateDialogOpen] = useState(false);
  const [viewMode, setViewMode] = useState<ViewMode>('table');
  const [deleteTarget, setDeleteTarget] = useState<Target | null>(null);
  const [snackbarMessage, setSnackbarMessage] = useState<string | null>(null);

  const handleTargetCreated = () => {
    // Invalidate the targets cache to refresh the list
    mutate('/api/proxy/compounds/targets/');
  };

  const handleDeletionComplete = (summary: string) => {
    setSnackbarMessage(summary);
    mutate('/api/proxy/compounds/targets/');
    // Compound/template listings may include target names that are now stale.
    mutate(
      (key) =>
        typeof key === 'string' &&
        (key.includes('/compounds/') ||
          key.includes('/compound-templates/') ||
          key.includes('/hypotheses/')),
      undefined,
      { revalidate: true }
    );
  };

  const handleViewModeChange = (_: React.MouseEvent<HTMLElement>, newMode: ViewMode | null) => {
    if (newMode !== null) {
      setViewMode(newMode);
    }
  };

  const columns: Column<Target>[] = [
    {
      key: 'name',
      label: 'Target Name',
      sortable: true,
      searchable: true,
      render: (value, row) => (
        <Box sx={{ display: 'flex', alignItems: 'center', gap: 1 }}>
          <Science fontSize="small" color="primary" />
          <Typography fontWeight={500}>{value}</Typography>
          {(row.has_recent_compounds || row.has_recent_assays) && (
            <Tooltip title="New activity in the last 7 days">
              <FiberNew fontSize="small" color="secondary" />
            </Tooltip>
          )}
        </Box>
      ),
    },
    {
      key: 'compound_count',
      label: 'Compounds',
      sortable: true,
      width: 120,
      render: (value, row) =>
        value ? (
          <Chip
            label={value}
            size="small"
            color="primary"
            variant="outlined"
            onClick={(e) => {
              e.stopPropagation();
              router.push(routes.registry.compounds({ target: row.id }));
            }}
            sx={{ cursor: 'pointer' }}
          />
        ) : (
          '-'
        ),
    },
    {
      key: 'assay_count',
      label: 'Assays',
      sortable: true,
      width: 120,
      render: (value, row) =>
        value ? (
          <Chip
            icon={<Biotech fontSize="small" />}
            label={value}
            size="small"
            color="secondary"
            variant="outlined"
            onClick={(e) => {
              e.stopPropagation();
              router.push(routes.registry.target(row.id) + '?tab=assays');
            }}
            sx={{ cursor: 'pointer' }}
          />
        ) : (
          '-'
        ),
    },
    {
      key: 'created_at',
      label: 'Created',
      sortable: true,
      width: 120,
      hiddenOnMobile: true,
      render: (value) =>
        value ? new Date(value).toLocaleDateString() : '-',
    },
    {
      key: 'actions',
      label: '',
      width: 56,
      render: (_value, row) => (
        <Tooltip
          title={canAdminister ? 'Delete target' : 'Admin operating level required'}
          arrow
        >
          <span>
            <IconButton
              size="small"
              color="error"
              disabled={!canAdminister}
              onClick={(e) => {
                e.stopPropagation();
                setDeleteTarget(row);
              }}
            >
              <DeleteForever fontSize="small" />
            </IconButton>
          </span>
        </Tooltip>
      ),
    },
  ];

  return (
    <Box sx={{ display: 'flex', flexDirection: 'column', height: '100vh', overflow: 'hidden' }}>
      <Container maxWidth="lg" sx={{ flex: 1, display: 'flex', flexDirection: 'column', overflow: 'hidden', py: 2 }}>
        <Box sx={{ flexShrink: 0 }}>
          <PageHeader
            breadcrumbs={[
              { label: 'Home', href: routes.home(), icon: 'home' },
              { label: 'Registry', href: routes.registry.targets() },
              { label: 'Targets', icon: 'target' },
            ]}
          />

          <Box sx={{ mb: 2, display: 'flex', justifyContent: 'space-between', alignItems: 'flex-start' }}>
            <Box>
              <Typography variant="h4" gutterBottom>
                Targets
              </Typography>
              <Typography color="text.secondary">
                Drug discovery targets and campaigns
              </Typography>
            </Box>
            <Box sx={{ display: 'flex', gap: 1, alignItems: 'center' }}>
              {/* View mode toggle */}
              <ToggleButtonGroup
                value={viewMode}
                exclusive
                onChange={handleViewModeChange}
                size="small"
              >
                <ToggleButton value="grid">
                  <Tooltip title="Card view">
                    <ViewModule />
                  </Tooltip>
                </ToggleButton>
                <ToggleButton value="table">
                  <Tooltip title="Table view">
                    <ViewList />
                  </Tooltip>
                </ToggleButton>
              </ToggleButtonGroup>

              <Button
                component={Link}
                href="/"
                variant="outlined"
                startIcon={<Apps />}
              >
                All Apps
              </Button>
              <Button
                component={Link}
                href={routes.registry.suppliers()}
                variant="outlined"
                startIcon={<LocalShipping />}
              >
                Suppliers
              </Button>
              <Button
                component={Link}
                href={routes.registry.scaffoldExtensions()}
                variant="outlined"
                startIcon={<Hub />}
              >
                Scaffolds
              </Button>
              <Tooltip
                title={canContribute ? '' : 'Requires Contributor or Admin operating level'}
                arrow
              >
                <span>
                  <Button
                    variant="contained"
                    startIcon={<Add />}
                    onClick={() => setCreateDialogOpen(true)}
                    disabled={!canContribute}
                  >
                    New Target
                  </Button>
                </span>
              </Tooltip>
            </Box>
          </Box>
        </Box>

        <Box sx={{ flex: 1, overflow: 'hidden', minHeight: 0 }}>
          {viewMode === 'table' ? (
            <DataTable
              data={targets}
              columns={columns}
              loading={isLoading}
              onRowClick={(target) => router.push(routes.registry.target(target.id))}
              getRowKey={(row) => row.id}
              emptyMessage="No targets found"
              comfortable
              fillHeight
            />
          ) : (
            <TargetCardGrid
              targets={targets}
              loading={isLoading}
              onTargetClick={(target) => router.push(routes.registry.target(target.id))}
              emptyMessage="No targets found"
            />
          )}
        </Box>
      </Container>

      <TargetCreateDialog
        open={createDialogOpen}
        onClose={() => setCreateDialogOpen(false)}
        onCreated={handleTargetCreated}
      />

      <TargetDeletionDialog
        open={deleteTarget !== null}
        target={deleteTarget}
        onClose={() => setDeleteTarget(null)}
        onDeleted={handleDeletionComplete}
      />

      <Snackbar
        open={snackbarMessage !== null}
        autoHideDuration={4000}
        onClose={() => setSnackbarMessage(null)}
        message={snackbarMessage ?? ''}
      />
    </Box>
  );
}
