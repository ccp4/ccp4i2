'use client';

import { use, useState, useCallback } from 'react';
import { useRouter } from 'next/navigation';
import {
  Container,
  Typography,
  Box,
  Paper,
  Grid2 as Grid,
  Chip,
  Skeleton,
  Divider,
  Button,
  Tooltip,
  IconButton,
} from '@mui/material';
import { Add, ChevronLeft, ChevronRight, ContentCopy, Check, Edit, Inventory, Medication, Science, TableChart } from '@mui/icons-material';
import Link from 'next/link';
import { PageHeader } from '@/components/compounds/PageHeader';
import { DataTable, Column } from '@/components/compounds/DataTable';
import { MoleculeView } from '@/components/compounds/MoleculeView';
import { BatchCreateDialog } from '@/components/compounds/BatchCreateDialog';
import { CompoundEditDialog } from '@/components/compounds/CompoundEditDialog';
import { AliasEditor } from '@/components/compounds/AliasEditor';
import { useCompoundsApi } from '@/lib/compounds/api';
import { useAuth } from '@/lib/compounds/auth-context';
import { routes } from '@/lib/compounds/routes';
import { Compound, Batch, Target } from '@/types/compounds/models';

interface PageProps {
  params: Promise<{ id: string }>;
}

function InfoRow({ label, value }: { label: string; value: React.ReactNode }) {
  return (
    <Box sx={{ display: 'flex', py: 0.5 }}>
      <Typography
        color="text.secondary"
        sx={{ minWidth: 140, fontWeight: 500 }}
      >
        {label}:
      </Typography>
      <Typography component="div">{value ?? '-'}</Typography>
    </Box>
  );
}

export default function CompoundDetailPage({ params }: PageProps) {
  const { id } = use(params);
  const router = useRouter();
  const api = useCompoundsApi();
  const { canContribute, canAdminister } = useAuth();

  const [batchDialogOpen, setBatchDialogOpen] = useState(false);
  const [editDialogOpen, setEditDialogOpen] = useState(false);
  const [smilesCopied, setSmilesCopied] = useState(false);

  const { data: compound, isLoading: compoundLoading, mutate: mutateCompound } = api.get<Compound>(
    `compounds/${id}/`
  );
  const { data: batches, isLoading: batchesLoading, mutate: mutateBatches } = api.get<Batch[]>(
    `batches/?compound=${id}`
  );
  const { data: target } = api.get<Target>(
    compound?.target ? `targets/${compound.target}/` : null
  );
  const { data: adjacent } = api.get<{
    previous: { id: string; formatted_id: string } | null;
    next: { id: string; formatted_id: string } | null;
  }>(`compounds/${id}/adjacent/`);

  const handleCopySmiles = useCallback(async () => {
    if (!compound?.smiles) return;
    try {
      await navigator.clipboard.writeText(compound.smiles);
      setSmilesCopied(true);
      setTimeout(() => setSmilesCopied(false), 2000);
    } catch (err) {
      console.error('Failed to copy:', err);
    }
  }, [compound?.smiles]);

  const columns: Column<Batch>[] = [
    {
      key: 'batch_number',
      label: 'Batch',
      sortable: true,
      width: 80,
      render: (value, row) => (
        <Box sx={{ display: 'flex', alignItems: 'center', gap: 1 }}>
          <Inventory fontSize="small" color="info" />
          <Typography fontWeight={500}>#{value}</Typography>
        </Box>
      ),
    },
    {
      key: 'supplier_name',
      label: 'Supplier',
      sortable: true,
      searchable: true,
      render: (value) => value || '-',
    },
    {
      key: 'supplier_ref',
      label: 'Supplier Ref',
      searchable: true,
      render: (value) =>
        value ? (
          <Typography fontFamily="monospace" fontSize="0.85rem">
            {value}
          </Typography>
        ) : (
          '-'
        ),
    },
    {
      key: 'amount',
      label: 'Amount (mg)',
      sortable: true,
      width: 100,
      render: (value) => (value ? parseFloat(value).toFixed(2) : '-'),
    },
    {
      key: 'qc_file_count',
      label: 'QC Files',
      sortable: true,
      width: 80,
      render: (value) =>
        value ? (
          <Chip label={value} size="small" color="success" variant="outlined" />
        ) : (
          '-'
        ),
    },
    {
      key: 'registered_at',
      label: 'Date',
      sortable: true,
      width: 100,
      render: (value) =>
        value ? new Date(value).toLocaleDateString() : '-',
    },
  ];

  return (
    <Container maxWidth="lg" sx={{ py: 3 }}>
      <PageHeader
        breadcrumbs={[
          { label: 'Home', href: routes.home(), icon: 'home' },
          { label: 'Targets', href: routes.registry.targets(), icon: 'target' },
          {
            label: target?.name || 'Target',
            href: compound?.target
              ? routes.registry.target(compound.target)
              : undefined,
            icon: 'target',
          },
          { label: compound?.formatted_id || 'Loading...', icon: 'compound' },
        ]}
      />

      {/* Compound header */}
      <Paper sx={{ p: 3, mb: 3, position: 'relative' }}>
        {/* Navigation arrows - fixed at top to avoid jumping when content changes */}
        {adjacent && (
          <>
            <Tooltip title={adjacent.previous ? `Previous: ${adjacent.previous.formatted_id}` : 'No previous compound'}>
              <span style={{ position: 'absolute', left: 8, top: 24 }}>
                <IconButton
                  onClick={() => adjacent.previous && router.push(routes.registry.compound(adjacent.previous.id))}
                  disabled={!adjacent.previous}
                  size="large"
                  sx={{
                    bgcolor: 'action.hover',
                    '&:hover': { bgcolor: 'action.selected' },
                    '&.Mui-disabled': { bgcolor: 'transparent' }
                  }}
                >
                  <ChevronLeft />
                </IconButton>
              </span>
            </Tooltip>
            <Tooltip title={adjacent.next ? `Next: ${adjacent.next.formatted_id}` : 'No next compound'}>
              <span style={{ position: 'absolute', right: 8, top: 24 }}>
                <IconButton
                  onClick={() => adjacent.next && router.push(routes.registry.compound(adjacent.next.id))}
                  disabled={!adjacent.next}
                  size="large"
                  sx={{
                    bgcolor: 'action.hover',
                    '&:hover': { bgcolor: 'action.selected' },
                    '&.Mui-disabled': { bgcolor: 'transparent' }
                  }}
                >
                  <ChevronRight />
                </IconButton>
              </span>
            </Tooltip>
          </>
        )}

        {compoundLoading ? (
          <>
            <Skeleton variant="text" width={300} height={40} />
            <Skeleton variant="text" width="100%" height={60} />
          </>
        ) : compound ? (
          <>
            <Box sx={{ display: 'flex', alignItems: 'center', justifyContent: 'space-between', mb: 2, mx: 5 }}>
              <Box sx={{ display: 'flex', alignItems: 'center', gap: 2 }}>
                <Medication sx={{ fontSize: 48, color: 'secondary.main' }} />
                <Box>
                  <Typography variant="h4" fontFamily="monospace">
                    {compound.formatted_id}
                  </Typography>
                  <Box sx={{ display: 'flex', gap: 1, mt: 0.5 }}>
                    {compound.stereo_comment &&
                      compound.stereo_comment !== 'unset' && (
                        <Chip
                          label={compound.stereo_comment}
                          size="small"
                          color="primary"
                          variant="outlined"
                        />
                      )}
                    {target && (
                      <Chip
                        icon={<Science fontSize="small" />}
                        label={target.name}
                        size="small"
                        onClick={() =>
                          router.push(routes.registry.target(target.id))
                        }
                      />
                    )}
                  </Box>
                </Box>
              </Box>
              <Box sx={{ display: 'flex', gap: 1 }}>
                {canAdminister && (
                  <Tooltip title="Edit compound details (Admin only)">
                    <Button
                      variant="outlined"
                      startIcon={<Edit />}
                      onClick={() => setEditDialogOpen(true)}
                    >
                      Edit
                    </Button>
                  </Tooltip>
                )}
                <Button
                  component={Link}
                  href={routes.assays.aggregate({ compound: compound.formatted_id })}
                  variant="outlined"
                  startIcon={<TableChart />}
                >
                  View Assay Data
                </Button>
              </Box>
            </Box>

            <Divider sx={{ my: 2 }} />

            <Grid container spacing={3}>
              <Grid size={{ xs: 12, md: 6 }}>
                <Typography variant="h6" gutterBottom>
                  Structure
                </Typography>
                <Box sx={{ display: 'flex', gap: 3, mb: 2 }}>
                  <MoleculeView smiles={compound.smiles} width={200} height={200} />
                  <Box sx={{ flex: 1 }}>
                    <InfoRow
                      label="SMILES"
                      value={
                        <Box sx={{ display: 'flex', alignItems: 'flex-start', gap: 0.5 }}>
                          <Typography
                            fontFamily="monospace"
                            fontSize="0.85rem"
                            sx={{ wordBreak: 'break-all' }}
                          >
                            {compound.smiles}
                          </Typography>
                          <Tooltip title={smilesCopied ? 'Copied!' : 'Copy SMILES'}>
                            <IconButton
                              size="small"
                              onClick={handleCopySmiles}
                              sx={{
                                p: 0.5,
                                color: smilesCopied ? 'success.main' : 'action.active',
                                flexShrink: 0,
                              }}
                            >
                              {smilesCopied ? <Check fontSize="small" /> : <ContentCopy fontSize="small" />}
                            </IconButton>
                          </Tooltip>
                        </Box>
                      }
                    />
                    <InfoRow
                      label="MW"
                      value={compound.molecular_weight?.toFixed(2)}
                    />
                    {compound.inchi && (
                      <InfoRow
                        label="InChI"
                        value={
                          <Typography
                            fontFamily="monospace"
                            fontSize="0.75rem"
                            sx={{ wordBreak: 'break-all' }}
                          >
                            {compound.inchi}
                          </Typography>
                        }
                      />
                    )}
                  </Box>
                </Box>
              </Grid>
              <Grid size={{ xs: 12, md: 6 }}>
                <Typography variant="h6" gutterBottom>
                  Provenance
                </Typography>
                <InfoRow label="Supplier" value={compound.supplier_name} />
                <InfoRow label="Supplier Ref" value={compound.supplier_ref} />
                <InfoRow label="Lab Book" value={compound.labbook_number} />
                <InfoRow label="Page" value={compound.page_number} />
                <InfoRow label="Barcode" value={compound.barcode} />
                <InfoRow
                  label="Registered"
                  value={
                    compound.registered_at
                      ? new Date(compound.registered_at).toLocaleString()
                      : null
                  }
                />
                <InfoRow
                  label="By"
                  value={
                    compound.registered_by_email ||
                    compound.legacy_registered_by
                  }
                />
              </Grid>
            </Grid>

            {compound.comments && (
              <>
                <Divider sx={{ my: 2 }} />
                <Typography variant="h6" gutterBottom>
                  Comments
                </Typography>
                <Typography>{compound.comments}</Typography>
              </>
            )}

            {/* Aliases section */}
            <Divider sx={{ my: 2 }} />
            <AliasEditor
              compoundId={compound.id}
              aliases={compound.aliases || []}
              canEdit={canContribute}
              onUpdate={() => mutateCompound()}
            />
          </>
        ) : (
          <Typography color="error">Compound not found</Typography>
        )}
      </Paper>

      {/* Batches table */}
      <DataTable
        data={batches}
        columns={columns}
        loading={batchesLoading}
        onRowClick={(batch) => router.push(routes.registry.batch(batch.id))}
        getRowKey={(row) => row.id}
        title={batches ? `${batches.length} batches` : undefined}
        emptyMessage="No batches registered for this compound"
        headerAction={
          compound && (
            <Tooltip title={canContribute ? '' : 'Requires Contributor or Admin operating level'} arrow>
              <span>
                <Button
                  variant="contained"
                  size="small"
                  startIcon={<Add />}
                  onClick={() => setBatchDialogOpen(true)}
                  disabled={!canContribute}
                >
                  Add Batch
                </Button>
              </span>
            </Tooltip>
          )
        }
      />

      {/* Batch creation dialog */}
      {compound && (
        <BatchCreateDialog
          open={batchDialogOpen}
          onClose={() => setBatchDialogOpen(false)}
          onCreated={(newBatch) => {
            mutateBatches();
            router.push(routes.registry.batch(newBatch.id));
          }}
          compoundId={compound.id}
          compoundFormattedId={compound.formatted_id}
        />
      )}

      {/* Compound edit dialog (admin only) */}
      {compound && (
        <CompoundEditDialog
          open={editDialogOpen}
          onClose={() => setEditDialogOpen(false)}
          onSaved={() => {
            mutateCompound();
          }}
          compound={compound}
        />
      )}
    </Container>
  );
}
