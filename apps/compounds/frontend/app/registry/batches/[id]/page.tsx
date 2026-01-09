'use client';

import { use } from 'react';
import { useRouter } from 'next/navigation';
import {
  Container,
  Typography,
  Box,
  Paper,
  Grid,
  Chip,
  Skeleton,
  Divider,
  Button,
  IconButton,
  Tooltip,
  Link as MuiLink,
} from '@mui/material';
import {
  Inventory,
  Medication,
  Science,
  Description,
  Download,
  OpenInNew,
} from '@mui/icons-material';
import { Breadcrumbs } from '@/components/Breadcrumbs';
import { DataTable, Column } from '@/components/DataTable';
import { useCompoundsApi } from '@/lib/api';
import { routes } from '@/lib/routes';
import { Batch, BatchQCFile, Compound, Target } from '@/types/models';

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

export default function BatchDetailPage({ params }: PageProps) {
  const { id } = use(params);
  const router = useRouter();
  const api = useCompoundsApi();

  const { data: batch, isLoading: batchLoading } = api.get<Batch>(
    `batches/${id}/`
  );
  const { data: qcFiles, isLoading: qcFilesLoading } = api.get<BatchQCFile[]>(
    `batch-qc-files/?batch=${id}`
  );
  const { data: compound } = api.get<Compound>(
    batch?.compound ? `compounds/${batch.compound}/` : null
  );
  const { data: target } = api.get<Target>(
    compound?.target ? `targets/${compound.target}/` : null
  );

  const columns: Column<BatchQCFile>[] = [
    {
      key: 'filename',
      label: 'File',
      sortable: true,
      searchable: true,
      render: (value, row) => (
        <Box sx={{ display: 'flex', alignItems: 'center', gap: 1 }}>
          <Description fontSize="small" color="action" />
          <Typography fontWeight={500}>{value || 'Unnamed file'}</Typography>
        </Box>
      ),
    },
    {
      key: 'comments',
      label: 'Comments',
      searchable: true,
      render: (value) =>
        value ? (
          <Typography
            sx={{
              maxWidth: 300,
              overflow: 'hidden',
              textOverflow: 'ellipsis',
              whiteSpace: 'nowrap',
            }}
            title={value}
          >
            {value}
          </Typography>
        ) : (
          '-'
        ),
    },
    {
      key: 'uploaded_at',
      label: 'Uploaded',
      sortable: true,
      width: 120,
      render: (value) =>
        value ? new Date(value).toLocaleDateString() : '-',
    },
    {
      key: 'file',
      label: 'Actions',
      width: 100,
      render: (value) =>
        value ? (
          <Box sx={{ display: 'flex', gap: 0.5 }}>
            <Tooltip title="Download">
              <IconButton
                size="small"
                href={value}
                target="_blank"
                onClick={(e) => e.stopPropagation()}
              >
                <Download fontSize="small" />
              </IconButton>
            </Tooltip>
            <Tooltip title="Open in new tab">
              <IconButton
                size="small"
                href={value}
                target="_blank"
                onClick={(e) => e.stopPropagation()}
              >
                <OpenInNew fontSize="small" />
              </IconButton>
            </Tooltip>
          </Box>
        ) : null,
    },
  ];

  return (
    <Container maxWidth="lg" sx={{ py: 3 }}>
      <Breadcrumbs
        items={[
          { label: 'Home', href: '/', icon: 'home' },
          { label: 'Targets', href: '/registry/targets', icon: 'target' },
          {
            label: target?.name || 'Target',
            href: compound?.target
              ? `/registry/targets/${compound.target}`
              : undefined,
            icon: 'target',
          },
          {
            label: compound?.formatted_id || 'Compound',
            href: batch?.compound
              ? `/registry/compounds/${batch.compound}`
              : undefined,
            icon: 'compound',
          },
          {
            label: batch ? `Batch #${batch.batch_number}` : 'Loading...',
            icon: 'batch',
          },
        ]}
      />

      {/* Batch header */}
      <Paper sx={{ p: 3, mb: 3 }}>
        {batchLoading ? (
          <>
            <Skeleton variant="text" width={300} height={40} />
            <Skeleton variant="text" width={200} />
          </>
        ) : batch ? (
          <>
            <Box sx={{ display: 'flex', alignItems: 'center', gap: 2, mb: 2 }}>
              <Inventory sx={{ fontSize: 48, color: 'info.main' }} />
              <Box>
                <Typography variant="h4">Batch #{batch.batch_number}</Typography>
                <Box sx={{ display: 'flex', gap: 1, mt: 0.5 }}>
                  {compound && (
                    <Chip
                      icon={<Medication fontSize="small" />}
                      label={compound.formatted_id}
                      size="small"
                      onClick={() =>
                        router.push(routes.registry.compound(compound.id))
                      }
                    />
                  )}
                  {target && (
                    <Chip
                      icon={<Science fontSize="small" />}
                      label={target.name}
                      size="small"
                      variant="outlined"
                      onClick={() =>
                        router.push(routes.registry.target(target.id))
                      }
                    />
                  )}
                </Box>
              </Box>
            </Box>

            <Divider sx={{ my: 2 }} />

            <Grid container spacing={3}>
              <Grid item xs={12} md={6}>
                <Typography variant="h6" gutterBottom>
                  Properties
                </Typography>
                <InfoRow
                  label="Amount"
                  value={
                    batch.amount ? `${parseFloat(batch.amount).toFixed(2)} mg` : null
                  }
                />
                <InfoRow label="Salt Code" value={batch.salt_code} />
                <InfoRow
                  label="MW (salt)"
                  value={batch.molecular_weight?.toFixed(2)}
                />
              </Grid>
              <Grid item xs={12} md={6}>
                <Typography variant="h6" gutterBottom>
                  Provenance
                </Typography>
                <InfoRow label="Supplier" value={batch.supplier_name} />
                <InfoRow label="Supplier Ref" value={batch.supplier_ref} />
                <InfoRow label="Lab Book" value={batch.labbook_number} />
                <InfoRow label="Page" value={batch.page_number} />
                <InfoRow
                  label="Registered"
                  value={
                    batch.registered_at
                      ? new Date(batch.registered_at).toLocaleString()
                      : null
                  }
                />
              </Grid>
            </Grid>

            {batch.comments && (
              <>
                <Divider sx={{ my: 2 }} />
                <Typography variant="h6" gutterBottom>
                  Comments
                </Typography>
                <Typography>{batch.comments}</Typography>
              </>
            )}
          </>
        ) : (
          <Typography color="error">Batch not found</Typography>
        )}
      </Paper>

      {/* QC Files table */}
      <DataTable
        data={qcFiles}
        columns={columns}
        loading={qcFilesLoading}
        getRowKey={(row) => row.id}
        title={qcFiles ? `${qcFiles.length} QC files` : undefined}
        emptyMessage="No QC files uploaded for this batch"
      />
    </Container>
  );
}
