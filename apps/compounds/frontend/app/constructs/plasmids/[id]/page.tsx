'use client';

import { use, useState } from 'react';
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
} from '@mui/material';
import { Download, Science, Link as LinkIcon, Delete } from '@mui/icons-material';
import Link from 'next/link';
import { PageHeader } from '@/components/compounds/PageHeader';
import { DataTable, Column } from '@/components/compounds/DataTable';
import { SeqVizViewer } from '@/components/compounds/SeqVizViewer';
import { ConfirmDialog } from '@/components/compounds/ConfirmDialog';
import { useCompoundsApi } from '@/lib/compounds/api';
import { routes } from '@/lib/compounds/routes';
import {
  PlasmidDetail,
  CassetteUseDetail,
  SequencingResult,
  GenbankContentResponse,
} from '@/types/compounds/constructs';

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

export default function PlasmidDetailPage({ params }: PageProps) {
  const { id } = use(params);
  const router = useRouter();
  const api = useCompoundsApi();
  const [deleteDialogOpen, setDeleteDialogOpen] = useState(false);
  const [deleting, setDeleting] = useState(false);

  // Fetch plasmid detail
  const { data: plasmid, isLoading: plasmidLoading } =
    api.get<PlasmidDetail>(`plasmids/${id}/`);

  // Fetch GenBank content (only if plasmid has a file)
  const { data: genbankResponse, isLoading: genbankLoading } =
    api.get<GenbankContentResponse>(
      plasmid?.genbank_file ? `plasmids/${id}/genbank_content/` : null
    );

  const handleDelete = async () => {
    setDeleting(true);
    try {
      const response = await fetch(`/api/proxy/compounds/plasmids/${id}/`, {
        method: 'DELETE',
      });
      if (response.ok) {
        router.push(routes.constructs.list());
      }
    } catch (error) {
      console.error('Failed to delete plasmid:', error);
    } finally {
      setDeleting(false);
    }
  };

  const cassetteColumns: Column<CassetteUseDetail>[] = [
    {
      key: 'cassette_display',
      label: 'Cassette',
      sortable: true,
      render: (value) => (
        <Typography fontFamily="monospace">{value || '-'}</Typography>
      ),
    },
    {
      key: 'expression_tags',
      label: 'Expression Tags',
      render: (_, row) => {
        if (!row.expression_tags || row.expression_tags.length === 0) {
          return '-';
        }
        return (
          <Box sx={{ display: 'flex', gap: 0.5, flexWrap: 'wrap' }}>
            {row.expression_tags.map((tag) => (
              <Chip
                key={tag.id}
                label={`${tag.expression_tag_type_name} (${tag.location_name})`}
                size="small"
                variant="outlined"
              />
            ))}
          </Box>
        );
      },
    },
    {
      key: 'alignment_file',
      label: 'Alignment',
      width: 100,
      render: (value) =>
        value ? (
          <Chip label="Has file" size="small" color="success" variant="outlined" />
        ) : (
          '-'
        ),
    },
  ];

  const sequencingColumns: Column<SequencingResult>[] = [
    {
      key: 'filename',
      label: 'File',
      sortable: true,
      render: (value, row) => (
        <Box sx={{ display: 'flex', alignItems: 'center', gap: 1 }}>
          <Typography>{value || 'Unknown file'}</Typography>
          {row.file && (
            <Button
              size="small"
              startIcon={<Download />}
              href={row.file}
              target="_blank"
              sx={{ minWidth: 'auto', p: 0.5 }}
            >
              Download
            </Button>
          )}
        </Box>
      ),
    },
    {
      key: 'cassette_display',
      label: 'Cassette',
      sortable: true,
    },
    {
      key: 'created_at',
      label: 'Date',
      sortable: true,
      width: 120,
      render: (value) =>
        value ? new Date(value).toLocaleDateString() : '-',
    },
  ];

  return (
    <Container maxWidth="lg" sx={{ py: 3 }}>
      <PageHeader
        breadcrumbs={[
          { label: 'Home', href: routes.home(), icon: 'home' },
          { label: 'Constructs', href: routes.constructs.list() },
          { label: plasmid?.formatted_id || 'Loading...' },
        ]}
      />

      {/* Plasmid header */}
      <Paper sx={{ p: 3, mb: 3 }}>
        {plasmidLoading ? (
          <>
            <Skeleton variant="text" width={300} height={40} />
            <Skeleton variant="text" width={200} />
          </>
        ) : plasmid ? (
          <>
            <Box
              sx={{
                display: 'flex',
                alignItems: 'center',
                justifyContent: 'space-between',
                mb: 2,
              }}
            >
              <Box sx={{ display: 'flex', alignItems: 'center', gap: 2 }}>
                <Science sx={{ fontSize: 48, color: 'primary.main' }} />
                <Box>
                  <Typography variant="h4" fontFamily="monospace">
                    {plasmid.formatted_id}
                  </Typography>
                  <Typography color="text.secondary">{plasmid.name}</Typography>
                </Box>
              </Box>
              <Box sx={{ display: 'flex', gap: 1 }}>
                {plasmid.genbank_file_url && (
                  <Button
                    variant="outlined"
                    startIcon={<Download />}
                    href={plasmid.genbank_file_url}
                    download
                  >
                    Download GenBank
                  </Button>
                )}
                <Button
                  variant="outlined"
                  color="error"
                  startIcon={<Delete />}
                  onClick={() => setDeleteDialogOpen(true)}
                >
                  Delete
                </Button>
              </Box>
            </Box>

            <Divider sx={{ my: 2 }} />

            <Grid container spacing={3}>
              <Grid size={{ xs: 12, md: 6 }}>
                <InfoRow label="Project" value={plasmid.project_name} />
                <InfoRow
                  label="Parent Plasmid"
                  value={
                    plasmid.parent ? (
                      <Link
                        href={routes.constructs.plasmid(plasmid.parent)}
                        style={{ color: 'inherit' }}
                      >
                        <Box
                          component="span"
                          sx={{
                            display: 'flex',
                            alignItems: 'center',
                            gap: 0.5,
                            color: 'primary.main',
                          }}
                        >
                          <LinkIcon fontSize="small" />
                          {plasmid.parent_formatted_id}
                        </Box>
                      </Link>
                    ) : null
                  }
                />
              </Grid>
              <Grid size={{ xs: 12, md: 6 }}>
                <InfoRow
                  label="Created"
                  value={new Date(plasmid.created_at).toLocaleString()}
                />
                <InfoRow label="By" value={plasmid.created_by_email} />
              </Grid>
            </Grid>
          </>
        ) : (
          <Typography color="error">Plasmid not found</Typography>
        )}
      </Paper>

      {/* Sequence Viewer */}
      <Box sx={{ mb: 3 }}>
        <Typography variant="h6" gutterBottom>
          Sequence Map
        </Typography>
        <SeqVizViewer
          genbankContent={genbankResponse?.content || null}
          loading={plasmidLoading || genbankLoading}
          height={500}
          filename={plasmid?.genbank_file?.split('/').pop()}
          downloadUrl={plasmid?.genbank_file_url || undefined}
        />
      </Box>

      {/* Cassettes Table */}
      <Box sx={{ mb: 3 }}>
        <Typography variant="h6" gutterBottom>
          Cassettes
        </Typography>
        <DataTable
          data={plasmid?.cassette_uses}
          columns={cassetteColumns}
          loading={plasmidLoading}
          getRowKey={(row) => row.id}
          emptyMessage="No cassettes registered for this plasmid"
        />
      </Box>

      {/* Sequencing Results Table */}
      <Box sx={{ mb: 3 }}>
        <Typography variant="h6" gutterBottom>
          Sequencing Results
        </Typography>
        <DataTable
          data={plasmid?.sequencing_results}
          columns={sequencingColumns}
          loading={plasmidLoading}
          getRowKey={(row) => row.id}
          emptyMessage="No sequencing results available"
        />
      </Box>

      <ConfirmDialog
        open={deleteDialogOpen}
        title="Delete Plasmid"
        message={
          <>
            Are you sure you want to delete{' '}
            <strong>{plasmid?.formatted_id}</strong> ({plasmid?.name})?
            <br />
            <br />
            This will also delete any associated GenBank files and cannot be undone.
          </>
        }
        confirmLabel="Delete"
        loading={deleting}
        onConfirm={handleDelete}
        onCancel={() => setDeleteDialogOpen(false)}
      />
    </Container>
  );
}
