'use client';

/**
 * Scaffold-extensions management page — exposes the chemist's saved
 * `ScaffoldExtension` rows (from slice 17's "Define this fragment"
 * affordance) alongside the curated seed catalog. Surfaces the raw
 * SMARTS as stored, the canonical aromatic SMARTS RDKit produces from
 * it (so the chemist can see whether their Kekulé SMILES paste was
 * normalised into a usable query), and a thumbnail render. Creator-
 * scoped delete; seed entries are read-only.
 */

import { useCallback, useState } from 'react';
import {
  Alert,
  Box,
  Chip,
  CircularProgress,
  Container,
  IconButton,
  Paper,
  Stack,
  Table,
  TableBody,
  TableCell,
  TableHead,
  TableRow,
  Tooltip,
  Typography,
} from '@mui/material';
import { DeleteOutline } from '@mui/icons-material';
import useSWR from 'swr';
import { PageHeader } from '@/components/compounds/PageHeader';
import { ConfirmDialog } from '@/components/compounds/ConfirmDialog';
import { MoleculeChip } from '@/components/compounds/MoleculeView';
import { routes } from '@/lib/compounds/routes';
import { useAuth } from '@/lib/compounds/auth-context';
import {
  SCAFFOLD_CATALOG_KEY,
  ScaffoldCatalogEntry,
  deleteScaffoldExtension,
  listScaffoldCatalog,
} from '@/lib/compounds/scaffold-extensions-api';

export default function ScaffoldExtensionsPage() {
  const { user } = useAuth();
  const { data, error, isLoading, mutate } = useSWR<ScaffoldCatalogEntry[]>(
    SCAFFOLD_CATALOG_KEY,
    () => listScaffoldCatalog(),
    { revalidateOnFocus: false },
  );

  const [confirmDelete, setConfirmDelete] = useState<ScaffoldCatalogEntry | null>(null);
  const [deleteError, setDeleteError] = useState<string | null>(null);

  const handleDelete = useCallback(async () => {
    if (!confirmDelete?.id) return;
    setDeleteError(null);
    try {
      await deleteScaffoldExtension(confirmDelete.id);
      setConfirmDelete(null);
      mutate();
    } catch (e) {
      setDeleteError(e instanceof Error ? e.message : String(e));
    }
  }, [confirmDelete, mutate]);

  return (
    <Box sx={{ display: 'flex', flexDirection: 'column', minHeight: '100vh' }}>
      <Container maxWidth="lg" sx={{ py: 2, flex: 1 }}>
        <PageHeader
          breadcrumbs={[
            { label: 'Home', href: routes.home(), icon: 'home' },
            { label: 'Registry', href: routes.registry.targets() },
            { label: 'Scaffold extensions' },
          ]}
        />

        <Typography variant="h5" sx={{ mt: 2, mb: 1 }}>
          Scaffold extensions
        </Typography>
        <Typography variant="body2" color="text.secondary" sx={{ mb: 2 }}>
          User-defined chemotypes available in NLP queries
          (<em>&ldquo;coloured by …&rdquo;</em>, <em>&ldquo;containing …&rdquo;</em>) and on
          the scatter view&apos;s chip strip. Raw SMARTS is stored verbatim;
          canonical SMARTS is what RDKit actually matches with after
          aromatic perception. Seed catalog rows are read-only.
        </Typography>

        {error && (
          <Alert severity="error" sx={{ mb: 2 }}>
            Failed to load catalog: {String(error)}
          </Alert>
        )}

        {isLoading ? (
          <Box sx={{ display: 'flex', justifyContent: 'center', py: 4 }}>
            <CircularProgress />
          </Box>
        ) : (
          <Paper variant="outlined">
            <Table size="small">
              <TableHead>
                <TableRow>
                  <TableCell>Name</TableCell>
                  <TableCell>Source</TableCell>
                  <TableCell>Structure</TableCell>
                  <TableCell>Raw SMARTS (as stored)</TableCell>
                  <TableCell>Canonical SMARTS (matched)</TableCell>
                  <TableCell>Added by</TableCell>
                  <TableCell sx={{ width: 56 }} />
                </TableRow>
              </TableHead>
              <TableBody>
                {(data ?? []).map((row) => {
                  const isExtension = row.source.startsWith('extension');
                  const ownedByMe =
                    isExtension && row.created_by != null && user?.username === row.created_by;
                  const canonical = row.canonical_smarts;
                  const rendered = canonical ?? row.smarts;
                  return (
                    <TableRow key={`${row.source}:${row.name}:${row.id ?? ''}`} hover>
                      <TableCell>
                        <Typography variant="body2" sx={{ fontWeight: 500 }}>
                          {row.name}
                        </Typography>
                        {row.aliases.length > 0 && (
                          <Typography variant="caption" color="text.secondary">
                            {row.aliases.join(', ')}
                          </Typography>
                        )}
                      </TableCell>
                      <TableCell>
                        <Chip
                          label={
                            row.source === 'seed'
                              ? 'Seed'
                              : row.source === 'extension:project'
                                ? `Project · ${row.target_name ?? ''}`
                                : 'Shared'
                          }
                          size="small"
                          variant={row.source === 'seed' ? 'outlined' : 'filled'}
                        />
                      </TableCell>
                      <TableCell sx={{ width: 120 }}>
                        <MoleculeChip smiles={rendered} size={100} disableHover />
                      </TableCell>
                      <TableCell>
                        <Box
                          component="code"
                          sx={{
                            fontFamily: 'monospace',
                            fontSize: '0.78rem',
                            wordBreak: 'break-all',
                          }}
                        >
                          {row.smarts}
                        </Box>
                      </TableCell>
                      <TableCell>
                        {canonical ? (
                          <Box
                            component="code"
                            sx={{
                              fontFamily: 'monospace',
                              fontSize: '0.78rem',
                              wordBreak: 'break-all',
                              color:
                                canonical === row.smarts
                                  ? 'text.secondary'
                                  : 'text.primary',
                            }}
                          >
                            {canonical}
                          </Box>
                        ) : (
                          <Typography variant="caption" color="error">
                            Unparseable — RDKit could not produce a canonical form
                          </Typography>
                        )}
                      </TableCell>
                      <TableCell>
                        {row.created_by ? (
                          <Stack direction="row" spacing={0.5} alignItems="center">
                            <Typography variant="caption">{row.created_by}</Typography>
                            {ownedByMe && (
                              <Chip label="me" size="small" variant="outlined" />
                            )}
                          </Stack>
                        ) : (
                          <Typography variant="caption" color="text.secondary">
                            —
                          </Typography>
                        )}
                      </TableCell>
                      <TableCell>
                        {ownedByMe && row.id && (
                          <Tooltip title="Delete this extension">
                            <IconButton
                              size="small"
                              onClick={() => setConfirmDelete(row)}
                            >
                              <DeleteOutline fontSize="small" />
                            </IconButton>
                          </Tooltip>
                        )}
                      </TableCell>
                    </TableRow>
                  );
                })}
              </TableBody>
            </Table>
          </Paper>
        )}

        <ConfirmDialog
          open={confirmDelete !== null}
          title="Delete scaffold extension?"
          message={
            <>
              Remove <strong>{confirmDelete?.name}</strong> from the catalog?
              This can&apos;t be undone, and any saved selections that
              referenced it by name will silently fail to colour by it
              after this.
            </>
          }
          onConfirm={handleDelete}
          onCancel={() => {
            setConfirmDelete(null);
            setDeleteError(null);
          }}
        />
        {deleteError && (
          <Alert severity="error" sx={{ mt: 2 }}>
            Delete failed: {deleteError}
          </Alert>
        )}
      </Container>
    </Box>
  );
}
