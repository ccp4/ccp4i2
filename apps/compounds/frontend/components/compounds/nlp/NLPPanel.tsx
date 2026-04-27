'use client';

import { useCallback, useEffect, useRef, useState } from 'react';
import {
  Alert,
  Box,
  Button,
  CircularProgress,
  Paper,
  Stack,
  TextField,
  Typography,
} from '@mui/material';
import { Send } from '@mui/icons-material';
import { usePathname, useRouter, useSearchParams } from 'next/navigation';
import { mutate } from 'swr';
import {
  applyClarifyPick,
  AssaySelector,
  CompoundSelector,
  NLPResponse,
  postNlpQuery,
  ProtocolCandidate,
  ScaffoldCandidate,
  SupplierCandidate,
  TargetCandidate,
  UnionCandidate,
  UserCandidate,
} from '@/lib/compounds/nlp-api';
import { SELECTIONS_LIST_KEY } from '@/lib/compounds/selections-api';
import { NLPResults } from './NLPResults';
import { SessionList } from './SessionList';

const EXAMPLE_PROMPTS = [
  'Show all CDK4 compounds',
  'CDK4 compounds with HTRF IC50 < 100 nM',
  'mEGFR compounds with mEGFR TR-FRET WT IC50 < 10 uM but TM IC50 > 1 uM',
  'MYC compounds tested in AlphaLISA',
];

export function NLPPanel() {
  const router = useRouter();
  const pathname = usePathname();
  const searchParams = useSearchParams();

  const [prompt, setPrompt] = useState('');
  const [response, setResponse] = useState<NLPResponse | null>(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);

  // Tracks the most-recently-replayed `?q=` value so the URL-state
  // restoration effect doesn't re-fire on every URL update we initiate.
  const lastReplayedQuery = useRef<string | null>(null);

  const submitPrompt = useCallback(async (text: string) => {
    if (!text.trim()) return;
    setLoading(true);
    setError(null);
    try {
      const result = await postNlpQuery({ prompt: text });
      setResponse(result);
      mutate(SELECTIONS_LIST_KEY);
    } catch (e) {
      setError(e instanceof Error ? e.message : String(e));
      setResponse(null);
    } finally {
      setLoading(false);
    }
  }, []);

  const writeQueryParam = useCallback(
    (text: string) => {
      lastReplayedQuery.current = text;
      const params = new URLSearchParams(searchParams?.toString() ?? '');
      params.set('q', text);
      router.replace(`${pathname}?${params.toString()}`);
    },
    [router, pathname, searchParams],
  );

  // URL-state restoration: if the user lands here with `?q=<prompt>`
  // (typically from the browser back-button after navigating to
  // /assays/aggregate), re-run that prompt so the result is restored.
  useEffect(() => {
    const q = searchParams?.get('q');
    if (!q || q === lastReplayedQuery.current) return;
    lastReplayedQuery.current = q;
    setPrompt(q);
    submitPrompt(q);
  }, [searchParams, submitPrompt]);

  const handleSubmit = useCallback(() => {
    const trimmed = prompt.trim();
    if (!trimmed) return;
    writeQueryParam(trimmed);
    submitPrompt(trimmed);
  }, [prompt, submitPrompt, writeQueryParam]);

  const handleClarifyPick = useCallback(
    async (
      partial: CompoundSelector | AssaySelector,
      field: string,
      candidate: TargetCandidate | ProtocolCandidate | UserCandidate | SupplierCandidate | UnionCandidate | ScaffoldCandidate,
      filterIndex?: number,
      scaffoldIndex?: number,
    ) => {
      const pinned = applyClarifyPick(partial, field, candidate, filterIndex, scaffoldIndex);
      setLoading(true);
      setError(null);
      try {
        const result = await postNlpQuery({ selector: pinned });
        setResponse(result);
        mutate(SELECTIONS_LIST_KEY);
      } catch (e) {
        setError(e instanceof Error ? e.message : String(e));
        setResponse(null);
      } finally {
        setLoading(false);
      }
    },
    [],
  );

  const handleExampleClick = useCallback(
    (text: string) => {
      setPrompt(text);
      writeQueryParam(text);
      submitPrompt(text);
    },
    [submitPrompt, writeQueryParam],
  );

  return (
    <Stack spacing={3}>
      <Paper elevation={1} sx={{ p: 3 }}>
        <Box sx={{ display: 'flex', gap: 1.5, alignItems: 'flex-start' }}>
          <TextField
            fullWidth
            multiline
            minRows={1}
            maxRows={3}
            placeholder="Describe the compounds you want (e.g. 'CDK4 compounds with HTRF IC50 < 100 nM')"
            value={prompt}
            onChange={(e) => setPrompt(e.target.value)}
            onKeyDown={(e) => {
              if (e.key === 'Enter' && !e.shiftKey) {
                e.preventDefault();
                handleSubmit();
              }
            }}
            disabled={loading}
            autoFocus
          />
          <Button
            variant="contained"
            onClick={handleSubmit}
            disabled={loading || !prompt.trim()}
            endIcon={loading ? <CircularProgress size={16} color="inherit" /> : <Send />}
            sx={{ minWidth: 110, alignSelf: 'stretch' }}
          >
            Ask
          </Button>
        </Box>

        {!response && !loading && (
          <Box sx={{ mt: 2 }}>
            <Typography variant="caption" color="text.secondary" sx={{ display: 'block', mb: 1 }}>
              Try one of these:
            </Typography>
            <Stack direction="row" spacing={1} sx={{ flexWrap: 'wrap', gap: 1 }}>
              {EXAMPLE_PROMPTS.map((example) => (
                <Button
                  key={example}
                  size="small"
                  variant="outlined"
                  onClick={() => handleExampleClick(example)}
                  sx={{ textTransform: 'none', fontSize: '0.8rem' }}
                >
                  {example}
                </Button>
              ))}
            </Stack>
          </Box>
        )}
      </Paper>

      <SessionList />

      {error && (
        <Alert severity="error">
          Failed to reach the NLP endpoint: {error}
        </Alert>
      )}

      {response && (
        <NLPResults response={response} onClarifyPick={handleClarifyPick} />
      )}
    </Stack>
  );
}
