'use client';

import { useCallback, useState } from 'react';
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
import {
  applyClarifyPick,
  CompoundSelector,
  NLPResponse,
  postNlpQuery,
} from '@/lib/compounds/nlp-api';
import { NLPResults } from './NLPResults';

const EXAMPLE_PROMPTS = [
  'Show all CDK4 compounds',
  'CDK4 compounds with HTRF IC50 < 100 nM',
  'mEGFR compounds with mEGFR TR-FRET WT IC50 < 10 uM but TM IC50 > 1 uM',
  'MYC compounds tested in AlphaLISA',
];

export function NLPPanel() {
  const [prompt, setPrompt] = useState('');
  const [response, setResponse] = useState<NLPResponse | null>(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);

  const submitPrompt = useCallback(async (text: string) => {
    if (!text.trim()) return;
    setLoading(true);
    setError(null);
    try {
      const result = await postNlpQuery({ prompt: text });
      setResponse(result);
    } catch (e) {
      setError(e instanceof Error ? e.message : String(e));
      setResponse(null);
    } finally {
      setLoading(false);
    }
  }, []);

  const handleSubmit = useCallback(
    () => submitPrompt(prompt),
    [prompt, submitPrompt],
  );

  const handleClarifyPick = useCallback(
    async (
      partial: CompoundSelector,
      field: string,
      pickedId: string,
      filterIndex?: number,
    ) => {
      const pinned = applyClarifyPick(partial, field, pickedId, filterIndex);
      setLoading(true);
      setError(null);
      try {
        const result = await postNlpQuery({ selector: pinned });
        setResponse(result);
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
      submitPrompt(text);
    },
    [submitPrompt],
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
