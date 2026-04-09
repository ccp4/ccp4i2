'use client';

import React, { useEffect, useState } from 'react';
import ReactMarkdown from 'react-markdown';
import remarkGfm from 'remark-gfm';
import {
  Box,
  Typography,
  Paper,
  CircularProgress,
  Alert,
  Table,
  TableBody,
  TableCell,
  TableContainer,
  TableHead,
  TableRow,
  Divider,
  useTheme,
} from '@mui/material';

interface MarkdownDocProps {
  /** Path to markdown file relative to public/docs/ */
  src: string;
  /** Optional title override (default: extracted from first # heading) */
  title?: string;
}

/**
 * Component that fetches and renders a markdown file.
 *
 * Markdown files should be placed in public/docs/ during the Docker build.
 * Supports GitHub-flavored markdown (tables, code blocks, etc).
 */
export function MarkdownDoc({ src, title }: MarkdownDocProps) {
  const theme = useTheme();
  const isDarkMode = theme.palette.mode === 'dark';
  const [content, setContent] = useState<string | null>(null);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);

  useEffect(() => {
    const fetchMarkdown = async () => {
      setLoading(true);
      setError(null);
      try {
        const response = await fetch(`/docs/${src}`);
        if (!response.ok) {
          throw new Error(`Failed to load documentation: ${response.statusText}`);
        }
        const text = await response.text();
        setContent(text);
      } catch (err) {
        setError(err instanceof Error ? err.message : 'Failed to load documentation');
      } finally {
        setLoading(false);
      }
    };

    fetchMarkdown();
  }, [src]);

  if (loading) {
    return (
      <Box sx={{ display: 'flex', justifyContent: 'center', py: 4 }}>
        <CircularProgress size={24} />
      </Box>
    );
  }

  if (error) {
    return (
      <Alert severity="error" sx={{ my: 2 }}>
        {error}
      </Alert>
    );
  }

  if (!content) {
    return null;
  }

  // Theme-aware colors for code blocks and table headers
  const codeBlockBg = isDarkMode ? 'grey.900' : 'grey.100';
  const tableHeaderBg = isDarkMode ? 'grey.800' : 'grey.50';
  const inlineCodeBg = isDarkMode ? 'rgba(255,255,255,0.1)' : 'rgba(0,0,0,0.06)';

  return (
    <Box
      sx={{
        '& h1': { display: 'none' }, // Hide h1 since we show title separately
        '& h2': {
          fontSize: '1.25rem',
          fontWeight: 600,
          mt: 3,
          mb: 1.5,
          pb: 0.5,
          borderBottom: '1px solid',
          borderColor: 'divider',
        },
        '& h3': {
          fontSize: '1rem',
          fontWeight: 600,
          mt: 2,
          mb: 1,
        },
        '& p': {
          fontSize: '0.875rem',
          lineHeight: 1.6,
          mb: 1.5,
        },
        '& ul, & ol': {
          pl: 3,
          mb: 1.5,
          '& li': {
            fontSize: '0.875rem',
            mb: 0.5,
          },
        },
        '& code': {
          fontFamily: 'monospace',
          fontSize: '0.8rem',
          bgcolor: inlineCodeBg,
          color: 'text.primary',
          px: 0.5,
          py: 0.25,
          borderRadius: 0.5,
        },
        '& pre': {
          bgcolor: codeBlockBg,
          color: 'text.primary',
          p: 1.5,
          borderRadius: 1,
          overflow: 'auto',
          mb: 2,
          '& code': {
            bgcolor: 'transparent',
            p: 0,
          },
        },
        '& hr': {
          my: 3,
          border: 'none',
          borderTop: '1px solid',
          borderColor: 'divider',
        },
        '& table': {
          width: '100%',
          mb: 2,
        },
        '& strong': {
          fontWeight: 600,
        },
      }}
    >
      <ReactMarkdown
        remarkPlugins={[remarkGfm]}
        components={{
          // Custom table rendering with MUI
          table: ({ children }) => (
            <TableContainer component={Paper} variant="outlined" sx={{ mb: 2 }}>
              <Table size="small">{children}</Table>
            </TableContainer>
          ),
          thead: ({ children }) => <TableHead>{children}</TableHead>,
          tbody: ({ children }) => <TableBody>{children}</TableBody>,
          tr: ({ children }) => <TableRow>{children}</TableRow>,
          th: ({ children }) => (
            <TableCell sx={{ fontWeight: 600, bgcolor: tableHeaderBg }}>
              {children}
            </TableCell>
          ),
          td: ({ children }) => <TableCell>{children}</TableCell>,
          // Custom heading rendering
          h2: ({ children }) => (
            <Typography variant="h6" component="h2" sx={{ mt: 3, mb: 1.5 }}>
              {children}
            </Typography>
          ),
          h3: ({ children }) => (
            <Typography variant="subtitle1" component="h3" sx={{ mt: 2, mb: 1, fontWeight: 600 }}>
              {children}
            </Typography>
          ),
          // Divider for hr
          hr: () => <Divider sx={{ my: 3 }} />,
        }}
      >
        {content}
      </ReactMarkdown>
    </Box>
  );
}

export default MarkdownDoc;
