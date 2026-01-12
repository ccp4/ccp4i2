'use client';

import { useState, useEffect, useMemo, Component, ReactNode } from 'react';
import {
  Box,
  Paper,
  Typography,
  CircularProgress,
  Alert,
  ToggleButtonGroup,
  ToggleButton,
  IconButton,
  Tooltip,
} from '@mui/material';
import { Download, Fullscreen, FullscreenExit, Code } from '@mui/icons-material';
import dynamic from 'next/dynamic';

// Error boundary to catch SeqViz loading/rendering errors
class SeqVizErrorBoundary extends Component<
  { children: ReactNode; fallback: ReactNode },
  { hasError: boolean }
> {
  constructor(props: { children: ReactNode; fallback: ReactNode }) {
    super(props);
    this.state = { hasError: false };
  }

  static getDerivedStateFromError() {
    return { hasError: true };
  }

  componentDidCatch(error: Error) {
    console.warn('SeqViz failed to load:', error.message);
  }

  render() {
    if (this.state.hasError) {
      return this.props.fallback;
    }
    return this.props.children;
  }
}

// Dynamic import of SeqViz with no SSR and error handling
const SeqVizComponent = dynamic(
  () =>
    import('seqviz')
      .then((mod) => mod.SeqViz)
      .catch((err) => {
        console.warn('SeqViz import failed:', err);
        // Return a placeholder component on import failure
        return () => null;
      }),
  {
    ssr: false,
    loading: () => (
      <Box sx={{ display: 'flex', justifyContent: 'center', alignItems: 'center', height: 400 }}>
        <CircularProgress />
      </Box>
    ),
  }
);

interface SeqVizViewerProps {
  /** GenBank file content as string */
  genbankContent: string | null;
  /** Loading state */
  loading?: boolean;
  /** Error message */
  error?: string;
  /** Default viewer type */
  defaultViewer?: 'linear' | 'circular' | 'both';
  /** Height of the viewer */
  height?: number;
  /** File name for download */
  filename?: string;
  /** URL for downloading the file */
  downloadUrl?: string;
}

type ViewerType = 'linear' | 'circular' | 'both' | 'raw';

/**
 * Fallback component that displays raw GenBank content
 */
function RawGenbankViewer({
  content,
  height,
}: {
  content: string;
  height: number;
}) {
  // Parse basic info from GenBank file
  const info = useMemo(() => {
    const lines = content.split('\n');
    const locusLine = lines.find((l) => l.startsWith('LOCUS'));
    const definitionLine = lines.find((l) => l.startsWith('DEFINITION'));
    const originIndex = lines.findIndex((l) => l.startsWith('ORIGIN'));

    let sequenceLength = 0;
    if (locusLine) {
      const match = locusLine.match(/(\d+)\s+bp/);
      if (match) sequenceLength = parseInt(match[1], 10);
    }

    return {
      locus: locusLine?.replace('LOCUS', '').trim().split(/\s+/)[0] || 'Unknown',
      definition: definitionLine?.replace('DEFINITION', '').trim() || '',
      length: sequenceLength,
      featureCount: (content.match(/^\s{5}\w+\s/gm) || []).length,
    };
  }, [content]);

  return (
    <Box sx={{ height, overflow: 'auto' }}>
      <Box sx={{ mb: 2, p: 2, bgcolor: 'grey.100', borderRadius: 1 }}>
        <Typography variant="subtitle2" gutterBottom>
          Sequence Information
        </Typography>
        <Typography variant="body2">
          <strong>Locus:</strong> {info.locus}
        </Typography>
        {info.definition && (
          <Typography variant="body2">
            <strong>Definition:</strong> {info.definition}
          </Typography>
        )}
        <Typography variant="body2">
          <strong>Length:</strong> {info.length.toLocaleString()} bp
        </Typography>
        <Typography variant="body2">
          <strong>Features:</strong> {info.featureCount}
        </Typography>
      </Box>
      <Box
        component="pre"
        sx={{
          p: 2,
          bgcolor: 'grey.50',
          borderRadius: 1,
          overflow: 'auto',
          fontSize: '0.75rem',
          fontFamily: 'monospace',
          whiteSpace: 'pre-wrap',
          wordBreak: 'break-all',
          maxHeight: height - 150,
        }}
      >
        {content}
      </Box>
    </Box>
  );
}

export function SeqVizViewer({
  genbankContent,
  loading = false,
  error,
  defaultViewer = 'both',
  height = 500,
  filename,
  downloadUrl,
}: SeqVizViewerProps) {
  const [viewerType, setViewerType] = useState<ViewerType>(defaultViewer);
  const [isFullscreen, setIsFullscreen] = useState(false);
  const [parseError, setParseError] = useState<string | null>(null);
  const [seqVizAvailable, setSeqVizAvailable] = useState(true);

  // Reset parse error when content changes
  useEffect(() => {
    setParseError(null);
  }, [genbankContent]);

  // Memoize the SeqViz props to prevent unnecessary re-renders
  const seqVizProps = useMemo(() => {
    if (!genbankContent || viewerType === 'raw') return null;

    return {
      file: genbankContent,
      viewer: viewerType as 'linear' | 'circular' | 'both',
      style: {
        height: isFullscreen ? 'calc(100vh - 100px)' : height,
        width: '100%',
      },
      showAnnotations: true,
      showPrimers: true,
      showIndex: true,
      translations: [],
      enzymes: [],
      primers: [],
      onError: (err: Error) => {
        console.error('SeqViz error:', err);
        setParseError(err.message);
      },
    };
  }, [genbankContent, viewerType, height, isFullscreen]);

  const handleViewerChange = (
    _: React.MouseEvent<HTMLElement>,
    newViewer: ViewerType | null
  ) => {
    if (newViewer !== null) {
      setViewerType(newViewer);
    }
  };

  const handleDownload = () => {
    if (downloadUrl) {
      window.open(downloadUrl, '_blank');
    } else if (genbankContent && filename) {
      const blob = new Blob([genbankContent], { type: 'text/plain' });
      const url = URL.createObjectURL(blob);
      const a = document.createElement('a');
      a.href = url;
      a.download = filename;
      a.click();
      URL.revokeObjectURL(url);
    }
  };

  // Fallback when SeqViz is not available
  const rawFallback = genbankContent ? (
    <RawGenbankViewer
      content={genbankContent}
      height={isFullscreen ? window.innerHeight - 100 : height}
    />
  ) : null;

  if (loading) {
    return (
      <Paper
        sx={{
          p: 4,
          display: 'flex',
          justifyContent: 'center',
          alignItems: 'center',
          height,
        }}
      >
        <CircularProgress />
      </Paper>
    );
  }

  if (error || parseError) {
    return (
      <Paper sx={{ p: 4, height }}>
        <Alert severity="error" sx={{ mb: 2 }}>
          {error || parseError}
        </Alert>
        {genbankContent && rawFallback}
      </Paper>
    );
  }

  if (!genbankContent) {
    return (
      <Paper
        sx={{
          p: 4,
          height,
          display: 'flex',
          justifyContent: 'center',
          alignItems: 'center',
        }}
      >
        <Typography color="text.secondary">
          No sequence file available
        </Typography>
      </Paper>
    );
  }

  return (
    <Paper
      sx={{
        p: 2,
        ...(isFullscreen && {
          position: 'fixed',
          top: 0,
          left: 0,
          right: 0,
          bottom: 0,
          zIndex: 1300,
          borderRadius: 0,
        }),
      }}
    >
      <Box
        sx={{
          mb: 2,
          display: 'flex',
          justifyContent: 'space-between',
          alignItems: 'center',
        }}
      >
        <Typography variant="subtitle2" color="text.secondary">
          Sequence Map
        </Typography>

        <Box sx={{ display: 'flex', alignItems: 'center', gap: 1 }}>
          <ToggleButtonGroup
            size="small"
            value={viewerType}
            exclusive
            onChange={handleViewerChange}
          >
            {seqVizAvailable && (
              <>
                <ToggleButton value="circular">Circular</ToggleButton>
                <ToggleButton value="linear">Linear</ToggleButton>
                <ToggleButton value="both">Both</ToggleButton>
              </>
            )}
            <ToggleButton value="raw">
              <Tooltip title="Raw GenBank">
                <Code fontSize="small" />
              </Tooltip>
            </ToggleButton>
          </ToggleButtonGroup>

          {(downloadUrl || (genbankContent && filename)) && (
            <Tooltip title="Download GenBank file">
              <IconButton onClick={handleDownload} size="small">
                <Download />
              </IconButton>
            </Tooltip>
          )}

          <Tooltip title={isFullscreen ? 'Exit fullscreen' : 'Fullscreen'}>
            <IconButton onClick={() => setIsFullscreen(!isFullscreen)} size="small">
              {isFullscreen ? <FullscreenExit /> : <Fullscreen />}
            </IconButton>
          </Tooltip>
        </Box>
      </Box>

      <Box sx={{ height: isFullscreen ? 'calc(100vh - 100px)' : height }}>
        {viewerType === 'raw' ? (
          rawFallback
        ) : (
          <SeqVizErrorBoundary
            fallback={
              <Box>
                <Alert severity="warning" sx={{ mb: 2 }}>
                  Interactive sequence viewer unavailable. Showing raw GenBank data.
                </Alert>
                {rawFallback}
              </Box>
            }
          >
            {seqVizProps && <SeqVizComponent {...seqVizProps} />}
          </SeqVizErrorBoundary>
        )}
      </Box>
    </Paper>
  );
}
