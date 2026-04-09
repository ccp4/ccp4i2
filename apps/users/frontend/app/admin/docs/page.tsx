'use client';

import React, { useEffect, useState } from 'react';
import {
  Container,
  Typography,
  Box,
  Paper,
  Button,
  Alert,
  CircularProgress,
  IconButton,
  Tooltip,
} from '@mui/material';
import {
  MenuBook,
  ArrowBack,
  DarkMode,
  LightMode,
} from '@mui/icons-material';
import Link from 'next/link';
import dynamic from 'next/dynamic';
import { apiGet } from '../../../lib/users/api';
import { useTheme } from '../../../theme/theme-provider';

// Dynamic import of MarkdownDoc to avoid SSR issues with react-markdown
const MarkdownDoc = dynamic(
  () => import('../../../components/compounds/docs/MarkdownDoc').then(mod => mod.MarkdownDoc),
  { ssr: false, loading: () => <CircularProgress size={24} /> }
);

interface CurrentUser {
  id: number;
  username: string;
  email: string;
  first_name: string;
  last_name: string;
  display_name: string;
  is_admin: boolean;
}

export default function DocsPage() {
  const { mode, toggleTheme } = useTheme();
  const [currentUser, setCurrentUser] = useState<CurrentUser | null>(null);
  const [loading, setLoading] = useState(true);

  const fetchCurrentUser = async (): Promise<CurrentUser | null> => {
    try {
      return await apiGet<CurrentUser>('me');
    } catch (err) {
      console.error('Error fetching current user:', err);
      return null;
    }
  };

  useEffect(() => {
    const loadData = async () => {
      setLoading(true);
      const user = await fetchCurrentUser();
      setCurrentUser(user);
      setLoading(false);
    };
    loadData();
  }, []);

  if (loading) {
    return (
      <Container maxWidth="lg" sx={{ py: 4 }}>
        <Box sx={{ display: 'flex', justifyContent: 'center', py: 8 }}>
          <CircularProgress />
        </Box>
      </Container>
    );
  }

  if (!currentUser?.is_admin) {
    return (
      <Container maxWidth="lg" sx={{ py: 4 }}>
        <Alert severity="warning">
          You need platform admin access to view admin documentation. Contact an administrator if you need access.
        </Alert>
        <Box sx={{ mt: 2 }}>
          <Button component={Link} href="/admin" startIcon={<ArrowBack />}>
            Back to Admin
          </Button>
        </Box>
      </Container>
    );
  }

  return (
    <Container maxWidth="lg" sx={{ py: 4 }}>
      <Box sx={{ display: 'flex', alignItems: 'center', gap: 2, mb: 4 }}>
        <MenuBook sx={{ fontSize: 40, color: 'grey.600' }} />
        <Box>
          <Typography variant="h4" component="h1">
            Admin Documentation
          </Typography>
          <Typography color="text.secondary">
            Platform administration guides and reference
          </Typography>
        </Box>
        <Box sx={{ ml: 'auto', display: 'flex', alignItems: 'center', gap: 1 }}>
          <Tooltip title={mode === 'light' ? 'Switch to dark mode' : 'Switch to light mode'}>
            <IconButton onClick={toggleTheme} sx={{ color: 'text.secondary' }}>
              {mode === 'light' ? <DarkMode /> : <LightMode />}
            </IconButton>
          </Tooltip>
          <Button component={Link} href="/admin" variant="outlined" startIcon={<ArrowBack />}>
            Back to Admin
          </Button>
        </Box>
      </Box>

      {/* Admin Guide Documentation */}
      <Paper sx={{ p: 3 }}>
        <MarkdownDoc src="admin-guide.md" />
      </Paper>
    </Container>
  );
}
