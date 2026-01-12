'use client';

import { useState } from 'react';
import {
  Box,
  IconButton,
  Menu,
  MenuItem,
  ListItemIcon,
  ListItemText,
  Divider,
  Tooltip,
  Chip,
} from '@mui/material';
import {
  Home,
  Science,
  Logout,
  Person,
} from '@mui/icons-material';
import { useRouter } from 'next/navigation';
import { Breadcrumbs, BreadcrumbItem } from './Breadcrumbs';
import { routes } from '@/lib/compounds/routes';

interface PageHeaderProps {
  /** Breadcrumb items for the current page */
  breadcrumbs: BreadcrumbItem[];
  /** Hide the navigation actions (home, targets, user menu) */
  hideActions?: boolean;
}

const REQUIRE_AUTH = process.env.NEXT_PUBLIC_REQUIRE_AUTH === 'true';

/**
 * Consistent page header for compounds app pages.
 * Includes breadcrumbs and navigation actions (home, targets, user menu with logout).
 *
 * Note: User authentication features (user menu, logout) are only available when
 * running in the Docker/ccp4i2 environment with MSAL configured. In standalone
 * development mode, auth features are disabled.
 */
export function PageHeader({ breadcrumbs, hideActions = false }: PageHeaderProps) {
  const router = useRouter();
  const [userMenuAnchor, setUserMenuAnchor] = useState<null | HTMLElement>(null);

  // Auth context is injected by the parent app (ccp4i2) when running in Docker
  // In standalone mode, we don't have access to MSAL
  const currentUser = null; // Populated by parent app context in production
  const handleLogout = () => {
    // Handled by parent app in production
    handleUserMenuClose();
  };

  const handleUserMenuOpen = (event: React.MouseEvent<HTMLElement>) => {
    setUserMenuAnchor(event.currentTarget);
  };

  const handleUserMenuClose = () => {
    setUserMenuAnchor(null);
  };

  const handleNavigateHome = () => {
    router.push('/');
  };

  const handleNavigateTargets = () => {
    router.push(routes.registry.targets());
  };

  return (
    <Box sx={{ mb: 2 }}>
      {/* Main row with breadcrumbs and actions */}
      <Box
        sx={{
          display: 'flex',
          alignItems: 'center',
          justifyContent: 'space-between',
          flexWrap: 'wrap',
          gap: 1,
        }}
      >
        {/* Breadcrumbs */}
        <Box sx={{ flexGrow: 1 }}>
          <Breadcrumbs items={breadcrumbs} />
        </Box>

        {/* Navigation actions */}
        {!hideActions && (
          <Box
            sx={{
              display: 'flex',
              alignItems: 'center',
              gap: 0.5,
            }}
          >
            {/* Home / App Selector */}
            <Tooltip title={REQUIRE_AUTH ? "Home / App Selector" : "Home"}>
              <IconButton
                size="small"
                onClick={handleNavigateHome}
                sx={{ color: 'text.secondary' }}
              >
                <Home fontSize="small" />
              </IconButton>
            </Tooltip>

            {/* Targets */}
            <Tooltip title="All Targets">
              <IconButton
                size="small"
                onClick={handleNavigateTargets}
                sx={{ color: 'text.secondary' }}
              >
                <Science fontSize="small" />
              </IconButton>
            </Tooltip>

            {/* User menu - only show when auth is required and user is available */}
            {REQUIRE_AUTH && currentUser && (
              <>
                <Divider orientation="vertical" flexItem sx={{ mx: 0.5 }} />
                <Tooltip title={(currentUser as { username?: string }).username || 'User'}>
                  <Chip
                    icon={<Person fontSize="small" />}
                    label={(currentUser as { name?: string }).name?.split(' ')[0] || 'User'}
                    size="small"
                    variant="outlined"
                    onClick={handleUserMenuOpen}
                    sx={{
                      cursor: 'pointer',
                      '&:hover': { bgcolor: 'action.hover' },
                    }}
                  />
                </Tooltip>
                <Menu
                  anchorEl={userMenuAnchor}
                  open={Boolean(userMenuAnchor)}
                  onClose={handleUserMenuClose}
                  anchorOrigin={{ vertical: 'bottom', horizontal: 'right' }}
                  transformOrigin={{ vertical: 'top', horizontal: 'right' }}
                >
                  <MenuItem disabled>
                    <ListItemText
                      primary={(currentUser as { name?: string }).name}
                      secondary={(currentUser as { username?: string }).username}
                      primaryTypographyProps={{ fontWeight: 500 }}
                    />
                  </MenuItem>
                  <Divider />
                  <MenuItem onClick={handleLogout}>
                    <ListItemIcon>
                      <Logout fontSize="small" />
                    </ListItemIcon>
                    <ListItemText>Sign Out</ListItemText>
                  </MenuItem>
                </Menu>
              </>
            )}
          </Box>
        )}
      </Box>
    </Box>
  );
}
