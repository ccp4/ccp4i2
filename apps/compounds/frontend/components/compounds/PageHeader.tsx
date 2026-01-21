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
  Select,
  SelectChangeEvent,
  Typography,
} from '@mui/material';
import {
  Home,
  Science,
  Logout,
  Person,
  DarkMode,
  LightMode,
  Visibility,
  Edit,
  AdminPanelSettings,
  Check,
} from '@mui/icons-material';
import { useRouter } from 'next/navigation';
import { Breadcrumbs, BreadcrumbItem } from './Breadcrumbs';
import { routes } from '@/lib/compounds/routes';
import { useTheme } from '@/lib/compounds/theme-provider';
import { useAuth, ROLE_LABELS, ROLE_DESCRIPTIONS } from '@/lib/compounds/auth-context';
import { UserRole } from '@/types/compounds/models';

interface PageHeaderProps {
  /** Breadcrumb items for the current page */
  breadcrumbs: BreadcrumbItem[];
  /** Hide the navigation actions (home, targets, user menu) */
  hideActions?: boolean;
}

const REQUIRE_AUTH = process.env.NEXT_PUBLIC_REQUIRE_AUTH === 'true';

// Icons for each role level
const RoleIcon = ({ role, fontSize = 'small' }: { role: UserRole; fontSize?: 'small' | 'inherit' }) => {
  switch (role) {
    case 'user':
      return <Visibility fontSize={fontSize} />;
    case 'contributor':
      return <Edit fontSize={fontSize} />;
    case 'admin':
      return <AdminPanelSettings fontSize={fontSize} />;
    default:
      return <Person fontSize={fontSize} />;
  }
};

// Color for each role level
const getRoleColor = (role: UserRole): 'default' | 'primary' | 'secondary' | 'error' | 'info' | 'success' | 'warning' => {
  switch (role) {
    case 'user':
      return 'default';
    case 'contributor':
      return 'primary';
    case 'admin':
      return 'error';
    default:
      return 'default';
  }
};

/**
 * Consistent page header for compounds app pages.
 * Includes breadcrumbs and navigation actions (home, targets, user menu with logout).
 *
 * Features:
 * - Operating level selector (user/contributor/admin) for role-based access control
 * - User menu with logout (when auth is enabled)
 * - Theme toggle (light/dark mode)
 * - Navigation shortcuts (home, targets)
 *
 * Note: User authentication features are only available when
 * running in the Docker/ccp4i2 environment with auth configured.
 */
export function PageHeader({ breadcrumbs, hideActions = false }: PageHeaderProps) {
  const router = useRouter();
  const { mode, toggleTheme } = useTheme();
  const {
    user,
    operatingLevel,
    availableLevels,
    canContribute,
    setOperatingLevel,
    isLoading,
  } = useAuth();

  const [userMenuAnchor, setUserMenuAnchor] = useState<null | HTMLElement>(null);
  const [roleMenuAnchor, setRoleMenuAnchor] = useState<null | HTMLElement>(null);
  const [isChangingRole, setIsChangingRole] = useState(false);

  const handleLogout = () => {
    // Handled by parent app in production
    handleUserMenuClose();
    // In integrated mode, this would call the MSAL logout
  };

  const handleUserMenuOpen = (event: React.MouseEvent<HTMLElement>) => {
    setUserMenuAnchor(event.currentTarget);
  };

  const handleUserMenuClose = () => {
    setUserMenuAnchor(null);
  };

  const handleRoleMenuOpen = (event: React.MouseEvent<HTMLElement>) => {
    setRoleMenuAnchor(event.currentTarget);
  };

  const handleRoleMenuClose = () => {
    setRoleMenuAnchor(null);
  };

  const handleRoleChange = async (newRole: UserRole) => {
    if (newRole === operatingLevel) {
      handleRoleMenuClose();
      return;
    }

    setIsChangingRole(true);
    try {
      await setOperatingLevel(newRole);
    } catch (err) {
      console.error('Failed to change operating level:', err);
    } finally {
      setIsChangingRole(false);
      handleRoleMenuClose();
    }
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
            {/* Operating Level Selector - Show when auth is required */}
            {REQUIRE_AUTH && user && availableLevels.length > 1 && (
              <>
                <Tooltip
                  title={
                    <Box>
                      <Typography variant="body2" fontWeight="bold">
                        Operating as: {ROLE_LABELS[operatingLevel]}
                      </Typography>
                      <Typography variant="caption">
                        {ROLE_DESCRIPTIONS[operatingLevel]}
                      </Typography>
                      <Typography variant="caption" display="block" sx={{ mt: 0.5, fontStyle: 'italic' }}>
                        Click to change
                      </Typography>
                    </Box>
                  }
                >
                  <Chip
                    icon={<RoleIcon role={operatingLevel} />}
                    label={ROLE_LABELS[operatingLevel]}
                    size="small"
                    color={getRoleColor(operatingLevel)}
                    variant={operatingLevel === 'user' ? 'outlined' : 'filled'}
                    onClick={handleRoleMenuOpen}
                    disabled={isChangingRole}
                    sx={{
                      cursor: 'pointer',
                      fontWeight: 500,
                      '&:hover': { opacity: 0.9 },
                    }}
                  />
                </Tooltip>
                <Menu
                  anchorEl={roleMenuAnchor}
                  open={Boolean(roleMenuAnchor)}
                  onClose={handleRoleMenuClose}
                  anchorOrigin={{ vertical: 'bottom', horizontal: 'right' }}
                  transformOrigin={{ vertical: 'top', horizontal: 'right' }}
                >
                  <MenuItem disabled sx={{ opacity: 1 }}>
                    <Typography variant="caption" color="text.secondary">
                      Operating Level
                    </Typography>
                  </MenuItem>
                  <Divider />
                  {(['user', 'contributor', 'admin'] as UserRole[])
                    .filter(role => availableLevels.includes(role))
                    .map(role => (
                      <MenuItem
                        key={role}
                        onClick={() => handleRoleChange(role)}
                        selected={role === operatingLevel}
                      >
                        <ListItemIcon>
                          <RoleIcon role={role} />
                        </ListItemIcon>
                        <ListItemText
                          primary={ROLE_LABELS[role]}
                          secondary={ROLE_DESCRIPTIONS[role]}
                          secondaryTypographyProps={{ variant: 'caption' }}
                        />
                        {role === operatingLevel && (
                          <Check fontSize="small" color="primary" sx={{ ml: 1 }} />
                        )}
                      </MenuItem>
                    ))}
                </Menu>
                <Divider orientation="vertical" flexItem sx={{ mx: 0.5 }} />
              </>
            )}

            {/* Read-only indicator when operating as user */}
            {REQUIRE_AUTH && user && !canContribute && (
              <>
                <Tooltip title="You are in read-only mode. Switch to Contributor or Admin to make changes.">
                  <Chip
                    icon={<Visibility fontSize="small" />}
                    label="Read-only"
                    size="small"
                    variant="outlined"
                    color="default"
                    sx={{ fontStyle: 'italic' }}
                  />
                </Tooltip>
                <Divider orientation="vertical" flexItem sx={{ mx: 0.5 }} />
              </>
            )}

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

            {/* Theme toggle */}
            <Tooltip title={mode === 'light' ? 'Switch to dark mode' : 'Switch to light mode'}>
              <IconButton
                size="small"
                onClick={toggleTheme}
                sx={{ color: 'text.secondary' }}
              >
                {mode === 'light' ? <DarkMode fontSize="small" /> : <LightMode fontSize="small" />}
              </IconButton>
            </Tooltip>

            {/* User menu - only show when auth is required and user is available */}
            {REQUIRE_AUTH && user && (
              <>
                <Divider orientation="vertical" flexItem sx={{ mx: 0.5 }} />
                <Tooltip title={user.email || user.username || 'User'}>
                  <Chip
                    icon={<Person fontSize="small" />}
                    label={user.display_name?.split(' ')[0] || 'User'}
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
                  <MenuItem disabled sx={{ opacity: 1 }}>
                    <ListItemText
                      primary={user.display_name}
                      secondary={
                        <Box component="span">
                          <Typography variant="caption" display="block">
                            {user.email}
                          </Typography>
                          <Typography variant="caption" display="block" color="text.secondary">
                            Role: {ROLE_LABELS[user.role]}
                          </Typography>
                        </Box>
                      }
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
