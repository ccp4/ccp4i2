'use client';

import React, { useEffect, useState } from 'react';
import {
  Container,
  Typography,
  Box,
  Paper,
  Table,
  TableBody,
  TableCell,
  TableContainer,
  TableHead,
  TableRow,
  Chip,
  Button,
  Alert,
  CircularProgress,
  IconButton,
  Tooltip,
  Menu,
  MenuItem,
  ListItemIcon,
  ListItemText,
  Divider,
} from '@mui/material';
import {
  AdminPanelSettings,
  PersonAdd,
  PersonRemove,
  Refresh,
  CloudUpload,
  MenuBook,
  DarkMode,
  LightMode,
  Visibility,
  Edit,
  MoreVert,
  Check,
  MailOutline,
  HourglassEmpty,
  ErrorOutline,
} from '@mui/icons-material';
import {
  Dialog,
  DialogTitle,
  DialogContent,
  DialogActions,
  TextField,
} from '@mui/material';
import Link from 'next/link';
import { apiDelete, apiGet, apiPost } from '../../lib/users/api';
import { useTheme } from '../../theme/theme-provider';

// Role types matching backend UserProfile.ROLE_CHOICES
type UserRole = 'user' | 'contributor' | 'admin';

const ROLE_LABELS: Record<UserRole, string> = {
  user: 'Viewer',
  contributor: 'Contributor',
  admin: 'Admin',
};

const ROLE_DESCRIPTIONS: Record<UserRole, string> = {
  user: 'Read-only access',
  contributor: 'Can add, edit, and delete data',
  admin: 'Full access including user management',
};

interface UserProfile {
  role: UserRole;
  is_platform_admin: boolean;
  legacy_username: string;
  legacy_display_name: string;
  imported_at: string | null;
  first_login_at: string | null;
  last_seen_at: string | null;
}

interface User {
  id: number;
  username: string;
  email: string;
  display_name: string;
  is_admin: boolean;
  is_active: boolean;
  first_name?: string;
  last_name?: string;
  profile?: UserProfile;
}

interface PendingInvite {
  id: number;
  email: string;
  note: string;
  status: 'pending' | 'sent' | 'failed';
  failure_reason: string;
  requested_by_email: string;
  requested_at: string;
  sent_by_email: string;
  sent_at: string | null;
  guest_object_id: string;
}

interface CurrentUser {
  id: number;
  username: string;
  email: string;
  first_name: string;
  last_name: string;
  display_name: string;
  is_admin: boolean;
  role: UserRole;
  operating_level: UserRole;
  can_contribute: boolean;
  can_administer: boolean;
  profile: UserProfile;
}

// Role icon component
const RoleIcon = ({ role }: { role: UserRole }) => {
  switch (role) {
    case 'user':
      return <Visibility fontSize="small" />;
    case 'contributor':
      return <Edit fontSize="small" />;
    case 'admin':
      return <AdminPanelSettings fontSize="small" />;
    default:
      return <Visibility fontSize="small" />;
  }
};

// Role color helper
const getRoleColor = (role: UserRole): 'default' | 'primary' | 'error' => {
  switch (role) {
    case 'user': return 'default';
    case 'contributor': return 'primary';
    case 'admin': return 'error';
    default: return 'default';
  }
};

export default function AdminPage() {
  const { mode, toggleTheme } = useTheme();
  const [currentUser, setCurrentUser] = useState<CurrentUser | null>(null);
  const [users, setUsers] = useState<User[]>([]);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);
  const [actionLoading, setActionLoading] = useState<number | null>(null);
  const [roleMenuAnchor, setRoleMenuAnchor] = useState<null | HTMLElement>(null);
  const [selectedUser, setSelectedUser] = useState<User | null>(null);
  const [invites, setInvites] = useState<PendingInvite[]>([]);
  const [inviteDialogOpen, setInviteDialogOpen] = useState(false);
  const [inviteEmail, setInviteEmail] = useState('');
  const [inviteNote, setInviteNote] = useState('');
  const [inviteSubmitting, setInviteSubmitting] = useState(false);

  const fetchCurrentUser = async (): Promise<CurrentUser | null> => {
    try {
      return await apiGet<CurrentUser>('me');
    } catch (err) {
      console.error('Error fetching current user:', err);
      return null;
    }
  };

  const fetchUsers = async (): Promise<User[]> => {
    try {
      return await apiGet<User[]>('users');
    } catch (err: any) {
      if (err.status === 403) {
        throw new Error('You do not have permission to view users. Admin access required.');
      }
      throw new Error('Failed to fetch users');
    }
  };

  const fetchInvites = async (): Promise<PendingInvite[]> => {
    try {
      return await apiGet<PendingInvite[]>('pending-invites');
    } catch {
      return [];
    }
  };

  const loadData = async () => {
    setLoading(true);
    setError(null);
    try {
      const user = await fetchCurrentUser();
      setCurrentUser(user);

      if (user?.is_admin) {
        const [userList, inviteList] = await Promise.all([
          fetchUsers(),
          fetchInvites(),
        ]);
        setUsers(userList);
        setInvites(inviteList);
      }
    } catch (err) {
      setError(err instanceof Error ? err.message : 'An error occurred');
    } finally {
      setLoading(false);
    }
  };

  const handleQueueInvite = async () => {
    const email = inviteEmail.trim();
    if (!email) return;
    setInviteSubmitting(true);
    try {
      await apiPost('pending-invites', { email, note: inviteNote.trim() });
      setInviteEmail('');
      setInviteNote('');
      setInviteDialogOpen(false);
      setInvites(await fetchInvites());
    } catch (err) {
      setError(err instanceof Error ? err.message : 'Failed to queue invite');
    } finally {
      setInviteSubmitting(false);
    }
  };

  const handleCancelInvite = async (inviteId: number) => {
    try {
      await apiDelete(`pending-invites/${inviteId}`);
      setInvites(await fetchInvites());
    } catch (err) {
      setError(err instanceof Error ? err.message : 'Failed to cancel invite');
    }
  };

  useEffect(() => {
    loadData();
  }, []);

  const handleGrantAdmin = async (userId: number) => {
    setActionLoading(userId);
    try {
      await apiPost(`users/${userId}/grant_admin`);
      await loadData();
    } catch (err) {
      setError(err instanceof Error ? err.message : 'Failed to grant admin');
    } finally {
      setActionLoading(null);
    }
  };

  const handleRevokeAdmin = async (userId: number) => {
    setActionLoading(userId);
    try {
      await apiPost(`users/${userId}/revoke_admin`);
      await loadData();
    } catch (err) {
      setError(err instanceof Error ? err.message : 'Failed to revoke admin');
    } finally {
      setActionLoading(null);
    }
  };

  // New role-based handlers
  const handleRoleMenuOpen = (event: React.MouseEvent<HTMLElement>, user: User) => {
    setRoleMenuAnchor(event.currentTarget);
    setSelectedUser(user);
  };

  const handleRoleMenuClose = () => {
    setRoleMenuAnchor(null);
    setSelectedUser(null);
  };

  const handleSetRole = async (newRole: UserRole) => {
    if (!selectedUser) return;

    setActionLoading(selectedUser.id);
    handleRoleMenuClose();

    try {
      await apiPost(`users/${selectedUser.id}/set_role`, { role: newRole });
      await loadData();
    } catch (err) {
      setError(err instanceof Error ? err.message : 'Failed to update role');
    } finally {
      setActionLoading(null);
    }
  };

  if (loading) {
    return (
      <Container maxWidth="lg" sx={{ py: 4 }}>
        <Box sx={{ display: 'flex', justifyContent: 'center', py: 8 }}>
          <CircularProgress />
        </Box>
      </Container>
    );
  }

  return (
    <Container maxWidth="lg" sx={{ py: 4 }}>
      <Box sx={{ display: 'flex', alignItems: 'center', gap: 2, mb: 4 }}>
        <AdminPanelSettings sx={{ fontSize: 40, color: 'grey.600' }} />
        <Box>
          <Typography variant="h4" component="h1">
            Platform Admin
          </Typography>
          <Typography color="text.secondary">
            User management and platform settings
          </Typography>
        </Box>
        <Box sx={{ ml: 'auto', display: 'flex', alignItems: 'center', gap: 1 }}>
          <Tooltip title={mode === 'light' ? 'Switch to dark mode' : 'Switch to light mode'}>
            <IconButton onClick={toggleTheme} sx={{ color: 'text.secondary' }}>
              {mode === 'light' ? <DarkMode /> : <LightMode />}
            </IconButton>
          </Tooltip>
          <Button component={Link} href="/" variant="outlined">
            Back to Apps
          </Button>
        </Box>
      </Box>

      {error && (
        <Alert severity="error" sx={{ mb: 3 }} onClose={() => setError(null)}>
          {error}
        </Alert>
      )}

      {/* Current User Info */}
      <Paper sx={{ p: 3, mb: 3 }}>
        <Typography variant="h6" gutterBottom>
          Current User
        </Typography>
        {currentUser ? (
          <Box sx={{ display: 'flex', alignItems: 'center', gap: 2 }}>
            <Typography>
              <strong>{currentUser.display_name}</strong> ({currentUser.email || currentUser.username})
            </Typography>
            {currentUser.is_admin && (
              <Chip label="Platform Admin" color="primary" size="small" />
            )}
          </Box>
        ) : (
          <Typography color="text.secondary">Not authenticated</Typography>
        )}
      </Paper>

      {/* Quick Links for Admin */}
      {currentUser?.is_admin && (
        <Paper sx={{ p: 3, mb: 3 }}>
          <Typography variant="h6" gutterBottom>
            Admin Tools
          </Typography>
          <Box sx={{ display: 'flex', gap: 2, flexWrap: 'wrap' }}>
            <Button
              component={Link}
              href="/admin/import"
              variant="outlined"
              startIcon={<CloudUpload />}
            >
              Data Import
            </Button>
            <Button
              component={Link}
              href="/admin/docs"
              variant="outlined"
              startIcon={<MenuBook />}
            >
              Documentation
            </Button>
          </Box>
        </Paper>
      )}

      {/* Invite queue */}
      {currentUser?.is_admin && (
        <Paper sx={{ p: 3, mb: 3 }}>
          <Box sx={{ display: 'flex', alignItems: 'center', mb: 2, gap: 1 }}>
            <Typography variant="h6">Invite Queue</Typography>
            <Tooltip title="Refresh">
              <IconButton onClick={loadData} size="small">
                <Refresh fontSize="small" />
              </IconButton>
            </Tooltip>
            <Button
              sx={{ ml: 'auto' }}
              variant="contained"
              startIcon={<PersonAdd />}
              onClick={() => setInviteDialogOpen(true)}
            >
              Invite User
            </Button>
          </Box>

          <Typography variant="body2" color="text.secondary" sx={{ mb: 2 }}>
            Queue external collaborators here. A platform operator drains the queue
            out-of-band via <code>invite-user.sh</code> (see operations docs) — the
            deployed app is not permitted to send B2B invitations directly.
          </Typography>

          <TableContainer>
            <Table size="small">
              <TableHead>
                <TableRow>
                  <TableCell>Email</TableCell>
                  <TableCell>Note</TableCell>
                  <TableCell>Status</TableCell>
                  <TableCell>Requested</TableCell>
                  <TableCell>Sent</TableCell>
                  <TableCell align="right">Actions</TableCell>
                </TableRow>
              </TableHead>
              <TableBody>
                {invites.map((inv) => {
                  const statusChip = inv.status === 'pending' ? (
                    <Chip icon={<HourglassEmpty fontSize="small" />} label="Pending" size="small" color="warning" />
                  ) : inv.status === 'sent' ? (
                    <Chip icon={<MailOutline fontSize="small" />} label="Sent" size="small" color="success" />
                  ) : (
                    <Tooltip title={inv.failure_reason || 'Failed'}>
                      <Chip icon={<ErrorOutline fontSize="small" />} label="Failed" size="small" color="error" />
                    </Tooltip>
                  );
                  return (
                    <TableRow key={inv.id}>
                      <TableCell>{inv.email}</TableCell>
                      <TableCell>{inv.note}</TableCell>
                      <TableCell>{statusChip}</TableCell>
                      <TableCell>
                        {new Date(inv.requested_at).toLocaleDateString()}
                        {inv.requested_by_email && (
                          <Typography variant="caption" display="block" color="text.secondary">
                            by {inv.requested_by_email}
                          </Typography>
                        )}
                      </TableCell>
                      <TableCell>
                        {inv.sent_at ? new Date(inv.sent_at).toLocaleDateString() : '—'}
                      </TableCell>
                      <TableCell align="right">
                        {inv.status === 'pending' && (
                          <Tooltip title="Cancel">
                            <IconButton size="small" onClick={() => handleCancelInvite(inv.id)}>
                              <PersonRemove fontSize="small" />
                            </IconButton>
                          </Tooltip>
                        )}
                      </TableCell>
                    </TableRow>
                  );
                })}
                {invites.length === 0 && (
                  <TableRow>
                    <TableCell colSpan={6} align="center">
                      <Typography color="text.secondary" variant="body2">
                        No invites queued
                      </Typography>
                    </TableCell>
                  </TableRow>
                )}
              </TableBody>
            </Table>
          </TableContainer>
        </Paper>
      )}

      {/* User Management */}
      {currentUser?.is_admin ? (
        <Paper sx={{ p: 3 }}>
          <Box sx={{ display: 'flex', alignItems: 'center', mb: 2 }}>
            <Typography variant="h6">User Management</Typography>
            <Tooltip title="Refresh">
              <IconButton onClick={loadData} sx={{ ml: 1 }}>
                <Refresh />
              </IconButton>
            </Tooltip>
          </Box>

          <TableContainer>
            <Table>
              <TableHead>
                <TableRow>
                  <TableCell>Name</TableCell>
                  <TableCell>Email</TableCell>
                  <TableCell>Status</TableCell>
                  <TableCell>Role</TableCell>
                  <TableCell>Last Seen</TableCell>
                  <TableCell align="right">Actions</TableCell>
                </TableRow>
              </TableHead>
              <TableBody>
                {users.map((user) => (
                  <TableRow key={user.id}>
                    <TableCell>{user.display_name || user.username}</TableCell>
                    <TableCell>{user.email}</TableCell>
                    <TableCell>
                      <Chip
                        label={user.is_active ? 'Active' : 'Inactive'}
                        color={user.is_active ? 'success' : 'default'}
                        size="small"
                      />
                    </TableCell>
                    <TableCell>
                      {(() => {
                        const userRole = user.profile?.role || (user.is_admin ? 'admin' : 'contributor');
                        return (
                          <Chip
                            icon={<RoleIcon role={userRole} />}
                            label={ROLE_LABELS[userRole]}
                            color={getRoleColor(userRole)}
                            variant={userRole === 'user' ? 'outlined' : 'filled'}
                            size="small"
                          />
                        );
                      })()}
                    </TableCell>
                    <TableCell>
                      {user.profile?.last_seen_at
                        ? new Date(user.profile.last_seen_at).toLocaleDateString()
                        : 'Never'}
                    </TableCell>
                    <TableCell align="right">
                      {actionLoading === user.id ? (
                        <CircularProgress size={24} />
                      ) : user.id === currentUser?.id ? (
                        <Tooltip title="Cannot modify your own role">
                          <span>
                            <IconButton size="small" disabled>
                              <MoreVert fontSize="small" />
                            </IconButton>
                          </span>
                        </Tooltip>
                      ) : (
                        <Tooltip title="Change role">
                          <IconButton
                            size="small"
                            onClick={(e) => handleRoleMenuOpen(e, user)}
                          >
                            <MoreVert fontSize="small" />
                          </IconButton>
                        </Tooltip>
                      )}
                    </TableCell>
                  </TableRow>
                ))}
                {users.length === 0 && (
                  <TableRow>
                    <TableCell colSpan={6} align="center">
                      <Typography color="text.secondary">No users found</Typography>
                    </TableCell>
                  </TableRow>
                )}
              </TableBody>
            </Table>
          </TableContainer>
        </Paper>
      ) : (
        <Alert severity="warning">
          You need platform admin access to manage users. Contact an administrator if you need access.
        </Alert>
      )}

      {/* Invite dialog */}
      <Dialog open={inviteDialogOpen} onClose={() => setInviteDialogOpen(false)} maxWidth="sm" fullWidth>
        <DialogTitle>Queue external invite</DialogTitle>
        <DialogContent>
          <Typography variant="body2" color="text.secondary" sx={{ mb: 2 }}>
            The invite is not sent immediately. It waits in the queue until a
            platform operator drains it.
          </Typography>
          <TextField
            autoFocus
            margin="dense"
            label="Email"
            type="email"
            fullWidth
            value={inviteEmail}
            onChange={(e) => setInviteEmail(e.target.value)}
            disabled={inviteSubmitting}
          />
          <TextField
            margin="dense"
            label="Note (optional)"
            fullWidth
            value={inviteNote}
            onChange={(e) => setInviteNote(e.target.value)}
            placeholder="e.g. Jane Smith, Diamond Light Source, collaborator on target X"
            disabled={inviteSubmitting}
          />
        </DialogContent>
        <DialogActions>
          <Button onClick={() => setInviteDialogOpen(false)} disabled={inviteSubmitting}>
            Cancel
          </Button>
          <Button
            variant="contained"
            onClick={handleQueueInvite}
            disabled={inviteSubmitting || !inviteEmail.trim()}
          >
            {inviteSubmitting ? <CircularProgress size={20} /> : 'Queue invite'}
          </Button>
        </DialogActions>
      </Dialog>

      {/* Role change menu */}
      <Menu
        anchorEl={roleMenuAnchor}
        open={Boolean(roleMenuAnchor)}
        onClose={handleRoleMenuClose}
        anchorOrigin={{ vertical: 'bottom', horizontal: 'right' }}
        transformOrigin={{ vertical: 'top', horizontal: 'right' }}
      >
        <MenuItem disabled sx={{ opacity: 1 }}>
          <Typography variant="caption" color="text.secondary">
            Set role for {selectedUser?.display_name || selectedUser?.username}
          </Typography>
        </MenuItem>
        <Divider />
        {(['user', 'contributor', 'admin'] as UserRole[]).map((role) => {
          const currentRole = selectedUser?.profile?.role ||
            (selectedUser?.is_admin ? 'admin' : 'contributor');
          const isSelected = role === currentRole;

          return (
            <MenuItem
              key={role}
              onClick={() => handleSetRole(role)}
              selected={isSelected}
            >
              <ListItemIcon>
                <RoleIcon role={role} />
              </ListItemIcon>
              <ListItemText
                primary={ROLE_LABELS[role]}
                secondary={ROLE_DESCRIPTIONS[role]}
                secondaryTypographyProps={{ variant: 'caption' }}
              />
              {isSelected && <Check fontSize="small" color="primary" sx={{ ml: 1 }} />}
            </MenuItem>
          );
        })}
      </Menu>
    </Container>
  );
}
