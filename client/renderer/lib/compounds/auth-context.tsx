'use client';

/**
 * Stub auth context for local development (Electron/desktop mode).
 *
 * In Docker builds, this file is REPLACED by the full auth-context.tsx from
 * apps/compounds/frontend/lib/compounds/auth-context.tsx which fetches
 * user info from Django's /api/users/me/ endpoint.
 *
 * This stub provides a mock admin user with full access for local development.
 */

import React, { createContext, useContext, ReactNode } from 'react';

// Types matching the compounds auth-context
export type UserRole = 'user' | 'contributor' | 'admin';

export interface CurrentUser {
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
  profile: {
    role: UserRole;
    is_platform_admin: boolean;
    legacy_username: string;
    legacy_display_name: string;
    imported_at: string | null;
    first_login_at: string | null;
    last_seen_at: string | null;
  };
}

interface AuthContextValue {
  user: CurrentUser | null;
  isLoading: boolean;
  error: Error | null;
  operatingLevel: UserRole;
  availableLevels: UserRole[];
  canContribute: boolean;
  canAdminister: boolean;
  setOperatingLevel: (level: UserRole) => Promise<void>;
  refreshUser: () => Promise<void>;
}

// Mock user with full admin access for local development
const mockUser: CurrentUser = {
  id: 0,
  username: 'local_user',
  email: '',
  first_name: 'Local',
  last_name: 'User',
  display_name: 'Local User',
  is_admin: true,
  role: 'admin',
  operating_level: 'admin',
  can_contribute: true,
  can_administer: true,
  profile: {
    role: 'admin',
    is_platform_admin: true,
    legacy_username: '',
    legacy_display_name: '',
    imported_at: null,
    first_login_at: null,
    last_seen_at: null,
  },
};

const defaultContext: AuthContextValue = {
  user: mockUser,
  isLoading: false,
  error: null,
  operatingLevel: 'admin',
  availableLevels: ['user', 'contributor', 'admin'],
  canContribute: true,
  canAdminister: true,
  setOperatingLevel: async () => {},
  refreshUser: async () => {},
};

const AuthContext = createContext<AuthContextValue>(defaultContext);

interface AuthProviderProps {
  children: ReactNode;
}

/**
 * Stub AuthProvider for local development.
 * Provides a mock admin user - no API calls are made.
 */
export function AuthProvider({ children }: AuthProviderProps) {
  return (
    <AuthContext.Provider value={defaultContext}>
      {children}
    </AuthContext.Provider>
  );
}

/**
 * Hook to access auth context.
 */
export function useAuth(): AuthContextValue {
  const context = useContext(AuthContext);
  if (!context) {
    throw new Error('useAuth must be used within an AuthProvider');
  }
  return context;
}

/**
 * Human-readable labels for role levels
 */
export const ROLE_LABELS: Record<UserRole, string> = {
  user: 'Viewer',
  contributor: 'Contributor',
  admin: 'Admin',
};

/**
 * Descriptions for each role level
 */
export const ROLE_DESCRIPTIONS: Record<UserRole, string> = {
  user: 'Read-only access',
  contributor: 'Can add, edit, and delete data',
  admin: 'Full access including user management',
};

/**
 * Icon names for each role (MUI icon names)
 */
export const ROLE_ICONS: Record<UserRole, 'Visibility' | 'Edit' | 'AdminPanelSettings'> = {
  user: 'Visibility',
  contributor: 'Edit',
  admin: 'AdminPanelSettings',
};
