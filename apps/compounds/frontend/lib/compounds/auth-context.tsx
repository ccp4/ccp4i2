'use client';

/**
 * Authentication context for the compounds frontend.
 *
 * Provides:
 * - Current user info (role, operating level, display name)
 * - Permission helpers (canContribute, canAdminister)
 * - Operating level management (switching between user/contributor/admin modes)
 *
 * In standalone development mode (NEXT_PUBLIC_REQUIRE_AUTH !== 'true'),
 * this provides a mock user with full admin access.
 */

import React, { createContext, useContext, useState, useCallback, useEffect } from 'react';
import useSWR from 'swr';
import { CurrentUser, UserRole, OperatingLevelInfo } from '@/types/compounds/models';

// =============================================================================
// Configuration
// =============================================================================

const REQUIRE_AUTH = process.env.NEXT_PUBLIC_REQUIRE_AUTH === 'true';

// API endpoints - uses /api/proxy/users to match proxy pattern
const USERS_API_BASE = '/api/proxy/users';

// =============================================================================
// Types
// =============================================================================

interface AuthContextValue {
  /** Current user info (null while loading) */
  user: CurrentUser | null;
  /** Whether auth data is still loading */
  isLoading: boolean;
  /** Error if auth fetch failed */
  error: Error | null;
  /** Current operating level (may be <= user's role) */
  operatingLevel: UserRole;
  /** Available operating levels for this user */
  availableLevels: UserRole[];
  /** Whether user can add/edit/delete at current operating level */
  canContribute: boolean;
  /** Whether user can perform admin actions at current operating level */
  canAdminister: boolean;
  /** Change operating level for this session */
  setOperatingLevel: (level: UserRole) => Promise<void>;
  /** Refresh user data */
  refreshUser: () => Promise<void>;
}

// Default context for server-side rendering and initial state
const defaultContext: AuthContextValue = {
  user: null,
  isLoading: true,
  error: null,
  operatingLevel: 'admin',
  availableLevels: ['admin'],
  canContribute: true,
  canAdminister: true,
  setOperatingLevel: async () => {},
  refreshUser: async () => {},
};

const AuthContext = createContext<AuthContextValue>(defaultContext);

// =============================================================================
// Auth Token Helpers
// =============================================================================

// Try to import auth helpers from ccp4i2 client (when integrated)
let getAccessToken: () => Promise<string | null>;

try {
  const authModule = require('../../utils/auth-token');
  getAccessToken = authModule.getAccessToken;
} catch {
  getAccessToken = async () => null;
}

/**
 * Fetch with authentication support
 */
async function authFetch(url: string, options: RequestInit = {}): Promise<Response> {
  const headers: Record<string, string> = {};

  const token = await getAccessToken();
  if (token) {
    headers['Authorization'] = `Bearer ${token}`;
  }

  if (options.headers) {
    Object.assign(headers, options.headers);
  }

  if (options.body && !(options.body instanceof FormData)) {
    headers['Content-Type'] = 'application/json';
  }

  return fetch(url, { ...options, headers });
}

/**
 * SWR fetcher for user endpoints
 */
async function userFetcher<T>(url: string): Promise<T> {
  const res = await authFetch(url);
  if (!res.ok) {
    const error = new Error(`Failed to fetch user data: ${res.status}`);
    (error as any).status = res.status;
    throw error;
  }
  return res.json();
}

// =============================================================================
// Mock User for Development
// =============================================================================

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

// =============================================================================
// Provider Component
// =============================================================================

interface AuthProviderProps {
  children: React.ReactNode;
}

export function AuthProvider({ children }: AuthProviderProps) {
  // In development mode without auth, use mock user
  const shouldFetch = REQUIRE_AUTH;

  // Fetch current user
  const {
    data: user,
    error,
    isLoading: userLoading,
    mutate: mutateUser,
  } = useSWR<CurrentUser>(
    shouldFetch ? `${USERS_API_BASE}/me/` : null,
    userFetcher,
    {
      revalidateOnFocus: false,
      revalidateOnReconnect: true,
      dedupingInterval: 30000, // 30 seconds
    }
  );

  // Local state for operating level (optimistic updates)
  const [localOperatingLevel, setLocalOperatingLevel] = useState<UserRole | null>(null);

  // Derive effective values
  const effectiveUser = shouldFetch ? user : mockUser;
  const isLoading = shouldFetch ? userLoading : false;

  // Operating level: use local state if set, otherwise from user data
  const operatingLevel = localOperatingLevel ?? effectiveUser?.operating_level ?? 'admin';

  // Available levels: In dev mode, provide all levels for testing
  // In production, base on user's assigned role
  const availableLevels: UserRole[] = shouldFetch
    ? (effectiveUser?.role === 'admin' ? ['user', 'contributor', 'admin'] :
       effectiveUser?.role === 'contributor' ? ['user', 'contributor'] : ['user'])
    : ['user', 'contributor', 'admin']; // All levels available in dev mode for testing

  // Permission helpers based on operating level
  const canContribute = operatingLevel === 'contributor' || operatingLevel === 'admin';
  const canAdminister = operatingLevel === 'admin';

  // Sync local operating level from server response
  useEffect(() => {
    if (user?.operating_level && localOperatingLevel === null) {
      setLocalOperatingLevel(user.operating_level);
    }
  }, [user?.operating_level, localOperatingLevel]);

  // Change operating level
  const setOperatingLevel = useCallback(async (level: UserRole) => {
    if (!REQUIRE_AUTH) {
      // In development mode, just update local state (no server call)
      setLocalOperatingLevel(level);
      return;
    }

    // Optimistic update
    setLocalOperatingLevel(level);

    try {
      const res = await authFetch(`${USERS_API_BASE}/me/operating-level/`, {
        method: 'POST',
        body: JSON.stringify({ level }),
      });

      if (!res.ok) {
        throw new Error('Failed to update operating level');
      }

      const data: OperatingLevelInfo = await res.json();
      setLocalOperatingLevel(data.operating_level);

      // Refresh user data to sync state
      mutateUser();
    } catch (err) {
      // Revert on error
      setLocalOperatingLevel(effectiveUser?.operating_level ?? 'admin');
      throw err;
    }
  }, [effectiveUser, mutateUser]);

  // Refresh user data
  const refreshUser = useCallback(async () => {
    await mutateUser();
  }, [mutateUser]);

  const contextValue: AuthContextValue = {
    user: effectiveUser ?? null,
    isLoading,
    error: error ?? null,
    operatingLevel,
    availableLevels,
    canContribute,
    canAdminister,
    setOperatingLevel,
    refreshUser,
  };

  return (
    <AuthContext.Provider value={contextValue}>
      {children}
    </AuthContext.Provider>
  );
}

// =============================================================================
// Hook
// =============================================================================

/**
 * Hook to access auth context.
 *
 * @example
 * ```tsx
 * function MyComponent() {
 *   const { user, canContribute, operatingLevel, setOperatingLevel } = useAuth();
 *
 *   if (!canContribute) {
 *     return <ReadOnlyMessage />;
 *   }
 *
 *   return <EditButton onClick={handleEdit} />;
 * }
 * ```
 */
export function useAuth(): AuthContextValue {
  const context = useContext(AuthContext);
  if (!context) {
    throw new Error('useAuth must be used within an AuthProvider');
  }
  return context;
}

// =============================================================================
// Role Display Helpers
// =============================================================================

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
