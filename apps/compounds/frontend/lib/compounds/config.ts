/**
 * Compounds configuration management.
 *
 * Fetches deployment-specific configuration from the backend,
 * allowing the UI to adapt to different compound ID prefixes.
 *
 * Usage:
 *   // In a component that needs dynamic config:
 *   const { config, isLoading } = useCompoundConfig();
 *   <input placeholder={`${config.compound_id_prefix}-00026...`} />
 *
 *   // For synchronous access (uses cached value):
 *   const prefix = getCompoundPrefix();
 *   const formattedId = formatCompoundId(26042);
 */

import useSWR from 'swr';

// API base path
const API_BASE = '/api/proxy/compounds';

/**
 * Configuration returned by the backend.
 */
export interface CompoundConfig {
  compound_id_prefix: string;
  compound_id_digits: number;
  compound_id_example: string;
}

/**
 * Default configuration (matches backend defaults).
 * Used when config hasn't been fetched yet or if fetch fails.
 */
const DEFAULT_CONFIG: CompoundConfig = {
  compound_id_prefix: 'NCL',
  compound_id_digits: 8,
  compound_id_example: 'NCL-00026042',
};

/**
 * Cached configuration for synchronous access.
 * Updated when useCompoundConfig fetches new data.
 */
let cachedConfig: CompoundConfig | null = null;

/**
 * Simple fetcher that doesn't require authentication.
 * The config endpoint is public.
 */
async function fetchConfig(url: string): Promise<CompoundConfig> {
  const res = await fetch(url);
  if (!res.ok) {
    throw new Error(`Failed to fetch config: ${res.status}`);
  }
  const config = await res.json();
  cachedConfig = config;
  return config;
}

/**
 * Hook to get compound configuration.
 *
 * Fetches configuration from the backend and caches it.
 * The configuration is cached aggressively since it rarely changes.
 *
 * @returns Object with config, isLoading, and error
 */
export function useCompoundConfig() {
  const { data, error, isLoading } = useSWR<CompoundConfig>(
    `${API_BASE}/config/`,
    fetchConfig,
    {
      revalidateOnFocus: false,
      revalidateOnReconnect: false,
      dedupingInterval: 300000, // Cache for 5 minutes
      errorRetryCount: 2,
      fallbackData: cachedConfig || DEFAULT_CONFIG,
    }
  );

  return {
    config: data || cachedConfig || DEFAULT_CONFIG,
    isLoading,
    error,
  };
}

/**
 * Get the compound ID prefix (synchronous, uses cached value).
 *
 * Returns the cached value from a previous useCompoundConfig call,
 * or the default "NCL" if not yet fetched.
 */
export function getCompoundPrefix(): string {
  return cachedConfig?.compound_id_prefix || DEFAULT_CONFIG.compound_id_prefix;
}

/**
 * Get the number of digits in compound IDs (synchronous).
 */
export function getCompoundDigits(): number {
  return cachedConfig?.compound_id_digits || DEFAULT_CONFIG.compound_id_digits;
}

/**
 * Format a registration number as a compound ID.
 *
 * @param regNumber - The integer registration number
 * @returns Formatted compound ID (e.g., "NCL-00026042")
 */
export function formatCompoundId(regNumber: number): string {
  const prefix = getCompoundPrefix();
  const digits = getCompoundDigits();
  return `${prefix}-${regNumber.toString().padStart(digits, '0')}`;
}

/**
 * Create a regex pattern for matching compound IDs.
 *
 * The pattern matches:
 * - PREFIX-00026042 (standard)
 * - PREFIX26042 (no dash)
 * - PREFIX 26042 (space)
 *
 * @returns RegExp with capturing group for the numeric portion
 */
export function getCompoundPattern(): RegExp {
  const prefix = getCompoundPrefix();
  // Escape any special regex characters in prefix
  const escaped = prefix.replace(/[.*+?^${}()|[\]\\]/g, '\\$&');
  return new RegExp(`${escaped}-?(\\d+)`, 'i');
}

/**
 * Extract registration number from a compound identifier string.
 *
 * @param identifier - String that may contain a compound ID
 * @returns The numeric registration number, or null if not found
 */
export function extractRegNumber(identifier: string): number | null {
  if (!identifier) return null;

  const pattern = getCompoundPattern();
  const match = identifier.match(pattern);
  if (match) {
    return parseInt(match[1], 10);
  }

  // Fall back to plain numeric
  if (/^\d+$/.test(identifier)) {
    return parseInt(identifier, 10);
  }

  return null;
}

/**
 * Get an example compound ID for use in placeholders.
 *
 * @returns Example like "NCL-00026042"
 */
export function getCompoundExample(): string {
  return cachedConfig?.compound_id_example || DEFAULT_CONFIG.compound_id_example;
}
