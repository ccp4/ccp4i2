/**
 * Centralized route definitions for the compounds apps.
 *
 * This module provides a single source of truth for all navigation routes.
 * When integrating into the CCP4i2 client, set NEXT_PUBLIC_ROUTE_PREFIX
 * to add a prefix to all routes (e.g., '/compounds' or leave empty for standalone).
 *
 * Usage:
 *   import { routes } from '@/lib/compounds/routes';
 *
 *   // In Link components:
 *   <Link href={routes.registry.targets()}>Targets</Link>
 *
 *   // In router.push:
 *   router.push(routes.registry.compound(compoundId));
 *
 *   // With query params:
 *   router.push(routes.assays.aggregate({ compound: 'NCL-00001' }));
 */

// Route prefix for integration - empty string for standalone, '/compounds' for integrated
const PREFIX = process.env.NEXT_PUBLIC_ROUTE_PREFIX || '';

/**
 * Build a route path with the configured prefix
 */
function route(path: string): string {
  return `${PREFIX}${path}`;
}

/**
 * Build a route with query parameters
 */
function routeWithQuery(path: string, params?: Record<string, string | number | undefined>): string {
  const base = route(path);
  if (!params) return base;

  const searchParams = new URLSearchParams();
  Object.entries(params).forEach(([key, value]) => {
    if (value !== undefined) {
      searchParams.set(key, String(value));
    }
  });

  const queryString = searchParams.toString();
  return queryString ? `${base}?${queryString}` : base;
}

/**
 * All application routes
 */
export const routes = {
  /** Home page */
  home: () => route('/'),

  /** Registry app routes */
  registry: {
    /** Targets list (registry home) */
    targets: () => route('/registry/targets'),

    /** Target detail page (dashboard) */
    target: (id: string | number) => route(`/registry/targets/${id}`),

    /** Target compounds list page */
    targetCompounds: (id: string | number) => route(`/registry/targets/${id}/compounds`),

    /** All compounds list with optional target filter */
    compounds: (params?: { target?: string }) =>
      routeWithQuery('/registry/compounds', params),

    /** Suppliers list */
    suppliers: () => route('/registry/suppliers'),

    /** Supplier detail page */
    supplier: (id: string | number) => route(`/registry/suppliers/${id}`),

    /** Compound detail page */
    compound: (id: string | number) => route(`/registry/compounds/${id}`),

    /** Batch detail page */
    batch: (id: string | number) => route(`/registry/batches/${id}`),

    /** Compound search */
    search: () => route('/registry/search'),

    /** Import compounds */
    import: () => route('/registry/import'),

    /** Register new compound */
    new: () => route('/registry/new'),
  },

  /** Assays app routes */
  assays: {
    /** Assays list (assays home) */
    list: () => route('/assays'),

    /** Assay detail page */
    detail: (id: string | number) => route(`/assays/${id}`),

    /** Protocols list */
    protocols: () => route('/assays/protocols'),

    /** Protocol detail page */
    protocol: (id: string | number) => route(`/assays/protocols/${id}`),

    /** Plate layouts list */
    plateLayouts: () => route('/assays/plate-layouts'),

    /** Plate layout detail/edit page */
    plateLayout: (id: string | number) => route(`/assays/plate-layouts/${id}`),

    /** Dilution series list */
    dilutionSeries: () => route('/assays/dilution-series'),

    /** Data series detail page */
    dataSeries: (id: string | number) => route(`/assays/data-series/${id}`),

    /** Import assay data */
    import: () => route('/assays/import'),

    /** Import Table of Values data */
    importTableOfValues: (params?: { protocol?: string }) =>
      routeWithQuery('/assays/import-tov', params),

    /** Import ADME data */
    importAdme: () => route('/assays/import-adme'),

    /** Aggregation view with optional filters */
    aggregate: (params?: {
      compound?: string;
      target?: string | number;
      protocol?: string | number;
    }) => routeWithQuery('/assays/aggregate', params),
  },

  /** Constructs app routes */
  constructs: {
    /** Plasmids list (constructs home) */
    list: () => route('/constructs'),

    /** Projects list */
    projects: () => route('/constructs/projects'),

    /** Project detail page */
    project: (id: string | number) => route(`/constructs/projects/${id}`),

    /** Plasmid detail page */
    plasmid: (id: string | number) => route(`/constructs/plasmids/${id}`),

    /** Proteins list */
    proteins: () => route('/constructs/proteins'),

    /** Protein detail page */
    protein: (id: string | number) => route(`/constructs/proteins/${id}`),
  },

  /** Admin routes (requires admin operating level) */
  admin: {
    /** Admin home */
    home: () => route('/admin'),

    /** User management */
    users: () => route('/admin/users'),
  },

  /** External app routes (for cross-app navigation when integrated) */
  external: {
    /** CCP4i2 home (only relevant when integrated) */
    ccp4i2: () => '/ccp4i2',

    /** CCP4i2 project */
    ccp4i2Project: (id: string | number) => `/ccp4i2/project/${id}`,
  },
};

/**
 * Get the current route prefix (useful for debugging)
 */
export function getRoutePrefix(): string {
  return PREFIX;
}

/**
 * Check if we're running in integrated mode
 */
export function isIntegratedMode(): boolean {
  return PREFIX !== '';
}
