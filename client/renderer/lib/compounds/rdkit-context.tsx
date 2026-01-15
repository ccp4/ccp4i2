/**
 * Re-export RDKit provider from the core providers directory.
 *
 * This file ensures compounds components that import from '@/lib/compounds/rdkit-context'
 * use the same React context as the root layout's RDKitProvider.
 *
 * The Docker build copies compounds components but NOT the compounds rdkit-context.tsx,
 * so this file is used instead, maintaining a single shared context.
 */
export { useRDKit, RDKitProvider } from '../../providers/rdkit-provider';
