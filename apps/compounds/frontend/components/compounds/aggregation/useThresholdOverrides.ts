'use client';

import { useCallback, useEffect, useMemo, useState } from 'react';
import { useAuth } from '@/lib/compounds/auth-context';
import {
  AggregationResponse,
  isCompactResponse,
  ProtocolInfo,
} from '@/types/compounds/aggregation';
import type { ProtocolThresholdUpdate } from '../ProtocolThresholdQuickEdit';

/**
 * Frontend-only patches over `data.protocols` so the threshold quick-edit
 * dialog can update colours live without re-fetching. The override map
 * resets when `data` changes (the next fetch reflects server state).
 *
 * Bundles the editing-protocol dialog state because it travels with the
 * same lifecycle. `handleEditProtocol` is a stable callback so it can sit
 * above any guard-clause early returns in the consumer (React #310).
 */
export function useThresholdOverrides(data: AggregationResponse | null | undefined) {
  const { canContribute } = useAuth();
  const [overrides, setOverrides] = useState<Record<string, ProtocolThresholdUpdate>>({});
  const [editingProtocol, setEditingProtocol] = useState<ProtocolInfo | null>(null);

  useEffect(() => {
    setOverrides({});
  }, [data]);

  const handleProtocolSaved = useCallback(
    (id: string, update: ProtocolThresholdUpdate) => {
      setOverrides((prev) => ({ ...prev, [id]: update }));
    },
    [],
  );

  const handleEditProtocol = useCallback(
    (protocol: ProtocolInfo) => setEditingProtocol(protocol),
    [],
  );

  const effectiveData = useMemo(() => {
    if (!data || !isCompactResponse(data)) return data;
    if (Object.keys(overrides).length === 0) return data;
    return {
      ...data,
      protocols: data.protocols.map((p) =>
        overrides[p.id] ? { ...p, ...overrides[p.id] } : p,
      ),
    };
  }, [data, overrides]);

  return {
    effectiveData,
    editingProtocol,
    setEditingProtocol,
    handleProtocolSaved,
    onEditProtocol: canContribute ? handleEditProtocol : undefined,
  };
}
