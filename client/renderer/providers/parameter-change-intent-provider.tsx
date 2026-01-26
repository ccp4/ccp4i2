import React, {
  createContext,
  useContext,
  useState,
  useCallback,
  useMemo,
  ReactNode,
} from "react";

export interface ParameterChangeIntent {
  jobId: number;
  parameterPath: string;
  reason: string;
  previousValue: any;
  timestamp: number;
}

interface ParameterChangeIntentContextType {
  // Legacy single-intent API (for backward compatibility)
  intent: ParameterChangeIntent | null;
  setIntent: (intent: Omit<ParameterChangeIntent, "timestamp">) => void;
  clearIntent: () => void;

  // New multi-intent API
  setIntentForPath: (
    intent: Omit<ParameterChangeIntent, "timestamp">
  ) => void;
  clearIntentForPath: (parameterPath: string) => void;
  getIntentForPath: (parameterPath: string) => ParameterChangeIntent | null;
  wasRecentlyChanged: (parameterPath: string, withinMs?: number) => boolean;
  getAllIntents: () => Map<string, ParameterChangeIntent>;
}

const ParameterChangeIntentContext = createContext<
  ParameterChangeIntentContextType | undefined
>(undefined);

// Default timeout for considering an intent "recent" (2 seconds)
const DEFAULT_RECENT_THRESHOLD_MS = 2000;

// Auto-cleanup threshold - intents older than this are removed (30 seconds)
// Needs to be long enough for file upload + server processing + digest calculation
const CLEANUP_THRESHOLD_MS = 30000;

export const ParameterChangeIntentProvider: React.FC<{
  children: ReactNode;
}> = ({ children }) => {
  // Map of parameterPath -> intent (supports multiple simultaneous edits)
  const [intents, setIntents] = useState<Map<string, ParameterChangeIntent>>(
    new Map()
  );

  // Legacy single intent (most recent one, for backward compatibility)
  const [lastIntent, setLastIntent] = useState<ParameterChangeIntent | null>(
    null
  );

  // Cleanup old intents periodically
  const cleanupOldIntents = useCallback(() => {
    const now = Date.now();
    setIntents((prev) => {
      const newMap = new Map(prev);
      let changed = false;
      for (const [path, intent] of newMap) {
        if (now - intent.timestamp > CLEANUP_THRESHOLD_MS) {
          newMap.delete(path);
          changed = true;
        }
      }
      return changed ? newMap : prev;
    });
  }, []);

  // Set intent for a specific path
  const setIntentForPath = useCallback(
    (intentWithoutTimestamp: Omit<ParameterChangeIntent, "timestamp">) => {
      const intent: ParameterChangeIntent = {
        ...intentWithoutTimestamp,
        timestamp: Date.now(),
      };

      setIntents((prev) => {
        const newMap = new Map(prev);
        newMap.set(intent.parameterPath, intent);
        return newMap;
      });

      // Also update legacy single intent
      setLastIntent(intent);

      // Schedule cleanup
      setTimeout(cleanupOldIntents, CLEANUP_THRESHOLD_MS + 100);
    },
    [cleanupOldIntents]
  );

  // Clear intent for a specific path
  const clearIntentForPath = useCallback((parameterPath: string) => {
    setIntents((prev) => {
      if (!prev.has(parameterPath)) return prev;
      const newMap = new Map(prev);
      newMap.delete(parameterPath);
      return newMap;
    });
  }, []);

  // Get intent for a specific path
  const getIntentForPath = useCallback(
    (parameterPath: string): ParameterChangeIntent | null => {
      return intents.get(parameterPath) || null;
    },
    [intents]
  );

  // Check if a path was recently changed by user
  const wasRecentlyChanged = useCallback(
    (
      parameterPath: string,
      withinMs: number = DEFAULT_RECENT_THRESHOLD_MS
    ): boolean => {
      const intent = intents.get(parameterPath);
      if (!intent) return false;
      return Date.now() - intent.timestamp < withinMs;
    },
    [intents]
  );

  // Get all current intents
  const getAllIntents = useCallback((): Map<string, ParameterChangeIntent> => {
    return new Map(intents);
  }, [intents]);

  // Legacy API: setIntent (updates both single and multi-intent)
  const setIntent = useCallback(
    (intentWithoutTimestamp: Omit<ParameterChangeIntent, "timestamp">) => {
      setIntentForPath(intentWithoutTimestamp);
    },
    [setIntentForPath]
  );

  // Legacy API: clearIntent (clears the last intent path)
  const clearIntent = useCallback(() => {
    if (lastIntent) {
      clearIntentForPath(lastIntent.parameterPath);
    }
    setLastIntent(null);
  }, [lastIntent, clearIntentForPath]);

  const contextValue = useMemo(
    () => ({
      // Legacy API
      intent: lastIntent,
      setIntent,
      clearIntent,
      // New multi-intent API
      setIntentForPath,
      clearIntentForPath,
      getIntentForPath,
      wasRecentlyChanged,
      getAllIntents,
    }),
    [
      lastIntent,
      setIntent,
      clearIntent,
      setIntentForPath,
      clearIntentForPath,
      getIntentForPath,
      wasRecentlyChanged,
      getAllIntents,
    ]
  );

  return (
    <ParameterChangeIntentContext.Provider value={contextValue}>
      {children}
    </ParameterChangeIntentContext.Provider>
  );
};

export const useParameterChangeIntent = () => {
  const context = useContext(ParameterChangeIntentContext);
  if (!context) {
    // Return safe defaults when used outside provider
    return {
      intent: null,
      setIntent: (intent: Omit<ParameterChangeIntent, "timestamp">) => {
        console.warn(
          `Attempt to set intent outside context: ${JSON.stringify(intent)}`
        );
      },
      clearIntent: () => {},
      setIntentForPath: (intent: Omit<ParameterChangeIntent, "timestamp">) => {
        console.warn(
          `Attempt to set intent outside context: ${JSON.stringify(intent)}`
        );
      },
      clearIntentForPath: (_parameterPath: string) => {},
      getIntentForPath: (_parameterPath: string) => null,
      wasRecentlyChanged: (_parameterPath: string, _withinMs?: number) => false,
      getAllIntents: () => new Map<string, ParameterChangeIntent>(),
    };
  }
  return context;
};
