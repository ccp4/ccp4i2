'use client';

import { createContext, useContext, useState, PropsWithChildren } from 'react';
import Script from 'next/script';

// RDKit module type (minimal typing for what we use)
interface RDKitModule {
  get_mol: (smiles: string) => {
    get_svg: (width?: number, height?: number) => string;
    delete: () => void;
  };
}

interface RDKitContextType {
  rdkitModule: RDKitModule | null;
  isLoading: boolean;
  error: string | null;
}

const RDKitContext = createContext<RDKitContextType>({
  rdkitModule: null,
  isLoading: true,
  error: null,
});

export const useRDKit = () => useContext(RDKitContext);

export const RDKitProvider: React.FC<PropsWithChildren> = ({ children }) => {
  const [rdkitModule, setRdkitModule] = useState<RDKitModule | null>(null);
  const [isLoading, setIsLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);

  const createArgs = {
    print(t: string) {
      console.log(['RDKit output', t]);
    },
    printErr(t: string) {
      console.error(['RDKit error', t]);
    },
    locateFile(path: string, prefix: string) {
      if (path.endsWith('RDKit_minimal.wasm')) return '/RDKit_minimal.wasm';
      return prefix + path;
    },
  };

  return (
    <RDKitContext.Provider value={{ rdkitModule, isLoading, error }}>
      <Script
        src="/RDKit_minimal.js"
        strategy="lazyOnload"
        id="rdkit-script"
        onLoad={async () => {
          try {
            // @ts-ignore - initRDKitModule is defined by the loaded script
            const module = await initRDKitModule(createArgs);
            setRdkitModule(module);
            setIsLoading(false);
            console.log('RDKit module loaded successfully');
          } catch (err) {
            console.error('Failed to initialize RDKit:', err);
            setError('Failed to load RDKit module');
            setIsLoading(false);
          }
        }}
        onError={() => {
          setError('Failed to load RDKit script');
          setIsLoading(false);
        }}
      />
      {children}
    </RDKitContext.Provider>
  );
};
