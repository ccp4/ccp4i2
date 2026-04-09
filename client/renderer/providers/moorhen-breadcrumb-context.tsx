"use client";

import React, { createContext, useContext, useState, useCallback, useMemo } from "react";

interface BreadcrumbItem {
  label: string;
  href: string;
}

interface MoorhenBreadcrumbContextValue {
  breadcrumbs: BreadcrumbItem[];
  setBreadcrumbs: (items: BreadcrumbItem[]) => void;
}

const MoorhenBreadcrumbContext =
  createContext<MoorhenBreadcrumbContextValue | null>(null);

export const MoorhenBreadcrumbProvider: React.FC<{
  children: React.ReactNode;
}> = ({ children }) => {
  const [breadcrumbs, setBreadcrumbs] = useState<BreadcrumbItem[]>([]);

  const setBreadcrumbsMemo = useCallback((items: BreadcrumbItem[]) => {
    setBreadcrumbs(items);
  }, []);

  const value = useMemo(
    () => ({ breadcrumbs, setBreadcrumbs: setBreadcrumbsMemo }),
    [breadcrumbs, setBreadcrumbsMemo]
  );

  return (
    <MoorhenBreadcrumbContext.Provider value={value}>
      {children}
    </MoorhenBreadcrumbContext.Provider>
  );
};

export const useMoorhenBreadcrumbs = (): MoorhenBreadcrumbContextValue => {
  const context = useContext(MoorhenBreadcrumbContext);
  if (!context) {
    throw new Error(
      "useMoorhenBreadcrumbs must be used within a MoorhenBreadcrumbProvider"
    );
  }
  return context;
};
