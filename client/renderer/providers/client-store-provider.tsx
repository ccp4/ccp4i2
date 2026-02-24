"use client";

import React, { ReactNode, useRef } from "react";
import { Provider } from "react-redux";
import { configureStore } from "@reduxjs/toolkit";
// @ts-expect-error - MoorhenStoreReducers exists at runtime but moorhen's .d.ts doesn't declare it
import { MoorhenStoreReducers } from "moorhen";

export function ClientStoreProvider({ children }: { children: ReactNode }) {
  const storeRef = useRef(
    configureStore({
      reducer: MoorhenStoreReducers,
      middleware: (getDefaultMiddleware) =>
        getDefaultMiddleware({
          serializableCheck: false,
        }),
    })
  );
  return <Provider store={storeRef.current}>{children}</Provider>;
}
