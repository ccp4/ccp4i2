import { configureStore } from "@reduxjs/toolkit";
// @ts-expect-error - MoorhenStoreReducers exists at runtime but moorhen's .d.ts doesn't declare it
import { MoorhenStoreReducers } from "moorhen";
import { createContext } from "react";

export const store = configureStore({
  reducer: MoorhenStoreReducers,
  middleware: (getDefaultMiddleware) =>
    getDefaultMiddleware({
      serializableCheck: false,
    }),
});

export type AppDispatch = typeof store.dispatch;
export type RootState = ReturnType<typeof store.getState>;
export const appStoreContext = createContext(null);
