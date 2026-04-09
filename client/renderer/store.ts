import { configureStore } from "@reduxjs/toolkit";
// @ts-ignore - MoorhenStoreReducers may lack .d.ts depending on build
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
