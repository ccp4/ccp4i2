/*
 * Copyright (C) 2025-2026 Newcastle University
 *
 * This file is part of CCP4i2.
 *
 * CCP4i2 is free software: you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License version 3,
 * modified in accordance with the provisions of the license to address
 * the requirements of UK law.
 *
 * See https://www.ccp4.ac.uk/ccp4license.php for details.
 */
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
