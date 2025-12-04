"use client";

import React, { ReactNode } from "react";
import { Provider } from "react-redux";
import { configureStore } from "@reduxjs/toolkit";
import { atomInfoCardsReducer, menusReducer, sliceNDiceReducer } from "moorhen";
import {
  moleculesReducer,
  mapsReducer,
  mouseSettingsReducer,
  backupSettingsReducer,
  shortcutSettingsReducer,
  labelSettingsReducer,
  sceneSettingsReducer,
  generalStatesReducer,
  hoveringStatesReducer,
  modalsReducer,
  mapContourSettingsReducer,
  moleculeMapUpdateReducer,
  sharedSessionReducer,
  refinementSettingsReducer,
  lhasaReducer,
  overlaysReducer,
  glRefSliceReducer,
  MoorhenReduxStore,
} from "moorhen";

export function ClientStoreProvider({ children }: { children: ReactNode }) {
  const store = configureStore({
    reducer: {
      molecules: moleculesReducer,
      maps: mapsReducer,
      mouseSettings: mouseSettingsReducer,
      backupSettings: backupSettingsReducer,
      shortcutSettings: shortcutSettingsReducer,
      labelSettings: labelSettingsReducer,
      sceneSettings: sceneSettingsReducer,
      generalStates: generalStatesReducer,
      hoveringStates: hoveringStatesReducer,
      modals: modalsReducer,
      mapContourSettings: mapContourSettingsReducer,
      moleculeMapUpdate: moleculeMapUpdateReducer,
      sharedSession: sharedSessionReducer,
      refinementSettings: refinementSettingsReducer,
      lhasa: lhasaReducer,
      sliceNDice: sliceNDiceReducer,
      //jsonValidation: jsonValidationReducer,
      //mrParse: mrPar,
      glRef: glRefSliceReducer,
      overlays: overlaysReducer,
      menus: menusReducer,
      atomInfoCards: atomInfoCardsReducer,
    },
    middleware: (getDefaultMiddleware) =>
      getDefaultMiddleware({
        serializableCheck: false,
      }),
  });
  //Note here am discarding the store I just created in favour of staic imported one.
  //This is *not* good !  Imposes a single store for the whole app.
  return <Provider store={MoorhenReduxStore}>{children}</Provider>;
}
