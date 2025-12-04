import { configureStore } from "@reduxjs/toolkit";
import {
  moleculesReducer,
  mapsReducer,
  mouseSettingsReducer,
  backupSettingsReducer,
  shortcutSettingsReducer,
  labelSettingsReducer,
  sceneSettingsReducer,
  //miscAppSettingsReducer,
  generalStatesReducer,
  hoveringStatesReducer,
  modalsReducer,
  mapContourSettingsReducer,
  glRefSliceReducer,
  moleculeMapUpdateReducer,
  sharedSessionReducer,
  refinementSettingsReducer,
  lhasaReducer,
  overlaysReducer,
} from "moorhen";
import { createContext } from "react";

export const store = configureStore({
  reducer: {
    molecules: moleculesReducer,
    maps: mapsReducer,
    mouseSettings: mouseSettingsReducer,
    backupSettings: backupSettingsReducer,
    shortcutSettings: shortcutSettingsReducer,
    labelSettings: labelSettingsReducer,
    sceneSettings: sceneSettingsReducer,
    //miscAppSettings: miscAppSettingsReducer,
    generalStates: generalStatesReducer,
    hoveringStates: hoveringStatesReducer,
    modals: modalsReducer,
    mapContourSettings: mapContourSettingsReducer,
    moleculeMapUpdate: moleculeMapUpdateReducer,
    sharedSession: sharedSessionReducer,
    refinementSettings: refinementSettingsReducer,
    lhasa: lhasaReducer,
    overlays: overlaysReducer,
    glRef: glRefSliceReducer,
  },
  middleware: (getDefaultMiddleware) =>
    getDefaultMiddleware({
      serializableCheck: false, //{ignoredActions: [FLUSH, REHYDRATE, PAUSE, PERSIST, PURGE, REGISTER] },
    }),
});

export type AppDispatch = typeof store.dispatch;
export type RootState = ReturnType<typeof store.getState>;
export const appStoreContext = createContext(null);
