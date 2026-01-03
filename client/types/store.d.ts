// types/store.d.ts or inline in your file
export interface StoreSchema {
  CCP4Dir: string;              // Path to CCP4 installation (for $CLIBD, $CBIN, etc.)
  projectRoot: string;          // Path to ccp4i2 project (where .venv lives)
  zoomLevel: number;
  devMode: boolean;
  CCP4I2_PROJECTS_DIR: string;
  theme: "light" | "dark";
}
