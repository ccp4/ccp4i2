/*
 * Copyright (C) 2025 Newcastle University
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
export const isElectron = (): boolean => {
  if (typeof window === "undefined") {
    return process.env.NEXT_PUBLIC_IS_ELECTRON === "true";
  }

  return (
    (window as any).process?.type === "renderer" ||
    !!(window as any).electronAPI ||
    process.env.NEXT_PUBLIC_IS_ELECTRON === "true"
  );
};

export const getApiBaseUrl = (): string => {
  if (isElectron()) {
    return "/api/proxy"; // Electron's proxy
  }
  return process.env.NEXT_PUBLIC_API_BASE_URL || "/api";
};

export const getPlatform = (): "electron" | "web" => {
  return isElectron() ? "electron" : "web";
};
