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
