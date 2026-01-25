import Script from "next/script";
import { PropsWithChildren, useEffect, useRef, useState } from "react";
import { useCCP4i2Window } from "../app-context";

// Check if running in Electron (no COEP restrictions) vs web browser
const isElectron = typeof window !== "undefined" && !!(window as any).electronAPI;

// Detect if we should use 64-bit WASM
// Currently forced to 32-bit - the 64-bit build in moorhen 0.22.7 has a BigInt conversion bug:
// "Cannot convert a BigInt value to a number" during initialization
// This is an issue in the moorhen package itself that needs to be fixed upstream.
// Safari is blocked at moorhen-wrapper level anyway due to WASM threading issues.
const shouldUse64BitWasm = (): boolean => {
  // TODO: Re-enable 64-bit when moorhen fixes the BigInt conversion issue
  return false; // Force 32-bit until moorhen 64-bit build is fixed
};

export const CootProvider: React.FC<PropsWithChildren> = (props) => {
  const { setCootModule } = useCCP4i2Window();
  const scriptElement = useRef<HTMLElement | null | undefined>(null);
  const [use64Bit] = useState(() => shouldUse64BitWasm());
  const [scriptLoaded, setScriptLoaded] = useState(false);
  const [scriptBlob, setScriptBlob] = useState<Blob | null>(null);

  // In web browsers, use API route for moorhen files to ensure CORP headers are set
  // (required for COEP/SharedArrayBuffer support). In Electron, serve directly from public/.
  const apiPrefix = isElectron ? "" : "/api/moorhen";
  const moorhenScript = use64Bit ? "/moorhen64.js" : "/moorhen.js";
  const moorhenWasm = use64Bit ? `${apiPrefix}/moorhen64.wasm` : `${apiPrefix}/moorhen.wasm`;

  // For web mode, we need to fetch and execute the script manually to work with COEP
  // The <Script> tag approach doesn't work because browsers block scripts that don't
  // have proper CORP headers when loaded into a COEP-enabled page
  useEffect(() => {
    if (isElectron || scriptLoaded) return;

    const loadScriptViaFetch = async () => {
      try {
        const scriptUrl = use64Bit ? "/api/moorhen/moorhen64.js" : "/api/moorhen/moorhen.js";
        console.log(`[Moorhen] Fetching script from ${scriptUrl}`);
        const response = await fetch(scriptUrl);
        if (!response.ok) {
          throw new Error(`Failed to fetch moorhen script: ${response.status}`);
        }
        const scriptText = await response.text();

        // Create blob for both script execution AND for Emscripten worker spawning
        const blob = new Blob([scriptText], { type: "application/javascript" });
        setScriptBlob(blob); // Store for later use with Emscripten

        const blobUrl = URL.createObjectURL(blob);

        const script = document.createElement("script");
        script.src = blobUrl;
        script.onload = () => {
          console.log("[Moorhen] Script loaded successfully");
          // Don't revoke the URL yet - we need it for workers
          setScriptLoaded(true);
        };
        script.onerror = (err) => {
          console.error("[Moorhen] Failed to execute script:", err);
          URL.revokeObjectURL(blobUrl);
        };
        document.head.appendChild(script);
        scriptElement.current = script;
      } catch (err) {
        console.error("[Moorhen] Failed to load script:", err);
      }
    };

    loadScriptViaFetch();
  }, [use64Bit, scriptLoaded]);

  useEffect(() => {
    console.log(`[Moorhen] Using ${use64Bit ? "64-bit" : "32-bit"} WASM build, isElectron=${isElectron}`);
    return () => {
      if (scriptElement.current) {
        scriptElement.current.parentElement?.removeChild(scriptElement.current);
      }
    };
  }, [use64Bit]);

  const createArgs = {
    print(t: string) {
      console.log(["output", t]);
    },
    printErr(t: string) {
      console.error(["output", t]);
    },
    // Pass the script blob directly to Emscripten for spawning workers
    // This is needed because we load the script via blob URL
    mainScriptUrlOrBlob: isElectron ? undefined : scriptBlob,
    locateFile(path: string, prefix: string) {
      // Route through API to ensure CORP headers are set for COEP compatibility (web only)
      if (path.endsWith("moorhen.wasm") || path.endsWith("moorhen64.wasm")) {
        return moorhenWasm;
      }
      // For other files, route through API if in web mode
      if (!isElectron && (path.endsWith(".wasm") || path.endsWith(".js") || path.endsWith(".data"))) {
        return `${apiPrefix}/${path}`;
      }
      // In Electron or for other files, use the default prefix + path
      return prefix + path;
    },
  };

  // Initialize the module once the script is loaded (web mode only)
  useEffect(() => {
    if (isElectron || !scriptLoaded || !scriptBlob) return;

    const initModule = async () => {
      const createModule = use64Bit
        ? (window as any).createCoot64Module
        : (window as any).createCootModule;

      if (typeof createModule !== "function") {
        console.error(`[Moorhen] ${use64Bit ? "createCoot64Module" : "createCootModule"} not found`);
        return;
      }

      try {
        console.log("[Moorhen] Initializing module...");
        const module = await (createModule as (args: typeof createArgs) => Promise<any>)(createArgs);
        console.log("[Moorhen] Module initialized successfully");
        setCootModule?.(module);
      } catch (err) {
        console.error("[Moorhen] Failed to initialize module:", err);
      }
    };

    initModule();
  }, [scriptLoaded, scriptBlob, use64Bit, setCootModule, moorhenWasm, apiPrefix]);

  // In Electron mode, use the standard Script component approach
  if (isElectron) {
    return (
      <>
        <Script
          src={moorhenScript}
          strategy="lazyOnload"
          id="moorhen-script-element"
          onLoad={async () => {
            // moorhen.js defines createCootModule, moorhen64.js defines createCoot64Module
            const createModule = use64Bit
              ? (window as any).createCoot64Module
              : (window as any).createCootModule;

            if (typeof createModule !== "function") {
              console.error(`[Moorhen] ${use64Bit ? "createCoot64Module" : "createCootModule"} not found`);
              return;
            }

            const module = await (createModule as (args: typeof createArgs) => Promise<any>)(createArgs);
            setCootModule?.(module);
            scriptElement.current = Array.from(
              document.getElementsByTagName("script")
            ).find((htmlElement: HTMLElement) => {
              return htmlElement.getAttribute("src") === moorhenScript;
            });
          }}
        />
        {props.children}
      </>
    );
  }

  // In web mode, script is loaded via fetch/blob in useEffect
  return <>{props.children}</>;
};
