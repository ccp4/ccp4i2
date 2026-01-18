import Script from "next/script";
import { PropsWithChildren, useEffect, useRef, useState } from "react";
import { useCCP4i2Window } from "../app-context";

// Check if running in Electron (no COEP restrictions) vs web browser (needs API route for CORP headers)
const isElectron = typeof window !== "undefined" && !!(window as any).electronAPI;

// Cache for worker blob URLs to avoid refetching
const workerBlobCache: Map<string, string> = new Map();

/**
 * Pre-fetch and cache worker scripts for COEP compatibility.
 * Workers in COEP contexts can only be loaded from blob: URLs.
 * We pre-fetch the scripts and create blob URLs that can be used synchronously.
 *
 * CootWorker.js uses importScripts("./moorhen.js") internally, which fails in
 * COEP contexts even with absolute URLs (browser cancels the requests).
 * Solution: inline the moorhen.js content directly into CootWorker.js so no
 * network requests are needed from within the worker.
 */
const prefetchWorkerScripts = async () => {
  if (isElectron || typeof window === "undefined") return;

  const cootWorkerUrl = "/api/moorhen/CootWorker.js";
  if (workerBlobCache.has(cootWorkerUrl)) return;

  try {
    // Fetch both CootWorker.js and moorhen.js (we use 32-bit since 64-bit has issues)
    const [cootWorkerRes, moorhenRes] = await Promise.all([
      fetch("/api/moorhen/CootWorker.js"),
      fetch("/api/moorhen/moorhen.js"),
    ]);

    if (!cootWorkerRes.ok || !moorhenRes.ok) {
      console.error("[Moorhen] Failed to fetch worker scripts");
      return;
    }

    let cootWorkerText = await cootWorkerRes.text();
    const moorhenText = await moorhenRes.text();

    // Replace the entire if/try/catch/else block that loads moorhen scripts.
    // Since we're prepending moorhen.js, createCootModule is already defined.
    // The exact minified pattern in CootWorker.js is:
    //   if(t&&!o)try{importScripts("./moorhen64.js"),r=createCoot64Module,a="moorhen64.js"}catch(e){...}else importScripts("./moorhen.js"),r=createCootModule,a="moorhen.js"
    // Replace with just: r=createCootModule,a="moorhen.js"
    cootWorkerText = cootWorkerText.replace(
      /if\(t&&!o\)try\{importScripts\("\.\/moorhen64\.js"\),r=createCoot64Module,a="moorhen64\.js"\}catch\(e\)\{[^}]+\}else importScripts\("\.\/moorhen\.js"\),r=createCootModule,a="moorhen\.js"/g,
      'r=createCootModule,a="moorhen.js"'
    );

    // Create a blob URL for moorhen.js that can be used by Emscripten to spawn pthread workers
    const moorhenBlob = new Blob([moorhenText], { type: "application/javascript" });
    const moorhenBlobUrl = URL.createObjectURL(moorhenBlob);

    // Add a wrapper that intercepts createCootModule to:
    // 1. Inject the correct mainScriptUrlOrBlob for pthread workers (blob URL)
    // 2. Provide a locateFile that returns absolute URLs (required in blob worker context)
    const workerWrapper = `
// Store the blob URL for moorhen.js so pthread workers can be created
var __moorhenBlobUrl = "${moorhenBlobUrl}";
// Base URL for resolving resources (from main page origin)
var __baseUrl = "${typeof window !== "undefined" ? window.location.origin : ""}/api/moorhen/";

if (typeof createCootModule === "function") {
  var originalCreateCootModule = createCootModule;
  createCootModule = function(args) {
    // Force mainScriptUrlOrBlob to use our blob URL for spawning pthread workers
    // And provide a locateFile that returns absolute URLs
    var wrappedArgs = Object.assign({}, args, {
      mainScriptUrlOrBlob: __moorhenBlobUrl,
      locateFile: function(path, prefix) {
        // In blob worker context, we need absolute URLs
        // Route through the API for CORP headers
        if (path.endsWith(".wasm") || path.endsWith(".data")) {
          return __baseUrl + path;
        }
        // For other files, try the original locateFile if available
        if (args && args.locateFile) {
          return args.locateFile(path, prefix);
        }
        // Fallback to absolute URL
        return __baseUrl + path;
      }
    });
    return originalCreateCootModule(wrappedArgs);
  };
}
`;

    // Combine: moorhen.js (defines createCootModule), wrapper (patches it), then patched CootWorker
    const combinedScript = moorhenText + "\n\n" + workerWrapper + "\n\n// === CootWorker.js (patched) ===\n\n" + cootWorkerText;

    const blob = new Blob([combinedScript], { type: "application/javascript" });
    const blobUrl = URL.createObjectURL(blob);
    workerBlobCache.set(cootWorkerUrl, blobUrl);
    workerBlobCache.set("/api/moorhen/CootWorker.js", blobUrl);
  } catch (err) {
    console.error("[Moorhen] Failed to prepare worker scripts:", err);
  }
};

/**
 * Patch the global Worker constructor for COEP compatibility in web mode.
 *
 * In a COEP-enabled page, Workers can only be loaded from:
 * 1. Same-origin blob: URLs
 * 2. data: URLs
 *
 * Regular URL-based workers (like `new Worker('/CootWorker.js')`) are blocked.
 * This patch intercepts Worker construction and converts URL-based workers to blob-based ones.
 */
const patchWorkerForCOEP = () => {
  if (isElectron || typeof window === "undefined") return;

  // Only patch once
  if ((window as any).__workerPatchedForCOEP) return;
  (window as any).__workerPatchedForCOEP = true;

  const OriginalWorker = window.Worker;

  // Create a patched Worker constructor
  function PatchedWorker(
    this: Worker,
    scriptURL: string | URL,
    options?: WorkerOptions
  ): Worker {
    const urlString = scriptURL.toString();

    // If it's already a blob URL or data URL, use it directly
    if (urlString.startsWith("blob:") || urlString.startsWith("data:")) {
      return new OriginalWorker(scriptURL, options);
    }

    // Resolve the URL to check cache
    const absoluteUrl = new URL(urlString, window.location.href).href;
    const pathname = new URL(absoluteUrl).pathname;

    // Check cache for various forms of the URL
    let cachedBlobUrl = workerBlobCache.get(urlString) ||
      workerBlobCache.get(absoluteUrl) ||
      workerBlobCache.get(pathname);

    if (cachedBlobUrl) {
      return new OriginalWorker(cachedBlobUrl, options);
    }

    // If not cached, the worker likely won't work due to COEP restrictions
    console.warn(`[Moorhen] Worker script not pre-cached: ${urlString}. May fail due to COEP.`);
    return new OriginalWorker(scriptURL, options);
  }

  // Copy static properties and prototype
  PatchedWorker.prototype = OriginalWorker.prototype;
  Object.setPrototypeOf(PatchedWorker, OriginalWorker);

  // Replace the global Worker constructor
  (window as any).Worker = PatchedWorker;

  // Start prefetching worker scripts
  prefetchWorkerScripts();
};

// Apply the patch immediately when this module loads (before moorhen)
if (!isElectron && typeof window !== "undefined") {
  patchWorkerForCOEP();
}

// Detect if we should use 64-bit WASM
// Currently forced to 32-bit - the 64-bit build in moorhen 0.22.7 has a BigInt conversion bug:
// "Cannot convert a BigInt value to a number" during initialization
// This is an issue in the moorhen package itself that needs to be fixed upstream.
// Safari is blocked at moorhen-wrapper level anyway due to WASM threading issues.
const shouldUse64BitWasm = (): boolean => {
  // TODO: Re-enable 64-bit when moorhen fixes the BigInt conversion issue
  // if (typeof navigator === "undefined") return false;
  // const ua = navigator.userAgent;
  // const isSafari = /^((?!chrome|android).)*safari/i.test(ua);
  // const isIOS = /iPad|iPhone|iPod/.test(ua);
  // return !isSafari && !isIOS;
  return false; // Force 32-bit until moorhen 64-bit build is fixed
};

export const CootProvider: React.FC<PropsWithChildren> = (props) => {
  const { setCootModule } = useCCP4i2Window();
  const scriptElement = useRef<HTMLElement | null | undefined>(null);
  const [use64Bit] = useState(() => shouldUse64BitWasm());
  const [scriptLoaded, setScriptLoaded] = useState(false);
  // Store the script blob for passing to Emscripten for worker creation
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
    return () => {
      if (scriptElement.current) {
        scriptElement.current.parentElement?.removeChild(scriptElement.current);
      }
    };
  }, []);

  const createArgs = {
    print(t: string) {
      console.log(["output", t]);
    },
    printErr(t: string) {
      console.error(["output", t]);
    },
    // Pass the script blob directly to Emscripten for spawning workers
    // This is needed because we load the script via blob URL, and API routes
    // don't work properly as worker script sources
    mainScriptUrlOrBlob: isElectron ? undefined : scriptBlob,
    locateFile(path: string, prefix: string) {
      // Route through API to ensure CORP headers are set for COEP compatibility (web only)
      if (path.endsWith("moorhen.wasm") || path.endsWith("moorhen64.wasm")) {
        return moorhenWasm;
      }
      // Handle CootWorker.js specifically - it needs absolute URL in web mode
      if (!isElectron && path === "CootWorker.js") {
        return `/api/moorhen/CootWorker.js`;
      }
      // For other files, route through API if in web mode
      if (!isElectron && (path.endsWith(".wasm") || path.endsWith(".js") || path.endsWith(".data"))) {
        return `/api/moorhen/${path}`;
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
        const module = await (createModule as (args: typeof createArgs) => Promise<any>)(createArgs);
        setCootModule?.(module);
      } catch (err) {
        console.error("[Moorhen] Failed to initialize module:", err);
      }
    };

    initModule();
  }, [scriptLoaded, scriptBlob, use64Bit, setCootModule, moorhenWasm]);

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
            // eslint-disable-next-line @typescript-eslint/no-explicit-any
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
