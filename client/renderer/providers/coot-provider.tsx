import Script from "next/script";
import { PropsWithChildren, useEffect, useRef, useState } from "react";
import { useCCP4i2Window } from "../app-context";

// Check if running in Electron (no COEP restrictions) vs web browser (needs API route for CORP headers)
const isElectron = typeof window !== "undefined" && !!(window as any).electronAPI;

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

  // In web browsers, use API route for moorhen files to ensure CORP headers are set
  // (required for COEP/SharedArrayBuffer support). In Electron, serve directly from public/.
  const apiPrefix = isElectron ? "" : "/api/moorhen";
  const moorhenScript = use64Bit ? `${apiPrefix}/moorhen64.js` : `${apiPrefix}/moorhen.js`;
  const moorhenWasm = use64Bit ? `${apiPrefix}/moorhen64.wasm` : `${apiPrefix}/moorhen.wasm`;

  useEffect(() => {
    console.log(`[Moorhen] Using ${use64Bit ? "64-bit" : "32-bit"} WASM build`);
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
    locateFile(path: string, prefix: string) {
      // Route through API to ensure CORP headers are set for COEP compatibility (web only)
      if (path.endsWith("moorhen.wasm") || path.endsWith("moorhen64.wasm")) {
        return moorhenWasm;
      }
      // For other files, route through API if in web mode
      if (!isElectron && (path.endsWith(".wasm") || path.endsWith(".js") || path.endsWith(".data"))) {
        return `/api/moorhen/${path}`;
      }
      // In Electron or for other files, use the default prefix + path
      return prefix + path;
    },
  };

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
};
