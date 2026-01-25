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

  // In web browsers, use API route for moorhen files to ensure CORP headers are set
  // (required for COEP/SharedArrayBuffer support). In Electron, serve directly from public/.
  const apiPrefix = isElectron ? "" : "/api/moorhen";
  const moorhenScript = use64Bit ? `${apiPrefix}/moorhen64.js` : `${apiPrefix}/moorhen.js`;
  const moorhenWasm = use64Bit ? `${apiPrefix}/moorhen64.wasm` : `${apiPrefix}/moorhen.wasm`;

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
    locateFile(path: string, prefix: string) {
      // Use the appropriate WASM file based on browser detection
      if (path.endsWith("moorhen.wasm") || path.endsWith("moorhen64.wasm")) {
        return moorhenWasm;
      }
      // In web mode, route files through API for CORP headers
      if (!isElectron && (path.endsWith(".wasm") || path.endsWith(".js") || path.endsWith(".data"))) {
        return `${apiPrefix}/${path}`;
      }
      // otherwise, use the default, the prefix (JS file's dir) + the path
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
