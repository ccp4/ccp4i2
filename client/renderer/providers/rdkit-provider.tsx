import Script from "next/script";
import { PropsWithChildren, useEffect, useRef } from "react";
import { useCCP4i2Window } from "../app-context";

export const RdkitProvider: React.FC<PropsWithChildren> = (props) => {
  const { setRdkitModule } = useCCP4i2Window();
  const scriptElement = useRef<HTMLElement | null | undefined>(null);

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
    locateFile(path: string, prefix: string) {
      // if it's rdkit.wasm, use a custom dir
      if (path.endsWith("RDKit_minimal.wasm")) return "/RDKit_minimal.wasm";
      if (path.endsWith("RDKit_minimal.wasm")) return prefix + path;
      // otherwise, use the default, the prefix (JS file's dir) + the path
      return prefix + path;
    },
  };

  return (
    <>
      <Script
        src="/RDKit_minimal.js"
        strategy="lazyOnload"
        id="rdkit-script-element"
        onLoad={async (arg) => {
          const cootModule =
            //@ts-ignore
            initRDKitModule(createArgs).then((module: any) => {
              //@ts-ignore
              setRdkitModule(module);
              scriptElement.current = Array.from(
                document.getElementsByTagName("script")
              ).find((htmlElement: HTMLElement) => {
                return htmlElement.getAttribute("src") === "/RDKit_minimal.js";
              });
            });
        }}
      />
      {props.children}
    </>
  );
};
