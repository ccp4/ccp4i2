import Script from "next/script";
import { PropsWithChildren, useEffect, useRef } from "react";
import { useCCP4i2Window } from "../app-context";

export const CootProvider: React.FC<PropsWithChildren> = (props) => {
  const { setCootModule } = useCCP4i2Window();
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
      // if it's moorhen.wasm, use a custom dir
      if (path.endsWith("moorhen.wasm")) return "/moorhen.wasm";
      if (path.endsWith("mtz.wasm")) return prefix + path;
      // otherwise, use the default, the prefix (JS file's dir) + the path
      return prefix + path;
    },
  };

  return (
    <>
      <Script
        src="/moorhen.js"
        strategy="lazyOnload"
        id="moorhen-script-element"
        onLoad={async (arg) => {
          const cootModule =
            //@ts-ignore
            createCootModule(createArgs).then((module: any) => {
              //@ts-ignore
              setCootModule(module);
              scriptElement.current = Array.from(
                document.getElementsByTagName("script")
              ).find((htmlElement: HTMLElement) => {
                return htmlElement.getAttribute("src") === "/moorhen.js";
              });
            });
        }}
      />
      {props.children}
    </>
  );
};
