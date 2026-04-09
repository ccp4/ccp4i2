import { useMemo, useRef, useEffect } from "react";
import { CCP4i2ReportElementProps } from "./CCP4i2ReportElement";
import { parse } from "acorn";

// Add a type declaration for window.main to satisfy TypeScript
declare global {
  interface Window {
    main?: () => void;
  }
}

export const CCP4i2ReportSVG: React.FC<CCP4i2ReportElementProps> = (props) => {
  const svgContainerRef = useRef<HTMLDivElement>(null);

  // Clean up SVG namespaces for rendering
  const cleanedInnerHTML = useMemo(() => {
    let html = props.item.innerHTML;
    // Remove all namespace prefixes like ns0:
    html = html.replace(/<\s*\/?\s*ns\d*:/g, (match) =>
      match.replace(/ns\d*:/, "")
    );
    // Remove xmlns:ns0="..." attributes
    html = html.replace(/xmlns:ns\d+="[^"]*"/g, "");
    return html;
  }, [props.item.innerHTML]);

  const interactiveSVG: string | null = useMemo(() => {
    try {
      // Parse the SVG string
      const parser = new DOMParser();
      const doc = parser.parseFromString(cleanedInnerHTML, "image/svg+xml");

      // Extract the script content
      const script = doc.querySelector("script");
      if (!script) return null; // No script found, return null

      let scriptContent = script.textContent || "";
      script.parentNode?.removeChild(script); // Remove script from SVG

      // Serialize the cleaned SVG back to a string
      const cleanedSVG = new XMLSerializer().serializeToString(doc);

      // Execute the script (if any)
      if (scriptContent) {
        // WARNING: Only do this if you trust the SVG source!
        // eslint-disable-next-line no-new-func
        // Parse the script
        const ast = parse(scriptContent, { ecmaVersion: 2020 });
        // Collect top-level symbols
        const symbols: string[] = [];
        for (const node of ast.body) {
          if (node.type === "VariableDeclaration") {
            for (const decl of node.declarations) {
              if (decl.id.type === "Identifier") {
                symbols.push(decl.id.name);
              }
            }
          } else if (node.type === "FunctionDeclaration" && node.id) {
            symbols.push(node.id.name);
          }
        }
        for (const symbol of symbols) {
          scriptContent = scriptContent + `\nwindow.${symbol} = ${symbol};\n`;
        }

        // Do NOT call window.main() here!
        // scriptContent = scriptContent + `\nwindow.main();\n`;

        // Evaluate the script to define functions/variables on window
        new Function(scriptContent)();
      }
      return cleanedSVG;
    } catch (error) {
      console.error("Error processing SVG:", error);
      return null;
    }
  }, [cleanedInnerHTML]);

  // Call window.main() only after SVG is in the DOM
  useEffect(() => {
    if (interactiveSVG && typeof window.main === "function") {
      window.main();
    }
  }, [interactiveSVG]);

  return (
    <div
      ref={svgContainerRef}
      style={{ marginTop: "1em" }}
      dangerouslySetInnerHTML={{ __html: interactiveSVG || cleanedInnerHTML }}
    />
  );
};
