/**
 * The `style="..."` strings the Python report API writes on its elements
 * (`addText(text=..., style="color:red;")`), turned into a React style
 * object. An html-era convention: the React elements should grow explicit
 * properties (a tone, an emphasis) and the 586 call sites migrate to them;
 * until then this is how a report's bold or red reaches the screen.
 */
export function cssToDict(cssText: string): Record<string, string> {
  const regex = /([\w-]*)\s*:\s*([^;]*)/g;
  const properties: Record<string, string> = {};
  let match: RegExpExecArray | null;
  while ((match = regex.exec(cssText))) {
    const camelCase = match[1].trim().replace(/-(.)/g, (_m, p) => p.toUpperCase());
    if (camelCase) properties[camelCase] = match[2].trim();
  }
  return properties;
}
