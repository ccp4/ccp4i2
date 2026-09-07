/**
 * Which children a container draws.
 *
 * Two qualifiers say "not this one": `hidden` (from a PHIL `.style`, or a
 * def.xml) and an `expertLevel` above the level the user chose. A container
 * is worth drawing only if some leaf beneath it is — most container scopes
 * carry no expertLevel at all, so testing the container alone leaves a stack
 * of empty accordions where the expert parameters used to be. This mirrors
 * the recursion the Qt UI did in `nestedAutoGenerate`.
 *
 * Kept apart from ccontainer.tsx so it can be tested without importing the
 * element graph.
 */

/**
 * An item's expert level, wherever it happens to live: PHIL-generated
 * containers set it directly as a qualifier; def.xml puts it inside
 * guiDefinition, which is where the great majority of tasks keep it.
 */
export const expertLevelOf = (item: any): number | undefined => {
  const qualifiers = item?._qualifiers;
  const level =
    qualifiers?.expertLevel ?? qualifiers?.guiDefinition?.expertLevel;
  return typeof level === "number" ? level : undefined;
};

export const isHidden = (item: any): boolean =>
  item?._qualifiers?.hidden === true;

export const passesExpertLevel = (
  item: any,
  maxExpertLevel: number
): boolean => {
  const level = expertLevelOf(item);
  if (level !== undefined && level > maxExpertLevel) return false;
  if (item?._baseClass !== "CContainer") return true;
  const children = item?._value;
  if (!children || typeof children !== "object") return false;
  return Object.values(children).some(
    (child) => !isHidden(child) && passesExpertLevel(child, maxExpertLevel)
  );
};

/** Whether a container should draw `item` at all. */
export const shouldRenderChild = (
  item: any,
  maxExpertLevel: number | undefined
): boolean => {
  if (isHidden(item)) return false;
  if (maxExpertLevel === undefined) return true;
  return passesExpertLevel(item, maxExpertLevel);
};
