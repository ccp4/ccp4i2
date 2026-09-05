/**
 * Which registry entry renders an item: its exact class if the registry
 * knows it, else its base kind. A server-side subclass the client has never
 * heard of -- the container class generated for each repeatable PHIL scope,
 * a CString with extra validation -- still renders as what it is.
 *
 * Kept apart from task-element.tsx so it can be tested without importing
 * the element graph.
 */
export const resolveComponentEntry = <T>(
  registry: Record<string, T>,
  item: { _class?: string; _baseClass?: string } | undefined
): T | undefined => {
  if (!item) return undefined;
  const byClass = item._class ? registry[item._class] : undefined;
  if (byClass) return byClass;
  return item._baseClass ? registry[item._baseClass] : undefined;
};
