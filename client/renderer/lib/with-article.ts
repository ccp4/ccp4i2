/** "a map", "an atomic model": a noun with its indefinite article. */
export function withArticle(noun: string): string {
  return `${/^[aeiou]/i.test(noun) ? "an" : "a"} ${noun}`;
}
