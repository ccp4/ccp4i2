// Markdown files imported for their raw text (webpack `asset/source`, see
// next.config.ts). Used to embed the scene grammar into the LLM authoring
// prompt scaffold.
declare module "*.md" {
  const content: string;
  export default content;
}
