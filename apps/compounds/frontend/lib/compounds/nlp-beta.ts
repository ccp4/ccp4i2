/**
 * Temporary beta gating for the NLP command bar.
 *
 * Slice 9 shipped a significant pivot (QuerySpec → CompoundSelector) and
 * the prompt-engineering hasn't yet been validated against real-chemist
 * phrasing at scale. Until the golden-set online eval is run and any
 * prompt-tuning falls out, we hide the "Ask" tile from the landing page
 * for users other than the ones in the allowlist below.
 *
 * This is **UI gating, not security** — the endpoint itself is still
 * protected by:
 *   1. The `COMPOUNDS_NLP_ENABLED` env var (per-instance feature flag)
 *   2. `[IsAuthenticated]` DRF permission
 *   3. The per-user daily cap
 *
 * So a beta user can hit `/nlp` directly by URL or the endpoint with
 * curl + token; this only controls discoverability from the landing
 * tiles.
 *
 * To widen: add emails to the set below.
 * To remove entirely: delete this file, drop the `isNlpBetaUser` calls
 * from `app/page.tsx` and `app/app-selector/page.tsx`.
 */

const NLP_BETA_EMAILS = new Set<string>([
  'martin.noble@newcastle.ac.uk',
  'martin.noble@ncl.ac.uk',
  'nmemn@newcastle.ac.uk',
]);

export function isNlpBetaUser(email: string | null | undefined): boolean {
  if (!email) return false;
  return NLP_BETA_EMAILS.has(email.toLowerCase());
}
