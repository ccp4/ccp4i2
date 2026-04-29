/**
 * Compounds-frontend binding of the shared api-fetch wrapper.
 *
 * Mirrors client/renderer/api-fetch.ts but with the compounds proxy
 * prefix (``/api/proxy/compounds/``) and the optional ``X-User-Email``
 * header injection that compounds-side flows rely on as a display
 * fallback when the JWT lacks an email claim.
 *
 * Generic implementation lives in @ccp4/ccp4i2-auth's createApiFetch
 * factory; this file binds it once and re-exports the helpers under
 * their legacy names so existing call sites
 * (``import { authFetch } from "../lib/compounds/api"``) keep working
 * unchanged after the local authFetch implementations are deleted.
 */

import {
  AUTH_ERROR_EVENT,
  type AuthErrorDetail,
  createApiFetch,
  getUserEmail,
} from "@ccp4/ccp4i2-auth";

export { AUTH_ERROR_EVENT };
export type { AuthErrorDetail };

const fetcher = createApiFetch({
  baseUrl: "/api/proxy/compounds/",
  injectUserEmail: getUserEmail,
});

// `authFetch` was the legacy compounds name for the auth-injecting
// fetch wrapper. It now aliases the shared apiFetch — note that this
// version THROWS on non-OK responses (the renderer-side semantic),
// so callers' existing `if (!res.ok) throw` blocks are dead code.
// They are harmless to keep; cleanup is a separate concern.
export const authFetch = fetcher.apiFetch;

export const apiFetch = fetcher.apiFetch;
export const apiJson = fetcher.apiJson;
export const apiText = fetcher.apiText;
export const apiBlob = fetcher.apiBlob;
export const apiArrayBuffer = fetcher.apiArrayBuffer;
export const apiPost = fetcher.apiPost;
export const apiPut = fetcher.apiPut;
export const apiPatch = fetcher.apiPatch;
export const apiDelete = fetcher.apiDelete;
export const apiGet = fetcher.apiGet;
export const apiUpload = fetcher.apiUpload;
export const swrFetcher = fetcher.swrFetcher;
export const swrPostFetcher = fetcher.swrPostFetcher;
