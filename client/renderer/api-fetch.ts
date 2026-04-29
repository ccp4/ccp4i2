/**
 * CCP4i2-specific binding of the shared api-fetch wrapper.
 *
 * The generic implementation lives in @ccp4/ccp4i2-auth's
 * ``createApiFetch`` factory; this file binds it to the
 * ``/api/proxy/ccp4i2/`` proxy prefix used by the desktop renderer and
 * the cloud web app, and re-exports the bound helpers under their
 * legacy names so existing import sites
 * (``import { apiFetch, apiJson, ... } from "../api-fetch"``) keep
 * working unchanged.
 *
 * AUTH_ERROR_EVENT and AuthErrorDetail are universal (the event name is
 * the same across consumers) and so come straight from the package.
 */

import {
  AUTH_ERROR_EVENT,
  type AuthErrorDetail,
  createApiFetch,
} from "@ccp4/ccp4i2-auth";

export { AUTH_ERROR_EVENT };
export type { AuthErrorDetail };

const fetcher = createApiFetch({ baseUrl: "/api/proxy/ccp4i2/" });

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
