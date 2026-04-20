import "@testing-library/jest-dom/vitest";
import { afterEach, vi } from "vitest";
import { cleanup } from "@testing-library/react";

// next/font/google is not executable in vitest — stub every font factory
// because the theme provider module is transitively imported by utils → api.
vi.mock("next/font/google", () => {
  const fontFactory = () => ({
    className: "",
    style: { fontFamily: "sans-serif" },
    variable: "--font-stub",
  });
  return new Proxy(
    { __esModule: true },
    {
      get: (target: any, key: string) => (key in target ? target[key] : fontFactory),
    }
  );
});

// next/navigation: imported transitively by auth/token utilities.
vi.mock("next/navigation", () => ({
  useRouter: () => ({
    push: () => {},
    replace: () => {},
    back: () => {},
    forward: () => {},
    refresh: () => {},
    prefetch: () => {},
  }),
  usePathname: () => "/",
  useSearchParams: () => new URLSearchParams(),
  useParams: () => ({}),
  redirect: () => {},
  notFound: () => {},
}));

// Stub network calls globally — there is no backend in the test env. Returns a
// minimal Response so SWR fetchers resolve cleanly instead of retrying with
// noisy AbortSignal errors from undici.
const stubResponse = () =>
  Promise.resolve(
    new Response("{}", {
      status: 200,
      headers: { "Content-Type": "application/json" },
    })
  );
vi.stubGlobal("fetch", stubResponse);

afterEach(() => {
  cleanup();
});
