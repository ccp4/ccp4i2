/**
 * #613: the Validation tab's XML toggle crashed the page. The parsed
 * validation (get_validation) and the prettified XML (get_pretty_endpoint_xml)
 * fetch the same endpoint with different fetchers; under one SWR key, whichever
 * ran first filled the cache for both, and the XML editor was handed the
 * parsed object. The two must be cached under different keys.
 */
import { describe, it, expect, vi } from "vitest";

const keys: unknown[] = [];
vi.mock("swr", () => ({
  default: (key: unknown) => {
    keys.push(key);
    return { data: undefined, error: undefined, mutate: vi.fn() };
  },
}));

import { useApi } from "../api";

describe("validation SWR keys", () => {
  it("parsed validation and pretty XML do not share a cache entry", () => {
    const api = useApi();
    const ef = { type: "jobs", id: 42, endpoint: "validation" };
    api.get_validation(ef);
    api.get_pretty_endpoint_xml(ef);
    expect(keys).toHaveLength(2);
    expect(JSON.stringify(keys[0])).not.toBe(JSON.stringify(keys[1]));
    expect(keys[1]).toBe("jobs/42/validation");
  });

  it("a null endpoint still disables the validation fetch", () => {
    keys.length = 0;
    useApi().get_validation(null);
    expect(keys).toEqual([null]);
  });
});
