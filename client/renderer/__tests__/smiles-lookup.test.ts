import { describe, expect, it, vi } from "vitest";
import { fetchSmilesBatch } from "../lib/smiles-lookup";

describe("fetchSmilesBatch", () => {
  it("asks for one batch with a page size that fits it, and follows next pages", async () => {
    const fetchJson = vi.fn()
      .mockResolvedValueOnce({
        results: [{ reg_number: 101, smiles: "CCO" }],
        next: "https://materia.example/api/compounds/compounds/?page=2&reg_number__in=101,102,103",
      })
      .mockResolvedValueOnce({ results: [{ reg_number: 103, smiles: "c1ccccc1" }], next: null });
    const found = await fetchSmilesBatch([101, 102, 103], fetchJson);
    expect(fetchJson.mock.calls[0][0]).toBe(
      "/api/proxy/compounds/compounds/?reg_number__in=101,102,103&page_size=3");
    expect(fetchJson.mock.calls[1][0]).toBe(
      "/api/proxy/compounds/compounds/?page=2&reg_number__in=101,102,103");
    expect(Object.fromEntries(found)).toEqual({ 101: "CCO", 102: null, 103: "c1ccccc1" });
  });

  it("takes an unpaginated list as it comes", async () => {
    const found = await fetchSmilesBatch([7], vi.fn().mockResolvedValue([{ reg_number: 7, smiles: "N" }]));
    expect(found.get(7)).toBe("N");
  });
});
