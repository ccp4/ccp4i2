/**
 * A directory the user made in advance is accepted while it is empty (#686).
 * Qt-era users habitually make the folder first and were refused for it.
 */
import { describe, it, expect } from "vitest";
import { describeExistingPath } from "../components/new-project-content";

describe("describeExistingPath", () => {
  it("allows a path nothing is at", () => {
    expect(
      describeExistingPath({
        exists: false,
        isDirectory: false,
        isEmpty: false,
      }),
    ).toBe("");
  });

  it("allows an empty directory made in advance", () => {
    expect(
      describeExistingPath({ exists: true, isDirectory: true, isEmpty: true }),
    ).toBe("");
  });

  it("refuses a directory with contents, saying so", () => {
    expect(
      describeExistingPath({ exists: true, isDirectory: true, isEmpty: false }),
    ).toBe("Directory already exists and is not empty");
  });

  it("refuses a plain file of that name", () => {
    expect(
      describeExistingPath({
        exists: true,
        isDirectory: false,
        isEmpty: false,
      }),
    ).toBe("A file with this name already exists");
  });
});
