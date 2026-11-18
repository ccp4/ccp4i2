/**
 * The atom-selection builder starts its chain boxes from the job's selection
 * when that selection names whole chains only (#876), so a job reopened with
 * "A/" shows chain A ticked rather than every chain.
 */
import { describe, expect, it } from "vitest";
import { chainsFromSelection } from "../lib/atom-selection";

const chains = ["A", "B", "C"];

describe("chainsFromSelection", () => {
  it.each([
    ["A", ["A"]],
    ["A/", ["A"]],
    ["  B/ ", ["B"]],
    ["{A or B}", ["A", "B"]],
    ["{A/ or C/}", ["A", "C"]],
    ["A or B", ["A", "B"]],
  ])("reads %j as chains %j", (selection, expected) => {
    expect([...(chainsFromSelection(selection, chains) ?? [])]).toEqual(expected);
  });

  it.each([
    [""],
    ["A/1-50"],
    ["A and protein"],
    ["{A} or A/401(ATP)"],
    ["Z/"],
    ["not A"],
  ])("leaves %j to be shown as typed", (selection) => {
    expect(chainsFromSelection(selection, chains)).toBeNull();
  });
});
