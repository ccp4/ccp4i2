/**
 * An unset numeric parameter arrives from the server with a sentinel _value
 * (0) and _valueState "NOT_SET". The text field must render it blank, show
 * the qualifier default as placeholder text, and send null when the user
 * clears a value so the server returns the parameter to unset.
 */
import React from "react";
import { describe, it, expect, vi, beforeEach } from "vitest";
import { render, screen, fireEvent, waitFor } from "@testing-library/react";

const mocks = vi.hoisted(() => ({
  items: {} as Record<string, any>,
  setParameter: vi.fn(async () => ({
    success: true,
    data: { updated_item: {} },
  })),
}));

vi.mock("../utils", () => ({
  useJob: () => ({
    useTaskItem: (name: string) => ({ item: mocks.items[name] }),
    getValidationColor: () => undefined,
    setParameter: mocks.setParameter,
    setParameterNoMutate: mocks.setParameter,
  }),
  valueOfItem: (item: any) => item?._value,
}));
vi.mock("../providers/task-provider", () => ({
  useTaskInterface: () => ({ inFlight: false, setInFlight: () => {} }),
}));
vi.mock("../providers/popcorn-provider", () => ({
  usePopcorn: () => ({ setMessage: () => {} }),
}));
vi.mock("../components/task/task-elements/error-info", () => ({
  ErrorTrigger: () => null,
}));

import { CSimpleElement } from "../components/task/task-elements/csimple";

const job: any = { id: 1, status: 1 };
const PATH = "task.container.controlParameters.FRAC";

const floatItem = (overrides: Record<string, any>) => ({
  _class: "CFloat",
  _baseClass: "CFloat",
  _value: 0,
  _qualifiers: { guiLabel: "Fraction", default: 0.05, min: 0, max: 1 },
  _CONTENTS_ORDER: [],
  _objectPath: PATH,
  ...overrides,
});

const renderField = () =>
  render(<CSimpleElement itemName="FRAC" job={job} type="float" />);

describe("CSimpleElement for an unset CFloat", () => {
  beforeEach(() => {
    mocks.setParameter.mockClear();
  });

  it("renders blank with the default as placeholder", () => {
    mocks.items.FRAC = floatItem({ _valueState: "NOT_SET" });
    renderField();
    const input = screen.getByRole("textbox");
    expect(input).toHaveValue("");
    expect(input).toHaveAttribute("placeholder", "0.05");
  });

  it("renders an explicitly set zero as 0", () => {
    mocks.items.FRAC = floatItem({ _valueState: "EXPLICITLY_SET" });
    renderField();
    expect(screen.getByRole("textbox")).toHaveValue("0");
  });

  it("does not write when an untouched unset field loses focus", () => {
    mocks.items.FRAC = floatItem({ _valueState: "NOT_SET" });
    renderField();
    fireEvent.blur(screen.getByRole("textbox"));
    expect(mocks.setParameter).not.toHaveBeenCalled();
  });

  it("sends null when a set value is cleared", async () => {
    mocks.items.FRAC = floatItem({ _value: 0.3, _valueState: "EXPLICITLY_SET" });
    renderField();
    const input = screen.getByRole("textbox");
    expect(input).toHaveValue("0.3");
    fireEvent.change(input, { target: { value: "" } });
    fireEvent.blur(input);
    await waitFor(() =>
      expect(mocks.setParameter).toHaveBeenCalledWith({
        object_path: PATH,
        value: null,
      })
    );
  });
});
