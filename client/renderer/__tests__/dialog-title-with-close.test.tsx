/**
 * Every viewing dialog carries a visible X (#675). Esc and a backdrop click
 * close a dialog too, but neither can be seen, and users asked how to get
 * out of the file preview.
 */
import React from "react";
import { describe, it, expect, vi } from "vitest";
import { render, screen, fireEvent } from "@testing-library/react";
import { Dialog } from "@mui/material";

import { DialogTitleWithClose } from "../components/dialog-title-with-close";

describe("DialogTitleWithClose", () => {
  it("shows the title and a Close button that calls onClose", () => {
    const onClose = vi.fn();
    render(
      <Dialog open onClose={onClose}>
        <DialogTitleWithClose onClose={onClose}>gamma_native.mtz</DialogTitleWithClose>
      </Dialog>
    );

    expect(screen.getByText("gamma_native.mtz")).toBeTruthy();
    fireEvent.click(screen.getByRole("button", { name: "Close" }));
    expect(onClose).toHaveBeenCalledTimes(1);
  });
});
