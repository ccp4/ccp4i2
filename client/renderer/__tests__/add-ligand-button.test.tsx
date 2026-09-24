/**
 * The "Add ligand here" button, rendered. The lookup rules are tested in
 * ligand-codes; these cover what only rendering shows: that the button is
 * disabled with a reason rather than absent, that the presence of some
 * other molecule does not enable it, and that two codes are offered rather
 * than the first being taken.
 */
import React from "react";
import { describe, it, expect, vi } from "vitest";
import { fireEvent, render, screen, waitFor } from "@testing-library/react";
import {
  AddLigandButton,
  NO_DICTIONARY_REASON,
  NO_MOLECULE_REASON,
} from "../components/moorhen/add-ligand-button";

const TARGET = { uniqueId: "/api/proxy/ccp4i2/files/42/download/" };
// Loaded through the project browser; last in the list, so the one Moorhen
// would treat as active.
const DECOY = { uniqueId: "/api/proxy/ccp4i2/files/77/download/" };

const button = () => screen.getByRole("button", { name: /add ligand here/i });

describe("AddLigandButton", () => {
  it("is disabled, and says why, with no dictionary", () => {
    render(
      <AddLigandButton ligandCodes={[]} molecules={[TARGET, DECOY]} memberCoordFileId={42} onAddLigand={vi.fn()} />
    );
    expect(button()).toBeDisabled();
    expect(screen.getByLabelText(NO_DICTIONARY_REASON)).toBeInTheDocument();
  });

  it("is disabled when only somebody else's molecule is loaded", () => {
    render(
      <AddLigandButton ligandCodes={["DRG"]} molecules={[DECOY]} memberCoordFileId={42} onAddLigand={vi.fn()} />
    );
    expect(button()).toBeDisabled();
    expect(screen.getByLabelText(NO_MOLECULE_REASON)).toBeInTheDocument();
  });

  it("places the one code directly when the tracked molecule is present", async () => {
    const onAddLigand = vi.fn(async () => {});
    render(
      <AddLigandButton
        ligandCodes={["DRG"]}
        molecules={[TARGET, DECOY]}
        memberCoordFileId={42}
        onAddLigand={onAddLigand}
      />
    );
    expect(button()).toBeEnabled();
    fireEvent.click(button());
    await waitFor(() => expect(onAddLigand).toHaveBeenCalledWith("DRG"));
    expect(onAddLigand).toHaveBeenCalledTimes(1);
  });

  it("offers a choice with two codes instead of taking the first", async () => {
    const onAddLigand = vi.fn(async () => {});
    render(
      <AddLigandButton
        ligandCodes={["F01", "F02"]}
        molecules={[TARGET, DECOY]}
        memberCoordFileId={42}
        onAddLigand={onAddLigand}
      />
    );
    fireEvent.click(button());
    expect(onAddLigand).not.toHaveBeenCalled();
    const items = await screen.findAllByRole("menuitem");
    expect(items.map((el) => el.textContent)).toEqual(["F01", "F02"]);
    fireEvent.click(screen.getByRole("menuitem", { name: "F02" }));
    await waitFor(() => expect(onAddLigand).toHaveBeenCalledWith("F02"));
    expect(onAddLigand).toHaveBeenCalledTimes(1);
  });
});
