"""Chains of a search model that are UNK residues: built, but not sequenced.

Model building can leave a chain as UNK residues (poly-Ala or poly-Gly with
no sequence). Such a model can still place, but it is seldom what someone
choosing a search model meant: Eleanor took one for molecular replacement
and found out from the Digest afterwards (#680). The molecular-replacement
tasks warn, and do not block.
"""
import os
import tempfile

import gemmi


def unknown_residue_chains(path):
    """[(chain, n_unk, n_residues)] for each polymer chain with UNK residues."""
    structure = gemmi.read_structure(str(path), format=gemmi.CoorFormat.Detect)
    if len(structure) == 0:
        return []
    found = []
    for chain in structure[0]:
        polymer = chain.get_polymer()
        if len(polymer) == 0:
            continue
        unk = sum(1 for residue in polymer if residue.name == "UNK")
        if unk:
            found.append((chain.name, unk, len(polymer)))
    return found


def search_model_unknown_residues(pdb_file):
    """unknown_residue_chains of a CPdbDataFile, as its selection leaves it."""
    if not pdb_file.isSet():
        return []
    if not pdb_file.isSelectionSet():
        return unknown_residue_chains(pdb_file.getFullPath())
    with tempfile.TemporaryDirectory() as work:
        return unknown_residue_chains(pdb_file.getSelectedAtomsFile("selected", work))


def describe(chains):
    """'chain A has 120 UNK of 130 residues; chain B ...'"""
    return "; ".join(f"chain {name} has {unk} UNK of {total} residues"
                     for name, unk, total in chains)


def warn_unknown_residues(error, task_name, models, code=180):
    """Append a warning to *error* for each search model with UNK chains.

    models: [(parameter path under container, CPdbDataFile)], for example
    ("inputData.XYZIN", inp.XYZIN). A model that cannot be read is left to
    the checks that own it.
    """
    from ccp4i2.core import CCP4ErrorHandling

    for param, pdb_file in models:
        try:
            chains = search_model_unknown_residues(pdb_file)
        except Exception:
            continue
        if chains:
            error.append(
                klass=task_name, code=code,
                details=(f"Search model {os.path.basename(str(pdb_file.getFullPath()))}: "
                         f"{describe(chains)}. UNK residues are built but not sequenced; "
                         "is this the model you meant to search with?"),
                name=f"{task_name}.container.{param}",
                severity=CCP4ErrorHandling.SEVERITY_WARNING)


def ensemble_models(input_data):
    """(parameter path, CPdbDataFile) for each structure in Phaser's ENSEMBLES."""
    models = []
    for i, ensemble in enumerate(getattr(input_data, "ENSEMBLES", None) or []):
        for j, item in enumerate(ensemble.pdbItemList):
            models.append((f"inputData.ENSEMBLES[{i}].pdbItemList[{j}].structure",
                           item.structure))
    return models
