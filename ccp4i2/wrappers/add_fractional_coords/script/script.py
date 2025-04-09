import argparse
import os

import gemmi
import pytest


def parse_args(args):
    description = (
        "Writes an mmCIF file that includes fractional coordinates. "
        "The input is a coordinate file "
        "in mmCIF, mmJSON or PDB format that may be gzipped. "
        "The output mmCIF file contains fractional coordinates in the "
        "_atom_site.fract_x, _atom_site.fract_y and _atom_site.fract_z fields."
    )
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("xyzin", help="input coordinate file")
    parser.add_argument("xyzout", help="output mmCIF file")
    return parser.parse_args(args)


def read_doc_and_structure(path):
    try:
        doc = gemmi.cif.read(path)
        structure = gemmi.make_structure_from_block(doc[0])
    except (RuntimeError, ValueError):
        try:
            doc = gemmi.cif.read_mmjson(path)
            structure = gemmi.make_structure_from_block(doc[0])
        except (RuntimeError, ValueError):
            structure = gemmi.read_structure(path)
            doc = structure.make_mmcif_document()
    return doc, structure


def add_fractional_coordinates(doc, structure):
    block = doc[0]
    category = block.get_mmcif_category("_atom_site.")
    category["fract_x"] = []
    category["fract_y"] = []
    category["fract_z"] = []
    for x, y, z in zip(category["Cartn_x"], category["Cartn_y"], category["Cartn_z"]):
        orthogonal = gemmi.Position(float(x), float(y), float(z))
        fractional = structure.cell.fractionalize(orthogonal)
        category["fract_x"].append(round(fractional.x, 3))
        category["fract_y"].append(round(fractional.y, 3))
        category["fract_z"].append(round(fractional.z, 3))
    block.set_mmcif_category("_atom_site.", category)


@pytest.mark.parametrize(
    "path",
    [
        os.path.join(os.environ["CCP4"], "examples", "data", "insulin.pdb"),
        os.path.join(os.environ["CCP4"], "examples", "toxd", "toxd.cif"),
    ],
)
def test(path):
    doc, structure = read_doc_and_structure(path)
    add_fractional_coordinates(doc, structure)
    category = doc[0].get_mmcif_category("_atom_site.")
    assert all(key in category for key in ("fract_x", "fract_y", "fract_z"))


def main(args=None):
    args = parse_args(args)
    doc, structure = read_doc_and_structure(args.xyzin)
    add_fractional_coordinates(doc, structure)
    doc.write_file(args.xyzout, gemmi.cif.Style.Aligned)


if __name__ == "__main__":
    main()
