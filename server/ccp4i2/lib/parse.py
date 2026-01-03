import logging

import gemmi

logger = logging.getLogger(f"ccp4i2:{__name__}")


def identify_data_type(path: str) -> str:
    path = str(path)
    for data_type_name, parse_function in (
        ("mtz", gemmi.read_mtz_file),
        ("map", gemmi.read_ccp4_map),
        ("sfcif", _parse_sfcif),
        ("model", _parse_structure),
        ("dictionary", _parse_chemical_component),
        ("sequence", _parse_pir_or_fasta),
    ):
        try:
            logger.info("Trying to parse as %s", data_type_name)
            parse_function(path)
            logger.info("Succeeded to parse as %s", data_type_name)
            return data_type_name
        except (RuntimeError, ValueError):
            logger.info("Failed to parse as %s", data_type_name)
            continue
    raise ValueError(f"Unknown file format: {path}")


def _parse_structure(path: str):
    structure = gemmi.read_structure(path, format=gemmi.CoorFormat.Detect)
    if len(structure) == 0 or not any(structure[0].all()):
        raise ValueError("No models in structure.")


def _parse_sfcif(path: str):
    rblocks = gemmi.as_refln_blocks(gemmi.cif.read(path))
    if not rblocks[0]:
        raise ValueError("No reflections in structure.")


def _parse_pir_or_fasta(path: str):
    with open(path, encoding="utf-8") as text:
        gemmi.read_pir_or_fasta(text.read())


def _parse_chemical_component(path: str):
    block = gemmi.cif.read(path)[-1]
    gemmi.make_chemcomp_from_block(block)
