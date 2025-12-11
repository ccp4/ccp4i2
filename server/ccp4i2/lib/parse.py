from os import PathLike
import gemmi
import logging

logger = logging.getLogger(f"ccp4x:{__name__}")


def identify_data_type(path: str):
    path = str(path)
    for data_type_name, parse_function in (
        ("mtz", gemmi.read_mtz_file),
        ("map", gemmi.read_ccp4_map),
        ("sfcif", parse_sfcif),
        ("model", parse_structure),
        ("dictionary", parse_chemical_component),
        ("sequence", parse_pir_or_fasta),
    ):
        try:
            logger.info("Trying to parse as %s", data_type_name)
            content = parse_function(path)
            result = {"data_type_name": data_type_name, "content": content}
            logger.info("Succeeded to parse as %s", data_type_name)
            return result
        except (RuntimeError, ValueError):
            logger.info("Failed to parse as %s", data_type_name)
            continue
    raise ValueError(f"Unknown file format: {path}")


def parse(path: str):
    path = str(path)
    for parse_function in (
        gemmi.read_mtz_file,
        gemmi.read_ccp4_map,
        parse_sfcif,
        parse_structure,
        parse_chemical_component,
        parse_pir_or_fasta,
    ):
        try:
            return parse_function(path)
        except (RuntimeError, ValueError):
            continue
    raise ValueError(f"Unknown file format: {path}")


def parse_structure(path: str) -> gemmi.Structure:
    structure = gemmi.read_structure(path, format=gemmi.CoorFormat.Detect)
    if len(structure) == 0 or not any(structure[0].all()):
        raise ValueError("No models in structure.")
    return structure


def parse_sfcif(path: str) -> gemmi.ReflnBlocks:
    rblocks = gemmi.as_refln_blocks(gemmi.cif.read(path))
    if not rblocks[0]:
        raise ValueError("No reflections in structure.")
    return rblocks


def parse_pir_or_fasta(path: str) -> list[gemmi.FastaSeq]:
    with open(path, encoding="utf-8") as text:
        return gemmi.read_pir_or_fasta(text.read())


def parse_chemical_component(path: str) -> gemmi.ChemComp:
    block = gemmi.cif.read(path)[-1]
    return gemmi.make_chemcomp_from_block(block)


def parse_mtz(file_path):
    mtz = gemmi.read_mtz_file(file_path)
    reflections = []

    # Extract column labels
    column_labels = [col.label for col in mtz.columns]

    for hkl, values in zip(mtz, mtz.array):
        reflection = {
            "h": hkl.h,
            "k": hkl.k,
            "l": hkl.l,
            **{label: value for label, value in zip(column_labels, values)},
        }
        reflections.append(reflection)

    return reflections
