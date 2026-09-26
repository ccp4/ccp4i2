import json
import logging
import gemmi
from typing import Dict, Type

from ccp4i2.core import CCP4File
from ccp4i2.core import CCP4XtalData
from ccp4i2.core import CCP4ModelData
from ccp4i2.core.CCP4Container import CContainer
from ccp4i2.core.base_object.cdata_file import CDataFile
from ccp4i2.core.base_object.cdata import CData
from ccp4i2.core.CCP4XtalData import CGenericReflDataFile, CMapDataFile, CMtzDataFile
from ccp4i2.core.CCP4ModelData import CPdbDataFile, CDictDataFile, CAsuDataFile
# Import stub class for isinstance checks - subclasses like CObsDataFile inherit from
# stubs (CMtzDataFile) not implementations (CMtzDataFile)
from ccp4i2.core.CCP4XtalData import CMtzDataFile, CUnmergedDataFile
# mmcifutils imported lazily in digest_cgenericrefldatafile_file_object to avoid
# numpy dependency at module load time (numpy in py-packages may be incompatible)
from ..containers.find_objects import find_objects
from ..containers.get_container import get_job_container
from ..containers.json_encoder import CCP4i2JsonEncoder
from ..plugins.plugin_context import get_plugin_with_context
from ..formats.cif_ligand import parse_cif_ligand_summary, extract_monomer_atoms_bonds, extract_all_monomers_atoms_bonds, generate_all_molblocks
from ..parameters.value_dict import value_dict_for_object
from ....db import models
from ...parse import identify_data_type
from ...json_safety import replace_non_finite

logger = logging.getLogger(f"ccp4i2:{__name__}")


def normalize_object_path(object_path: str) -> str:
    """
    Normalize object paths from frontend to match backend container structure.

    The frontend JSON encoder includes the full hierarchy path which includes
    `.container.` (e.g., `prosmart_refmac.container.inputData.XYZIN`), but the
    backend container structure doesn't have that extra level.

    This function strips the `.container.` segment if present after the task name.

    Args:
        object_path: Path like "prosmart_refmac.container.inputData.XYZIN"

    Returns:
        Normalized path like "prosmart_refmac.inputData.XYZIN"
    """
    # Split into parts
    parts = object_path.split('.')

    # If second element is 'container', remove it
    # e.g., ['prosmart_refmac', 'container', 'inputData', 'XYZIN']
    #    -> ['prosmart_refmac', 'inputData', 'XYZIN']
    if len(parts) >= 2 and parts[1] == 'container':
        parts = [parts[0]] + parts[2:]

    return '.'.join(parts)


def _class_named(full_class_name):
    """A file class by bare name, from the one CData registry (core/cdata_registry.py)."""
    from ccp4i2.core.cdata_registry import cdata_classes
    return cdata_classes().get(full_class_name)


# The mimetype <-> class tables are the ones the database is seeded from; this
# module used to carry its own copy, which drifted (it lacked text/plain and
# the PhaserTNG DAG type, so files of those types could not be digested).
from ccp4i2.db.ccp4i2_static_data import FILETYPES_CLASS, FILETYPES_TEXT


def is_basic_type(obj):
    return isinstance(obj, (str, int, float, bool, type(None)))


def flatten_instance(obj):
    if is_basic_type(obj):
        return obj
    elif isinstance(obj, list):
        return [flatten_instance(item) for item in obj]
    elif isinstance(obj, tuple):
        return tuple(flatten_instance(item) for item in obj)
    elif isinstance(obj, dict):
        return {key: flatten_instance(value) for key, value in obj.items()}
    elif hasattr(obj, "__dict__"):
        return {
            key: flatten_instance(value)
            for key, value in vars(obj).items()
            if not callable(value) and not key.startswith("_")
        }
    elif hasattr(obj, "__slots__"):
        return {
            slot: flatten_instance(getattr(obj, slot))
            for slot in obj.__slots__
            if hasattr(obj, slot)
        }
    else:
        # fallback: use repr to avoid exceptions
        return repr(obj)


def digest_file(the_file: models.File):
    mimetype = the_file.type.pk
    if mimetype not in FILETYPES_TEXT:
        return {"status": "Failed", "reason": "File type not supported for digest"}

    mimetype_index = FILETYPES_TEXT.index(mimetype)
    logger.debug("mimetype_index %s", mimetype_index)
    if mimetype_index >= len(FILETYPES_CLASS):
        return {"status": "Failed", "reason": "File type not supported for digest"}
    class_name = FILETYPES_CLASS[mimetype_index]
    logger.debug("class_name %s", class_name)

    # Use dynamic class registry to find the class
    full_class_name = f"C{class_name}"
    the_class = _class_named(full_class_name)

    if the_class is None:
        return {"status": "Failed", "reason": f"File type class not found: {full_class_name}"}
    logger.debug("the_class %s", the_class)
    try:
        file_object = the_class()
        file_object.setFullPath(str(the_file.path))
    except Exception as err:
        logger.exception("Error creating file object %s", the_file, exc_info=err)
        return {"status": "Failed", "reason": str(err), "digest": {}}
    return digest_file_object(file_object)


def digest_param_file(the_job, object_path):
    # Use plugin context for consistent container access (same as set_param/get_param)
    plugin_result = get_plugin_with_context(the_job)
    if not plugin_result.success:
        return {"status": "Failed", "reason": plugin_result.error, "digest": {}}

    plugin = plugin_result.data

    # Normalize path to strip .container. segment if present from frontend
    normalized_path = normalize_object_path(object_path)

    try:
        file_object: CDataFile = plugin.container.find_by_path(normalized_path, skip_first=True)
        return digest_file_object(file_object)
    except IndexError as err:
        logger.exception("Error finding object with path %s (normalized: %s)", object_path, normalized_path, exc_info=err)
        return {"status": "Failed", "reason": str(err), "digest": {}}
    except Exception as err:
        logger.exception("Other exception %s (normalized: %s)", object_path, normalized_path, exc_info=err)
        return {"status": "Failed", "reason": str(err), "digest": {}}


def json_safe(obj):
    """Replace non-finite floats (NaN, +/-inf) with None, recursively.

    A digest is sent as JSON, which has no NaN: one in the payload makes the
    renderer raise and the endpoint answer 500 instead of a digest. They do
    occur in real files -- gemmi's CifToMtz writes a NaN dataset wavelength
    when the structure-factor mmCIF records none -- and NaN is truthy, so an
    ``if value:`` guard does not keep it out.

    The implementation now lives in `ccp4i2.lib.json_safety`, which KPI
    handling and the JSON renderer share; this name is kept for its callers.
    """
    return replace_non_finite(obj)


def digest_file_object(file_object: CDataFile):
    return json_safe(_digest_file_object(file_object))


def _digest_file_object(file_object: CDataFile):
    if not isinstance(file_object, CCP4File.CDataFile):
        return {"status": "Failed", "reason": "Not a valid file object", "digest": {}}
    if not file_object.isSet():
        return {"status": "Failed", "reason": "File object is not set", "digest": {}}
    if isinstance(file_object, CCP4ModelData.CPdbDataFile):
        return digest_cpdbdata_file_object(file_object)
    if isinstance(file_object, CCP4XtalData.CGenericReflDataFile):
        return digest_cgenericrefldatafile_file_object(file_object)
    # CMtzDataFile inherits from CDataFile, not CGenericReflDataFile, so check separately
    # Use CMtzDataFile for isinstance check because subclasses inherit from stubs
    if isinstance(file_object, CMtzDataFile):
        return digest_cmtzdatafile_file_object(file_object)
    if isinstance(file_object, CMapDataFile):
        return digest_cmapdatafile_file_object(file_object)
    if isinstance(file_object, CCP4ModelData.CSeqAlignDataFile):
        return digest_cseqaligndata_file_object(file_object)
    if isinstance(file_object, CCP4ModelData.CSeqDataFile):
        return digest_cseqdata_file_object(file_object)
    if isinstance(file_object, (CCP4ModelData.CDictDataFile, CDictDataFile)):
        return digest_cdictdata_file_object(file_object)
    if isinstance(file_object, (CCP4ModelData.CAsuDataFile, CAsuDataFile)):
        return digest_casudatafile_file_object(file_object)
    # CUnmergedDataFile can hold MTZ, mmCIF, SCA, XDS etc.  For mmCIF
    # files, route through the generic refl digest to get rblock_infos.
    # For MTZ, delegate to the MTZ digest.  Others fall through.
    if isinstance(file_object, CUnmergedDataFile):
        path = str(file_object.fullPath) if file_object.fullPath else ""
        ext = path.rsplit(".", 1)[-1].lower() if path else ""
        if ext in ("cif", "mmcif", "ent"):
            generic_refl = CCP4XtalData.CGenericReflDataFile()
            generic_refl.setFullPath(path)
            return digest_cgenericrefldatafile_file_object(generic_refl)
        if ext == "mtz":
            return digest_cmtzdatafile_file_object(file_object)
    if type(file_object) is CCP4File.CDataFile:
        return digest_cdatafile_file_object(file_object)
    return digest_other_file_object(file_object)


def digest_other_file_object(file_object: CDataFile):
    try:
        file_object.loadFile()
        file_object.setContentFlag()
        contents = file_object.getFileContent()
        content_dict = value_dict_for_object(contents)

        # If this happens to be an mmCIF reflection file (eg a
        # CUnmergedDataFile pointing at a .cif/.mmcif), enrich the
        # digest with per-block info so UIs (aimless_pipe etc) can
        # offer a block selector.
        full_path = None
        try:
            full_path = str(file_object.fullPath) if file_object.fullPath else None
        except Exception:
            full_path = None
        if full_path and full_path.lower().endswith((".cif", ".mmcif", ".ent")):
            try:
                content_dict["rblock_infos"] = _mmcif_block_infos(full_path)
                content_dict["format"] = "mmcif"
            except Exception as err:
                logger.warning("Failed to scan mmcif blocks for %s: %s", full_path, err)
        return content_dict
    except Exception as err:
        logger.exception("Error digesting file %s", file_object, exc_info=err)
        return {
            "status": "Failed",
            "reason": f"Failed digesting CDataFile {err}",
            "digest": {},
        }


def _mmcif_block_infos(full_path: str):
    """Scan every reflection block in an mmCIF file and return a list
    of dicts describing it (cell, spacegroup, columns, FreeR status,
    and a `suitableForMerging` flag identifying unmerged-intensity
    blocks suitable as aimless input)."""
    from ccp4i2.pipelines.import_merged.script import mmcifutils

    mmcif = gemmi.cif.read_file(full_path)
    rblocks = gemmi.as_refln_blocks(mmcif)
    rblock_infos = []
    for rb in rblocks:
        blkinfo = mmcifutils.CifBlockInfo(rb)
        block_dict = flatten_instance(blkinfo)

        block_dict["hasFreeR"] = blkinfo.haveFreeR()
        block_dict["freerValid"] = blkinfo.validFreeR() or False
        block_dict["freerWarnings"] = blkinfo.freerWarning() or []

        if hasattr(blkinfo, "labelsets") and blkinfo.labelsets:
            block_dict["typeCodes"] = blkinfo.labelsets.getTypeCodes()
            block_dict["columnSetsText"] = blkinfo.labelsets.columnsetstext()

        suitable = False
        if getattr(blkinfo, "unmerged", False):
            try:
                status, _msg = blkinfo.columnsOK(formerged=False)
                suitable = status >= 0
            except Exception:
                suitable = False
        block_dict["suitableForMerging"] = suitable

        rblock_infos.append(block_dict)
    return rblock_infos


def _extract_chain_sequences(gemmi_structure):
    """Extract per-chain polymer sequences from a gemmi Structure.

    Returns a dict mapping chain ID to one-letter sequence string,
    for protein and nucleic acid polymer chains only.
    """
    import gemmi

    sequences = {}
    if len(gemmi_structure) == 0:
        return sequences

    model = gemmi_structure[0]
    for chain in model:
        polymer = chain.get_polymer()
        if polymer and len(polymer) > 0:
            seq = gemmi.one_letter_code([res.name for res in polymer])
            if seq:
                sequences[chain.name] = seq
    return sequences


def _extract_ligands(gemmi_structure):
    """Extract non-polymer, non-water residues (ligands and ions) from a gemmi
    Structure.

    Returns a list of dicts with the residue code, its chain, sequence number,
    and a Moorhen/coot CID (`/model/chain/seqnum(NAME)`) — the form a scene
    selection uses. This lets a scene author (a person, or an LLM) reference a
    ligand precisely instead of guessing where it sits. Classification uses
    gemmi's residue tables; anything not recognised as an amino acid, nucleic
    acid, or water is treated as a ligand (so novel compounds are included).
    """
    import gemmi

    ligands = []
    if len(gemmi_structure) == 0:
        return ligands

    model = gemmi_structure[0]
    for chain in model:
        for res in chain:
            info = gemmi.find_tabulated_residue(res.name)
            if info is not None and (
                info.is_amino_acid() or info.is_nucleic_acid() or info.is_water()
            ):
                continue
            seq_num = res.seqid.num
            ligands.append({
                "chainId": chain.name,
                "name": res.name,
                "seqNum": seq_num,
                "cid": f"/1/{chain.name}/{seq_num}({res.name})",
                "nAtoms": len(res),
            })
    return ligands


def digest_cpdbdata_file_object(file_object: CPdbDataFile):
    content_dict = {}
    if not isinstance(file_object, CCP4ModelData.CPdbDataFile):
        return {"status": "Failed", "reason": "Not a CPdbDataFile object", "digest": {}}
    if not file_object.isSet():
        return {"status": "Failed", "reason": "File object is not set", "digest": {}}
    try:
        file_object.loadFile()
        file_object.setContentFlag()
        contents = file_object.getFileContent()
        content_dict = value_dict_for_object(contents)
        # If value_dict returns None, return empty dict
        if content_dict is None:
            content_dict = {}

        # Explicitly extract composition and sequences from the content object.
        # These are @property accessors backed by underscore attributes (_composition,
        # _gemmi_structure) which value_dict_for_object may miss if earlier strategies
        # in handle_cdata populate the result before the known_content_props check.
        if contents is not None:
            composition = getattr(contents, 'composition', None)
            if composition is not None and 'composition' not in content_dict:
                content_dict['composition'] = value_dict_for_object(composition)

            gemmi_struct = getattr(contents, '_gemmi_structure', None)
            if gemmi_struct is not None and 'sequences' not in content_dict:
                content_dict['sequences'] = _extract_chain_sequences(gemmi_struct)
            if gemmi_struct is not None and 'ligands' not in content_dict:
                content_dict['ligands'] = _extract_ligands(gemmi_struct)

        return content_dict
    except Exception as err:
        logger.exception("Error digesting file %s", file_object, exc_info=err)
        return {
            "status": "Failed",
            "reason": f"Failed digesting CPdbDataFile {err}",
            "digest": {},
        }


_MAP_SUBTYPE_LABELS = {
    1: "Normal (electron density)",
    2: "Difference (Fo-Fc)",
    3: "Anomalous difference",
    4: "Mask",
    5: "Half map",  # CMapDataFile.SUBTYPE_HALFMAP (cryo-EM); label kept forward-compatible
}


def digest_cmapdatafile_file_object(file_object):
    """Digest a real-space CCP4/MRC map by reading its header (gemmi only).

    Reports the grid, full-cell sampling (MX/MY/MZ) and voxel spacing, unit cell,
    axis order, start (nxstart...), spacegroup, density statistics and the
    CMapDataFile subType -- the header information a map file otherwise never
    surfaces in the file preview. Reads only the header (the DMIN/DMAX/DMEAN/RMS
    stats come from the header words), so it does not load the whole grid.
    """
    if not isinstance(file_object, CMapDataFile):
        return {"status": "Failed", "reason": "Not a CMapDataFile object", "digest": {}}
    if not file_object.isSet():
        return {"status": "Failed", "reason": "File object is not set", "digest": {}}
    try:
        import gemmi

        path = str(file_object.fullPath)
        m = gemmi.read_ccp4_map(path)  # header only; no setup() so nothing expands

        cell = [round(m.header_float(i), 4) for i in range(11, 17)]
        mx, my, mz = m.header_i32(8), m.header_i32(9), m.header_i32(10)
        spacing = [
            round(cell[0] / mx, 4) if mx else None,
            round(cell[1] / my, 4) if my else None,
            round(cell[2] / mz, 4) if mz else None,
        ]
        angles_90 = all(abs(a - 90.0) < 1e-3 for a in cell[3:6])
        ispg = m.header_i32(23)

        sub_type = None
        try:
            if file_object.subType.isSet():
                sub_type = int(file_object.subType)
        except Exception:
            sub_type = None

        return {
            "format": "CCP4/MRC map",
            "mode": m.header_i32(4),
            "grid_sampling": [mx, my, mz],
            "grid_stored": [m.header_i32(1), m.header_i32(2), m.header_i32(3)],
            "start": [m.header_i32(5), m.header_i32(6), m.header_i32(7)],
            "axis_order": [m.header_i32(17), m.header_i32(18), m.header_i32(19)],
            "cell": {"a": cell[0], "b": cell[1], "c": cell[2],
                     "alpha": cell[3], "beta": cell[4], "gamma": cell[5]},
            "spacing": spacing,
            "spacegroup": ispg,
            "statistics": {
                "min": round(m.header_float(20), 5),
                "max": round(m.header_float(21), 5),
                "mean": round(m.header_float(22), 5),
                "rms": round(m.header_float(55), 5),
            },
            "sub_type": sub_type,
            "sub_type_label": _MAP_SUBTYPE_LABELS.get(sub_type),
            "likely_em": ispg in (0, 1) and angles_90,
        }
    except Exception as err:
        logger.exception("Error digesting map file %s", file_object, exc_info=err)
        return {
            "status": "Failed",
            "reason": f"Failed digesting CMapDataFile {err}",
            "digest": {},
        }


def digest_cseqdata_file_object(file_object: CPdbDataFile):
    content_dict = {}
    if not isinstance(file_object, CCP4ModelData.CSeqDataFile):
        return {"status": "Failed", "reason": "Not a CSeqDataFile object", "digest": {}}
    if not file_object.isSet():
        return {"status": "Failed", "reason": "File object is not set", "digest": {}}
    try:
        file_object.loadFile()
        file_object.setContentFlag()
        contents = file_object.getFileContent()
        content_dict = value_dict_for_object(contents)
        return content_dict
    except Exception as err:
        logger.exception("Error digesting file %s", file_object, exc_info=err)
        return {
            "status": "Failed",
            "reason": f"Failed digesting CSeqDataFile {err}",
            "digest": {},
        }


def digest_cseqaligndata_file_object(file_object):
    """
    Digest a CSeqAlignDataFile to provide sequence identifiers.

    Uses identifyFile() to detect the alignment format and extract
    the list of sequence IDs, which can be used to populate TARGETINDEX
    enumerators in task interfaces (chainsaw, sculptor).

    Returns a dict with:
    - format: Detected alignment format (clustal, fasta, phylip, stockholm)
    - identifiers: List of sequence ID strings from the alignment
    """
    if not isinstance(file_object, CCP4ModelData.CSeqAlignDataFile):
        return {"status": "Failed", "reason": "Not a CSeqAlignDataFile object", "digest": {}}
    if not file_object.isSet():
        return {"status": "Failed", "reason": "File object is not set", "digest": {}}
    try:
        fmt, id_list = file_object.identifyFile()
        return {
            "format": fmt,
            "identifiers": id_list,
        }
    except Exception as err:
        logger.exception("Error digesting CSeqAlignDataFile %s", file_object, exc_info=err)
        return {
            "status": "Failed",
            "reason": f"Failed digesting CSeqAlignDataFile: {err}",
            "digest": {},
        }


def digest_casudatafile_file_object(file_object: CAsuDataFile):
    """
    Digest a CAsuDataFile to provide sequence selection information.

    Returns a dict with:
    - sequences: List of sequence entries with index, name, polymerType, nCopies, selected
    - This allows the frontend to render checkboxes for sequence selection
    """
    if not isinstance(file_object, (CCP4ModelData.CAsuDataFile, CAsuDataFile)):
        return {"status": "Failed", "reason": "Not a CAsuDataFile object", "digest": {}}
    if not file_object.isSet():
        return {"status": "Failed", "reason": "File object is not set", "digest": {}}
    try:
        file_object.loadFile()
        contents = file_object.getFileContent()  # CAsuContent

        if contents is None:
            return {"status": "Failed", "reason": "No file content loaded", "digest": {}}

        # Get the selection CDict from the file object
        selection = file_object.selection if hasattr(file_object, 'selection') else None

        # Build sequence list with selection state
        sequences = []
        if hasattr(contents, 'seqList') and contents.seqList is not None:
            for idx, seq in enumerate(contents.seqList):
                # Get sequence name
                name = ""
                if hasattr(seq, 'name') and seq.name is not None:
                    if hasattr(seq.name, 'value') and seq.name.isSet():
                        name = str(seq.name.value)
                    elif seq.isSet('name'):
                        name = str(seq.name)
                if not name:
                    name = f"Sequence_{idx}"

                # Get polymer type
                polymer_type = "PROTEIN"
                if hasattr(seq, 'polymerType') and seq.polymerType is not None:
                    if hasattr(seq.polymerType, 'value') and seq.polymerType.isSet():
                        polymer_type = str(seq.polymerType.value)
                    elif seq.isSet('polymerType'):
                        polymer_type = str(seq.polymerType)

                # Get nCopies
                n_copies = 1
                if hasattr(seq, 'nCopies') and seq.nCopies is not None:
                    if hasattr(seq.nCopies, 'value') and seq.nCopies.isSet():
                        n_copies = int(seq.nCopies.value)
                    elif seq.isSet('nCopies'):
                        n_copies = int(seq.nCopies)

                # Get sequence string (truncated for display)
                sequence_str = ""
                if hasattr(seq, 'sequence') and seq.sequence is not None:
                    if hasattr(seq.sequence, 'value') and seq.sequence.isSet():
                        sequence_str = str(seq.sequence.value)
                    elif seq.isSet('sequence'):
                        sequence_str = str(seq.sequence)

                # Get description
                description = ""
                if hasattr(seq, 'description') and seq.description is not None:
                    if hasattr(seq.description, 'value') and seq.description.isSet():
                        description = str(seq.description.value)
                    elif seq.isSet('description'):
                        description = str(seq.description)

                # Determine selection state from CDict
                # Default is True (selected) if not explicitly set
                selected = True
                if selection is not None and hasattr(selection, 'get'):
                    selected = selection.get(name, True)

                sequences.append({
                    "index": idx,
                    "name": name,
                    "polymerType": polymer_type,
                    "nCopies": n_copies,
                    "sequenceLength": len(sequence_str),
                    "sequence": sequence_str,
                    "sequencePreview": sequence_str[:50] + "..." if len(sequence_str) > 50 else sequence_str,
                    "description": description,
                    "selected": selected,
                })

        return {
            "sequences": sequences,
            "sequenceCount": len(sequences),
        }

    except Exception as err:
        logger.exception("Error digesting CAsuDataFile %s", file_object, exc_info=err)
        return {
            "status": "Failed",
            "reason": f"Failed digesting CAsuDataFile: {err}",
            "digest": {},
        }


def digest_cdictdata_file_object(file_object: CPdbDataFile):
    if not isinstance(file_object, (CCP4ModelData.CDictDataFile, CDictDataFile)):
        return {
            "status": "Failed",
            "reason": "Not a CDictDataFile object",
            "digest": {},
        }
    file_path = file_object.fullPath.__str__()
    content_dict = parse_cif_ligand_summary(file_path)
    # Extract atoms/bonds for ALL monomers in the dictionary file
    monomers = extract_all_monomers_atoms_bonds(file_path)
    # Generate 2D molblocks for all monomers (best-effort, won't break digest on failure)
    try:
        molblocks = generate_all_molblocks(file_path)
    except Exception as e:
        logger.warning("Failed to generate molblocks for %s: %s", file_path, e)
        molblocks = {}
    return {
        "ligands": content_dict,
        "monomers": monomers,
        "molblocks": molblocks,
    }


def digest_cmtzdatafile_file_object(file_object):
    """
    Digest a CMtzDataFile by converting it to CGenericReflDataFile.

    CMtzDataFile inherits from CDataFile, not CGenericReflDataFile,
    so we create a CGenericReflDataFile with the same path to use
    the standard reflection data digestion.
    """
    if not file_object.isSet():
        return {"status": "Failed", "reason": "File object is not set", "digest": {}}
    try:
        # Create a CGenericReflDataFile with the same path
        generic_refl = CCP4XtalData.CGenericReflDataFile()
        generic_refl.setFullPath(str(file_object.fullPath))
        return digest_cgenericrefldatafile_file_object(generic_refl)
    except Exception as err:
        logger.exception("Error digesting CMtzDataFile %s", file_object, exc_info=err)
        return {
            "status": "Failed",
            "reason": f"Failed digesting CMtzDataFile {err}",
            "digest": {},
        }


def digest_cgenericrefldatafile_file_object(file_object: CGenericReflDataFile):
    content_dict = {}
    if not isinstance(file_object, CCP4XtalData.CGenericReflDataFile):
        return {
            "status": "Failed",
            "reason": "Not a CGenericReflDataFile object",
            "digest": {},
        }
    if not file_object.isSet():
        return {"status": "Failed", "reason": "File object is not set", "digest": {}}
    try:
        file_object.loadFile()
        file_object.setContentFlag()
        contents = file_object.getFileContent()
        content_dict = value_dict_for_object(contents)
        content_dict["format"] = file_object.getFormat()
        content_dict["merged"] = file_object.getMerged()

        # Content-based diagnosis (the single Python authority) as an ADDITIVE
        # `diagnosis` block: content-detected format, real merged/anomalous,
        # cell/SG/wavelength/resolution, StarAniso, and a `needs` list of
        # metadata the file lacks (SHELX cell/SG/dataType). The legacy top-level
        # `format` (extension-based) and `merged` (stub) keys are left untouched
        # for backward compatibility; the thin UI reads `diagnosis.*` instead.
        try:
            # Cached (path, mtime, size): the digest can be requested many times
            # for the same file across front-end renders, so the file read is
            # memoised rather than repeated on each call.
            from ccp4i2.lib.utils.files.reflection_diagnosis import (
                diagnose_reflection_file_cached,
            )
            content_dict["diagnosis"] = diagnose_reflection_file_cached(
                str(file_object.fullPath)
            )
        except Exception as diag_err:
            logger.warning(
                "diagnose_reflection_file failed for %s: %s",
                file_object.fullPath, diag_err,
            )

        # Initialize FreeR summary fields
        content_dict["hasFreeR"] = False
        content_dict["freerValid"] = False
        content_dict["freerWarnings"] = []

        if file_object.getFormat() == "mmcif":
            rblock_infos = _mmcif_block_infos(file_object.fullPath.__str__())
            content_dict["rblock_infos"] = rblock_infos

            # Set summary FreeR info from first block with merged data
            for block_info in rblock_infos:
                if block_info.get("hasFreeR"):
                    content_dict["hasFreeR"] = True
                    content_dict["freerValid"] = block_info.get("freerValid", False)
                    content_dict["freerWarnings"] = block_info.get("freerWarnings", [])
                    break

        # Add column groups for MTZ files
        if contents and hasattr(contents, 'getColumnGroups'):
            try:
                column_groups = contents.getColumnGroups()
                # Use value_dict_for_object for proper recursive CData serialization
                content_dict["columnGroups"] = [
                    value_dict_for_object(grp) for grp in column_groups
                ]

                # Extract FreeR info from MTZ column groups
                for grp in column_groups:
                    grp_type = str(grp.columnGroupType) if hasattr(grp, 'columnGroupType') else None
                    if grp_type == 'FreeR':
                        content_dict["hasFreeR"] = True
                        # For MTZ, FreeR is valid if the column exists
                        # More detailed validation would require reading the data
                        content_dict["freerValid"] = True
                        # Get the FreeR column label
                        if hasattr(grp, 'columnList') and len(grp.columnList) > 0:
                            freer_label = str(grp.columnList[0].columnLabel) if hasattr(grp.columnList[0], 'columnLabel') else None
                            content_dict["freerColumnLabel"] = freer_label
                        break

            except Exception as col_err:
                logger.warning("Failed to get column groups: %s", col_err)

        return content_dict
    except Exception as err:
        logger.exception("Error digesting file %s", file_object, exc_info=err)
        return {
            "status": "Failed",
            "reason": f"Failed digesting CGenericReflDataFile {err}",
            "digest": {},
        }


def digest_cdatafile_file_object(file_object: CDataFile):
    if not isinstance(file_object, CCP4File.CDataFile):
        return {"status": "Failed", "reason": "Not a CDataFile object", "digest": {}}
    if not file_object.isSet():
        return {"status": "Failed", "reason": "File object is not set", "digest": {}}
    try:
        data_type_name = identify_data_type(str(file_object.fullPath))
        if data_type_name in ["mtz", "sfcif"]:
            specific_object = CCP4XtalData.CGenericReflDataFile()
            specific_object.setFullPath(str(file_object.fullPath))
            return digest_cgenericrefldatafile_file_object(specific_object)
        elif data_type_name == "model":
            specific_object = CCP4ModelData.CPdbDataFile()
            specific_object.setFullPath(str(file_object.fullPath))
            return digest_cpdbdata_file_object(specific_object)
        elif data_type_name == "map":
            specific_object = CCP4XtalData.CMapDataFile()
            specific_object.setFullPath(str(file_object.fullPath))
            return digest_other_file_object(specific_object)
        elif data_type_name == "sequence":
            specific_object = CCP4ModelData.CSeqDataFile()
            specific_object.setFullPath(str(file_object.fullPath))
            return digest_cseqdata_file_object(specific_object)
        return digest_other_file_object(file_object)

    except Exception as err:
        logger.exception("Error digesting file %s", file_object, exc_info=err)
        return {
            "status": "Failed",
            "reason": f"Failed digesting CDataFile {err}",
            "digest": {},
        }
