import json
import logging
import gemmi
import importlib
from typing import Dict, Type

from core import CCP4File
from core import CCP4XtalData
from core import CCP4ModelData
from core.CCP4Container import CContainer
from core.base_object.cdata_file import CDataFile
from core.base_object.cdata import CData
from core.CCP4XtalData import CGenericReflDataFile, CMapDataFile, CMtzDataFile
from core.CCP4ModelData import CPdbDataFile, CDictDataFile
# Import stub class for isinstance checks - subclasses like CObsDataFile inherit from
# stubs (CMtzDataFileStub) not implementations (CMtzDataFile)
from core.cdata_stubs.CCP4XtalData import CMtzDataFileStub
from pipelines.import_merged.script import mmcifutils
from ..containers.find_objects import find_objects
from ..containers.get_container import get_job_container
from ..containers.json_encoder import CCP4i2JsonEncoder
from ..plugins.plugin_context import get_plugin_with_context
from ..formats.cif_ligand import parse_cif_ligand_summary
from ..parameters.value_dict import value_dict_for_object
from ....db import models
from ...parse import identify_data_type

logger = logging.getLogger(f"ccp4x:{__name__}")


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


def _build_class_registry() -> Dict[str, Type[CData]]:
    """
    Build registry of available CData classes by name.

    This follows the same pattern as core/task_manager/def_xml_handler.py
    to dynamically discover all CData implementation classes.
    """
    registry = {}

    # Import from implementation modules
    implementation_modules = [
        'CCP4File',
        'CCP4ModelData',
        'CCP4XtalData',
    ]

    for module_name in implementation_modules:
        try:
            module = importlib.import_module(f'core.{module_name}')
            for attr_name in dir(module):
                if attr_name.startswith('_'):
                    continue
                attr = getattr(module, attr_name)
                if (
                    isinstance(attr, type)
                    and issubclass(attr, CData)
                    and attr is not CData
                ):
                    registry[attr.__name__] = attr
        except ImportError as e:
            logger.warning(f"Could not import {module_name}: {e}")
            continue

    return registry


# Build class registry at module level for use in digest_file()
CLASS_REGISTRY = _build_class_registry()

FILETYPES_TEXT = [
    "Unknown",
    "application/CCP4-seq",
    "chemical/x-pdb",
    "MultiPDB",
    "application/CCP4-mtz",
    "application/CCP4-unmerged-mtz",
    "application/CCP4-unmerged-experimental",
    "application/CCP4-map",
    "application/refmac-dictionary",
    "application/refmac-TLS",
    "application/CCP4-mtz-freerflag",
    "application/CCP4-mtz-observed",
    "application/CCP4-mtz-phases",
    "application/CCP4-mtz-map",
    "Dummy",
    "application/CCP4-seqalign",
    "application/CCP4-mtz-mini",
    "application/coot-script",
    "application/refmac-external-restraints",
    "application/CCP4-scene",
    "application/CCP4-shelx-FA",
    "application/phaser-sol",
    "chemical/x-mdl-molfile",
    "application/iMosflm-xml",
    "application/CCP4-image",
    "application/CCP4-generic-reflections",
    "application/HHPred-alignments",
    "application/Blast-alignments",
    "chemical/x-pdb-ensemble",
    "application/CCP4-asu-content",
    "application/dials-jfile",
    "application/dials-pfile",
    "application/phaser-rfile",
    "application/refmac-keywords",
    "application/x-pdf",
    "application/postscript",
    "application/EBI-validation-xml",
    "chemical/x-cif",
]
FILETYPES_CLASS = [
    "DataFile",
    "SeqDataFile",
    "PdbDataFile",
    "",
    "MtzDataFile",
    "MtzDataFile",
    "UnmergedDataFile",
    "MapDataFile",
    "DictDataFile",
    "TLSDataFile",
    "FreeRDataFile",
    "ObsDataFile",
    "PhsDataFile",
    "MapCoeffsDataFile",
    "",
    "SeqAlignDataFile",
    "MiniMtzDataFile",
    "CootHistoryDataFile",
    "RefmacRestraintsDataFile",
    "SceneDataFile",
    "ShelxFADataFile",
    "PhaserSolDataFile",
    "MDLMolDataFile",
    "ImosflmXmlDataFile",
    "ImageFile",
    "GenericReflDataFile",
    "HhpredDataFile",
    "BlastDataFile",
    "EnsemblePdbDataFile",
    "AsuDataFile",
    "DialsJsonFile",
    "DialsPickleFile",
    "PhaserRFileDataFile",
    "RefmacKeywordFile",
    "PDFDataFile",
    "PostscriptDataFile",
    "EBIValidationXMLDataFile",
    "MmcifReflDataFile",
]


def is_basic_type(obj):
    return isinstance(obj, (str, int, float, bool, type(None)))


def is_custom_class_instance(obj):
    return not is_basic_type(obj) and not isinstance(obj, (list, dict, tuple))


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
    the_class = CLASS_REGISTRY.get(full_class_name)

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


def digest_file_object(file_object: CDataFile):
    if not isinstance(file_object, CCP4File.CDataFile):
        return {"status": "Failed", "reason": "Not a valid file object", "digest": {}}
    if not file_object.isSet():
        return {"status": "Failed", "reason": "File object is not set", "digest": {}}
    if isinstance(file_object, CCP4ModelData.CPdbDataFile):
        return digest_cpdbdata_file_object(file_object)
    if isinstance(file_object, CCP4XtalData.CGenericReflDataFile):
        return digest_cgenericrefldatafile_file_object(file_object)
    # CMtzDataFile inherits from CDataFile, not CGenericReflDataFile, so check separately
    # Use CMtzDataFileStub for isinstance check because subclasses inherit from stubs
    if isinstance(file_object, CMtzDataFileStub):
        return digest_cmtzdatafile_file_object(file_object)
    if isinstance(file_object, CCP4ModelData.CSeqDataFile):
        return digest_cseqdata_file_object(file_object)
    if isinstance(file_object, (CCP4ModelData.CDictDataFile, CDictDataFile)):
        return digest_cdictdata_file_object(file_object)
    if type(file_object) is CCP4File.CDataFile:
        return digest_cdatafile_file_object(file_object)
    return digest_other_file_object(file_object)


def digest_other_file_object(file_object: CDataFile):
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
            "reason": f"Failed digesting CDataFile {err}",
            "digest": {},
        }


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
        return content_dict
    except Exception as err:
        logger.exception("Error digesting file %s", file_object, exc_info=err)
        return {
            "status": "Failed",
            "reason": f"Failed digesting CPdbDataFile {err}",
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


def digest_cdictdata_file_object(file_object: CPdbDataFile):
    if not isinstance(file_object, (CCP4ModelData.CDictDataFile, CDictDataFile)):
        return {
            "status": "Failed",
            "reason": "Not a CDictDataFile object",
            "digest": {},
        }
    content_dict = parse_cif_ligand_summary(file_object.fullPath.__str__())
    return content_dict


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
        if file_object.getFormat() == "mmcif":
            mmcif = gemmi.cif.read_file(file_object.fullPath.__str__())
            rblocks = gemmi.as_refln_blocks(mmcif)
            rblock_infos = []
            for rb in rblocks:
                blkinfo = mmcifutils.CifBlockInfo(rb)
                rblock_infos.append(flatten_instance(blkinfo))
            content_dict["rblock_infos"] = rblock_infos
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
        result = identify_data_type(str(file_object.fullPath))
        if result["data_type_name"] in ["mtz", "sfcif"]:
            specific_object = CCP4XtalData.CGenericReflDataFile()
            specific_object.setFullPath(str(file_object.fullPath))
            return digest_cgenericrefldatafile_file_object(specific_object)
        elif result["data_type_name"] == "model":
            specific_object = CCP4ModelData.CPdbDataFile()
            specific_object.setFullPath(str(file_object.fullPath))
            return digest_cpdbdata_file_object(specific_object)
        elif result["data_type_name"] == "map":
            specific_object = CCP4XtalData.CMapDataFile()
            specific_object.setFullPath(str(file_object.fullPath))
            return digest_other_file_object(specific_object)
        elif result["data_type_name"] == "sequence":
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
