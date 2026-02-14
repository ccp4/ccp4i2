"""
This file provides lazy loading
of plugin classes using explicit import statements.
"""

import importlib
import logging
import traceback

_PATHS = {
    "AMPLE": "ccp4i2.wrappers.AMPLE.script.AMPLE:AMPLE",
    "AUSPEX": "ccp4i2.wrappers.AUSPEX.script.auspex_wrapper:AUSPEX",
    "AcedrgLink": "ccp4i2.wrappers.AcedrgLink.script.AcedrgLink:AcedrgLink",
    "AlternativeImportXIA2": "ccp4i2.wrappers.AlternativeImportXIA2.script.AlternativeImportXIA2:AlternativeImportXIA2",
    "Lidia": "ccp4i2.wrappers.Lidia.script.Lidia:Lidia",
    "LidiaAcedrgNew": "ccp4i2.pipelines.LidiaAcedrgNew.script.LidiaAcedrgNew:LidiaAcedrgNew",
    "MakeLink": "ccp4i2.pipelines.MakeLink.script.MakeLink:MakeLink",
    "MakeMonster": "ccp4i2.wrappers.MakeMonster.script.MakeMonster:MakeMonster",
    "MakeProjectsAndDoLigandPipeline": "ccp4i2.pipelines.MakeProjectsAndDoLigandPipeline.script.MakeProjectsAndDoLigandPipeline:MakeProjectsAndDoLigandPipeline",
    "Platonyzer": "ccp4i2.wrappers.Platonyzer.script.Platonyzer:Platonyzer",
    "PrepareDeposit": "ccp4i2.pipelines.PrepareDeposit.script.PrepareDeposit:PrepareDeposit",
    "ProvideAlignment": "ccp4i2.wrappers.ProvideAlignment.script.ProvideAlignment:ProvideAlignment",
    "ProvideAsuContents": "ccp4i2.wrappers.ProvideAsuContents.script.ProvideAsuContents:ProvideAsuContents",
    "ProvideSequence": "ccp4i2.wrappers.ProvideSequence.script.ProvideSequence:ProvideSequence",
    "ProvideTLS": "ccp4i2.wrappers.ProvideTLS.script.ProvideTLS:ProvideTLS",
    "SIMBAD": "ccp4i2.wrappers.SIMBAD.script.SIMBAD:SIMBAD",
    "ShelxCD": "ccp4i2.wrappers.ShelxCDE.script.ShelxCD:ShelxCD",
    "SubstituteLigand": "ccp4i2.pipelines.SubstituteLigand.script.SubstituteLigand:SubstituteLigand",
    "SubtractNative": "ccp4i2.wrappers.SubtractNative.script.SubtractNative:SubtractNative",
    "TestObsConversions": "ccp4i2.wrappers.TestObsConversions.script.TestObsConversions:TestObsConversions",
    "acedrg": "ccp4i2.wrappers.acedrg.script.acedrg:acedrg",
    "acedrgNew": "ccp4i2.wrappers.acedrgNew.script.acedrgNew:acedrgNew",
    "acorn": "ccp4i2.wrappers.acorn.script.acorn:acorn",
    "add_fractional_coords": "ccp4i2.wrappers.add_fractional_coords.script.add_fractional_coords:add_fractional_coords",
    "adding_stats_to_mmcif_i2": "ccp4i2.wrappers.adding_stats_to_mmcif_i2.script.adding_stats_to_mmcif_i2:adding_stats_to_mmcif_i2",
    "aimless": "ccp4i2.wrappers.aimless.script.aimless:aimless",
    "aimless_pipe": "ccp4i2.pipelines.aimless_pipe.script.aimless_pipe:aimless_pipe",
    "arcimboldo": "ccp4i2.wrappers.arcimboldo.script.arcimboldo:arcimboldo",
    "arp_warp_classic": "ccp4i2.wrappers.arp_warp_classic.script.arp_warp_classic:arp_warp_classic",
    "baverage": "ccp4i2.wrappers.baverage.script.baverage:baverage",
    "buster": "ccp4i2.wrappers.buster.script.buster:buster",
    "cad_copy_column": "ccp4i2.wrappers.cad_copy_column.script.cad_copy_column:cad_copy_column",
    "ccp4mg_edit_model": "ccp4i2.wrappers.ccp4mg_edit_model.script.ccp4mg_edit_model:ccp4mg_edit_model",
    "ccp4mg_edit_nomrbump": "ccp4i2.wrappers.ccp4mg_edit_nomrbump.script.ccp4mg_edit_nomrbump:ccp4mg_edit_nomrbump",
    "ccp4mg_general": "ccp4i2.wrappers.ccp4mg_general.script.ccp4mg_general:ccp4mg_general",
    "chainsaw": "ccp4i2.wrappers.chainsaw.script.chainsaw:chainsaw",
    "chltofom": "ccp4i2.wrappers.chltofom.script.chltofom:chltofom",
    "cif2mtz": "ccp4i2.wrappers.cif2mtz.script.cif2mtz:cif2mtz",
    "clustalw": "ccp4i2.wrappers.clustalw.script.clustalw:clustalw",
    "cmapcoeff": "ccp4i2.wrappers.cmapcoeff.script.cmapcoeff:cmapcoeff",
    "comit": "ccp4i2.wrappers.comit.script.comit:comit",
    "convert2mtz": "ccp4i2.wrappers.convert2mtz.script.convert2mtz:convert2mtz",
    "coordinate_selector": "ccp4i2.wrappers.coordinate_selector.script.coordinate_selector:coordinate_selector",
    "coot1": "ccp4i2.wrappers.coot1.script.coot1:coot1",
    "coot_find_ligand": "ccp4i2.wrappers.coot_find_ligand.script.coot_find_ligand:coot_find_ligand",
    "coot_find_waters": "ccp4i2.wrappers.coot_find_waters.script.coot_find_waters:coot_find_waters",
    "coot_rebuild": "ccp4i2.wrappers.coot_rebuild.script.coot_rebuild:coot_rebuild",
    "coot_rsr_morph": "ccp4i2.wrappers.coot_rsr_morph.script.coot_rsr_morph:coot_rsr_morph",
    "coot_script_lines": "ccp4i2.wrappers.coot_script_lines.script.coot_script_lines:coot_script_lines",
    "cpatterson": "ccp4i2.wrappers.cpatterson.script.cpatterson:cpatterson",
    "cphasematch": "ccp4i2.wrappers.cphasematch.script.cphasematch:cphasematch",
    "crank2": "ccp4i2.pipelines.crank2.script.crank2_script:crank2",
    "crank2_comb_phdmmb": "ccp4i2.pipelines.crank2.wrappers.crank2_comb_phdmmb.script.crank2_comb_phdmmb:crank2_comb_phdmmb",
    "crank2_createfree": "ccp4i2.pipelines.crank2.wrappers.crank2_createfree.script.crank2_createfree:crank2_createfree",
    "crank2_dmfull": "ccp4i2.pipelines.crank2.wrappers.crank2_dmfull.script.crank2_dmfull:crank2_dmfull",
    "crank2_faest": "ccp4i2.pipelines.crank2.wrappers.crank2_faest.script.crank2_faest:crank2_faest",
    "crank2_handdet": "ccp4i2.pipelines.crank2.wrappers.crank2_handdet.script.crank2_handdet:crank2_handdet",
    "crank2_mbref": "ccp4i2.pipelines.crank2.wrappers.crank2_mbref.script.crank2_mbref:crank2_mbref",
    "crank2_phas": "ccp4i2.pipelines.crank2.wrappers.crank2_phas.script.crank2_phas:crank2_phas",
    "crank2_phdmmb": "ccp4i2.pipelines.crank2.wrappers.crank2_phdmmb.script.crank2_phdmmb:crank2_phdmmb",
    "crank2_ref": "ccp4i2.pipelines.crank2.wrappers.crank2_ref.script.crank2_ref:crank2_ref",
    "crank2_refatompick": "ccp4i2.pipelines.crank2.wrappers.crank2_refatompick.script.crank2_refatompick:crank2_refatompick",
    "crank2_substrdet": "ccp4i2.pipelines.crank2.wrappers.crank2_substrdet.script.crank2_substrdet:crank2_substrdet",
    "csymmatch": "ccp4i2.wrappers.csymmatch.script.csymmatch:csymmatch",
    "ctruncate": "ccp4i2.wrappers.ctruncate.script.ctruncate:ctruncate",
    "density_calculator": "ccp4i2.wrappers.density_calculator.script.density_calculator:density_calculator",
    "dials_image": "ccp4i2.wrappers.dials_image.script.dials_image:dials_image",
    "dials_rlattice": "ccp4i2.wrappers.dials_rlattice.script.dials_rlattice:dials_rlattice",
    "dr_mr_modelbuild_pipeline": "ccp4i2.pipelines.dr_mr_modelbuild_pipeline.script.dr_mr_modelbuild_pipeline:dr_mr_modelbuild_pipeline",
    "dui": "ccp4i2.wrappers.dui.script.dui:dui",
    "editbfac": "ccp4i2.wrappers.editbfac.script.editbfac:editbfac",
    "edstats": "ccp4i2.wrappers.edstats.script.edstats:edstats",
    "fft": "ccp4i2.wrappers.fft.script.fft:fft",
    "findmyseq": "ccp4i2.wrappers.findmyseq.script.findmyseq:findmyseq",
    "freerflag": "ccp4i2.wrappers.freerflag.script.freerflag:freerflag",
    "gesamt": "ccp4i2.wrappers.gesamt.script.gesamt:gesamt",
    "hklin2cif": "ccp4i2.pipelines.PrepareDeposit.wrappers.hklin2cif.script.hklin2cif:hklin2cif",
    "i2Dimple": "ccp4i2.wrappers.i2Dimple.script.i2Dimple:i2Dimple",
    "imosflm": "ccp4i2.wrappers.imosflm.script.imosflm:imosflm",
    "import_merged": "ccp4i2.pipelines.import_merged.script.import_merged:import_merged",
    "import_mosflm": "ccp4i2.wrappers.import_mosflm.script.import_mosflm:import_mosflm",
    "import_serial": "ccp4i2.wrappers.import_serial.script.import_serial:import_serial",
    "import_serial_pipe": "ccp4i2.pipelines.import_serial_pipe.script.import_serial_pipe:import_serial_pipe",
    "import_xia2": "ccp4i2.pipelines.import_xia2.script.import_xia2:import_xia2",
    "lorestr_i2": "ccp4i2.wrappers.lorestr_i2.script.lorestr_i2:lorestr_i2",
    "mergeMtz": "ccp4i2.wrappers.mergeMtz.script.mergeMtz:mergeMtz",
    "metalCoord": "ccp4i2.wrappers.metalCoord.script.metalCoord:metalCoord",
    "modelASUCheck": "ccp4i2.wrappers.modelASUCheck.script.modelASUCheck:modelASUCheck",
    "modelcraft": "ccp4i2.wrappers.modelcraft.script.modelcraft:modelcraft",
    "molrep_den": "ccp4i2.wrappers.molrep_den.script.molrep_den:molrep_den",
    "molrep_mr": "ccp4i2.wrappers.molrep_mr.script.molrep_mr:molrep_mr",
    "molrep_pipe": "ccp4i2.pipelines.molrep_pipe.script.molrep_pipe:molrep_pipe",
    "molrep_selfrot": "ccp4i2.wrappers.molrep_selfrot.script.molrep_selfrot:molrep_selfrot",
    "morda_i2": "ccp4i2.wrappers.morda_i2.script.morda_i2:morda_i2",
    "mosflm": "ccp4i2.wrappers.mosflm.script.mosflm:mosflm",
    "mrbump_basic": "ccp4i2.wrappers.mrbump_basic.script.mrbump_basic:mrbump_basic",
    "mrbump_model_prep": "ccp4i2.pipelines.dr_mr_modelbuild_pipeline.wrappers.mrbump_model_prep.script.mrbump_model_prep:mrbump_model_prep",
    "mrparse": "ccp4i2.wrappers.mrparse.script.mrparse_wrapper:mrparse",
    "mrparse_simple": "ccp4i2.pipelines.dr_mr_modelbuild_pipeline.wrappers.mrparse_simple.script.mrparse_simple_wrapper:mrparse_simple",
    "mtzheader": "ccp4i2.wrappers.mtzheader.script.mtzheader:mtzheader",
    "mtzutils": "ccp4i2.wrappers.mtzutils.script.mtzutils:mtzutils",
    "pairef": "ccp4i2.wrappers.pairef.script.pairef:pairef",
    "parrot": "ccp4i2.wrappers.parrot.script.parrot:parrot",
    "pdb_extract_wrapper": "ccp4i2.pipelines.PrepareDeposit.wrappers.pdb_extract_wrapper.script.pdb_extract_wrapper:pdb_extract_wrapper",
    "pdb_redo_api": "ccp4i2.wrappers.pdb_redo_api.script.pdb_redo_api:pdb_redo_api",
    "pdbset_ui": "ccp4i2.wrappers.pdbset_ui.script.pdbset_ui:pdbset_ui",
    "pdbview_edit": "ccp4i2.wrappers.pdbview_edit.script.pdbview_edit:pdbview_edit",
    "phaser_EP": "ccp4i2.pipelines.phaser_ep.script.phaser_EP:phaser_EP",
    "phaser_EP_AUTO": "ccp4i2.pipelines.phaser_pipeline.wrappers.phaser_EP_AUTO.script.phaser_EP_AUTO:phaser_EP_AUTO",
    "phaser_EP_LLG": "ccp4i2.pipelines.phaser_pipeline.wrappers.phaser_EP_LLG.script.phaser_EP_LLG:phaser_EP_LLG",
    "phaser_MR": "ccp4i2.pipelines.phaser_pipeline.wrappers.phaser_MR.script.phaser_MR:phaser_MR",
    "phaser_MR_AUTO": "ccp4i2.pipelines.phaser_pipeline.wrappers.phaser_MR_AUTO.script.phaser_MR_AUTO:phaser_MR_AUTO",
    "phaser_MR_FRF": "ccp4i2.pipelines.phaser_pipeline.wrappers.phaser_MR_FRF.script.phaser_MR_FRF:phaser_MR_FRF",
    "phaser_MR_FTF": "ccp4i2.pipelines.phaser_pipeline.wrappers.phaser_MR_FTF.script.phaser_MR_FTF:phaser_MR_FTF",
    "phaser_MR_PAK": "ccp4i2.pipelines.phaser_pipeline.wrappers.phaser_MR_PAK.script.phaser_MR_PAK:phaser_MR_PAK",
    "phaser_MR_RNP": "ccp4i2.pipelines.phaser_pipeline.wrappers.phaser_MR_RNP.script.phaser_MR_RNP:phaser_MR_RNP",
    "phaser_analysis": "ccp4i2.wrappers.phaser_analysis.script.phaser_analysis:phaser_analysis",
    "phaser_ensembler": "ccp4i2.wrappers.phaser_ensembler.script.phaser_ensembler:phaser_ensembler",
    "phaser_phil": "ccp4i2.wrappers.phaser_phil.script.phaser_phil:phaser_phil",
    "phaser_pipeline": "ccp4i2.pipelines.phaser_pipeline.script.phaser_pipeline:phaser_pipeline",
    "phaser_rnp_pipeline": "ccp4i2.pipelines.phaser_rnp_pipeline.script.phaser_rnp_pipeline:phaser_rnp_pipeline",
    "phaser_simple": "ccp4i2.pipelines.phaser_simple.script.phaser_simple:phaser_simple",
    "phaser_singleMR": "ccp4i2.wrappers.phaser_singleMR.script.phaser_singleMR:phaser_singleMR",
    "pisa_analyse": "ccp4i2.pipelines.pisapipe.wrappers.pisa_analyse.script.pisa_analyse:pisa_analyse",
    "pisa_list": "ccp4i2.pipelines.pisapipe.wrappers.pisa_list.script.pisa_list:pisa_list",
    "pisa_xml": "ccp4i2.pipelines.pisapipe.wrappers.pisa_xml.script.pisa_xml:pisa_xml",
    "pisapipe": "ccp4i2.pipelines.pisapipe.script.pisapipe:pisapipe",
    "pointless": "ccp4i2.wrappers.pointless.script.pointless:pointless",
    "pointless_reindexToMatch": "ccp4i2.wrappers.pointless_reindexToMatch.script.pointless_reindexToMatch:pointless_reindexToMatch",
    "privateer": "ccp4i2.wrappers.privateer.script.privateer_wrapper:privateer",
    "prosmart": "ccp4i2.wrappers.prosmart.script.prosmart:prosmart",
    "prosmart_refmac": "ccp4i2.pipelines.prosmart_refmac.script.prosmart_refmac:prosmart_refmac",
    "pyphaser_mr": "ccp4i2.wrappers.pyphaser_mr.script.pyphaser_mr:pyphaser_mr",
    "qtpisa": "ccp4i2.wrappers.qtpisa.script.qtpisa:qtpisa",
    "refmac": "ccp4i2.wrappers.refmac.script.refmac:refmac",
    "scaleit": "ccp4i2.wrappers.scaleit.script.scaleit:scaleit",
    "scalepack2mtz": "ccp4i2.wrappers.scalepack2mtz.script.scalepack2mtz:scalepack2mtz",
    "sculptor": "ccp4i2.wrappers.sculptor.script.sculptor:sculptor",
    "servalcat": "ccp4i2.wrappers.servalcat.script.servalcat:servalcat",
    "servalcat_pipe": "ccp4i2.pipelines.servalcat_pipe.script.servalcat_pipe:servalcat_pipe",
    "sheetbend": "ccp4i2.wrappers.sheetbend.script.sheetbend:sheetbend",
    "shelx": "ccp4i2.pipelines.shelx.script.shelx_script:shelx",
    "shelxeMR": "ccp4i2.wrappers.shelxeMR.script.shelxeMR:shelxeMR",
    "slicendice": "ccp4i2.wrappers.slicendice.script.slicendice:slicendice",
    "splitMtz": "ccp4i2.wrappers.splitMtz.script.splitMtz:splitMtz",
    "tableone": "ccp4i2.pipelines.tableone.script.tableone:tableone",
    "validate_protein": "ccp4i2.wrappers.validate_protein.script.validate_protein:validate_protein",
    "x2mtz": "ccp4i2.wrappers.x2mtz.script.x2mtz:x2mtz",
    "xia2_aimless": "ccp4i2.pipelines.import_xia2.wrappers.xia2_aimless.script.xia2_aimless:xia2_aimless",
    "xia2_ctruncate": "ccp4i2.pipelines.import_xia2.wrappers.xia2_ctruncate.script.xia2_ctruncate:xia2_ctruncate",
    "xia2_dials": "ccp4i2.wrappers.xia2_dials.script.xia2_dials:xia2_dials",
    "xia2_integration": "ccp4i2.pipelines.import_xia2.wrappers.xia2_integration.script.xia2_integration:xia2_integration",
    "xia2_multiplex": "ccp4i2.wrappers.xia2_multiplex.script.xia2_multiplex:xia2_multiplex",
    "xia2_pointless": "ccp4i2.pipelines.import_xia2.wrappers.xia2_pointless.script.xia2_pointless:xia2_pointless",
    "xia2_run": "ccp4i2.pipelines.import_xia2.wrappers.xia2_run.script.xia2_run:xia2_run",
    "xia2_ssx_reduce": "ccp4i2.wrappers.xia2_ssx_reduce.script.xia2_ssx_reduce:xia2_ssx_reduce",
    "xia2_xds": "ccp4i2.wrappers.xia2_xds.script.xia2_xds:xia2_xds",
    "zanuda": "ccp4i2.wrappers.zanuda.script.zanuda:zanuda",
}

_CACHE = {}


def get_plugin_class(task_name: str):
    """
    Get a plugin class by name.

    The plugin is imported lazily on first access and cached.

    Args:
        task_name: Name of the task/plugin

    Returns:
        Plugin class, or None if not found
    """
    if task_name in _CACHE:
        return _CACHE[task_name]

    path = _PATHS.get(task_name)
    if path is None:
        return None

    try:
        module_name, class_name = path.split(":")
        module = importlib.import_module(module_name)
        _CACHE[task_name] = getattr(module, class_name)
        return _CACHE[task_name]
    except Exception as e:
        logger = logging.getLogger(f"ccp4i2:{__name__}")
        logger.error(f"Failed to import plugin {task_name}: {e}")
        logger.error(f"Traceback:\n{traceback.format_exc()}")
        return None
