"""
This file provides lazy loading
of report classes using explicit import statements.
"""

import importlib
import warnings

_PATHS = {
    "AMPLE": "ccp4i2.wrappers.AMPLE.script.AMPLE_report:AMPLE_report",
    "AUSPEX": "ccp4i2.wrappers.AUSPEX.script.AUSPEX_report:AUSPEX_report",
    "Acedrg": "ccp4i2.pipelines.LidiaAcedrgNew.script.LidiaAcedrgNew_report:acedrgNew_report",
    "AcedrgLink": "ccp4i2.wrappers.AcedrgLink.script.AcedrgLink_report:AcedrgLink_report",
    "AlternativeImportXIA2": "ccp4i2.wrappers.AlternativeImportXIA2.script.AlternativeImportXIA2_report:AlternativeImportXIA2_report",
    "LidiaAcedrgNew": "ccp4i2.pipelines.LidiaAcedrgNew.script.LidiaAcedrgNew_report:LidiaAcedrgNew_report",
    "MakeLink": "ccp4i2.pipelines.MakeLink.script.MakeLink_report:MakeLink_report",
    "MakeMonster": "ccp4i2.wrappers.MakeMonster.script.MakeMonster_report:MakeMonster_report",
    "MakeProjectsAndDoLigandPipeline": "ccp4i2.pipelines.MakeProjectsAndDoLigandPipeline.script.MakeProjectsAndDoLigandPipeline_report:MakeProjectsAndDoLigandPipeline_report",
    "PrepareDeposit": "ccp4i2.pipelines.PrepareDeposit.script.PrepareDeposit_report:PrepareDeposit_report",
    "ProvideAlignment": "ccp4i2.wrappers.ProvideAlignment.script.ProvideAlignment_report:ProvideAlignment_report",
    "ProvideAsuContents": "ccp4i2.wrappers.ProvideAsuContents.script.ProvideAsuContents_report:ProvideAsuContents_report",
    "ProvideSequence": "ccp4i2.wrappers.ProvideSequence.script.ProvideSequence_report:ProvideSequence_report",
    "ProvideTLS": "ccp4i2.wrappers.ProvideTLS.script.ProvideTLS_report:ProvideTLS_report",
    "RvapiReport": "ccp4i2.wrappers.morda_i2.script.morda_i2_report:RvapiReport",
    "SIMBAD": "ccp4i2.wrappers.SIMBAD.script.SIMBAD_report:SIMBAD_report",
    "ShelxCD": "ccp4i2.wrappers.ShelxCDE.script.ShelxCD_report:ShelxCD_report",
    "SubstituteLigand": "ccp4i2.pipelines.SubstituteLigand.script.SubstituteLigand_report:SubstituteLigand_report",
    "SubtractNative": "ccp4i2.wrappers.SubtractNative.script.SubtractNative_report:SubtractNative_report",
    "TestObsConversions": "ccp4i2.wrappers.TestObsConversions.script.TestObsConversions_report:TestObsConversions_report",
    "acorn": "ccp4i2.wrappers.acorn.script.acorn_report:acorn_report",
    "add_fractional_coords": "ccp4i2.wrappers.add_fractional_coords.script.add_fractional_coords_report:add_fractional_coords_report",
    "adding_stats_to_mmcif_i2": "ccp4i2.wrappers.adding_stats_to_mmcif_i2.script.adding_stats_to_mmcif_i2_report:adding_stats_to_mmcif_i2_report",
    "aimless": "ccp4i2.wrappers.adding_stats_to_mmcif_i2.script.adding_stats_to_mmcif_i2_report:aimless_report",
    "arcimboldo": "ccp4i2.wrappers.arcimboldo.script.arcimboldo_report:arcimboldo_report",
    "arp_warp_classic": "ccp4i2.wrappers.arp_warp_classic.script.arp_warp_classic_report:arp_warp_classic_report",
    "buster": "ccp4i2.wrappers.buster.script.buster_report:buster_report",
    "ccp4mg_edit_model": "ccp4i2.wrappers.ccp4mg_edit_model.script.ccp4mg_edit_model_report:ccp4mg_edit_model_report",
    "ccp4mg_edit_nomrbump": "ccp4i2.wrappers.ccp4mg_edit_nomrbump.script.ccp4mg_edit_nomrbump_report:ccp4mg_edit_no_mrbump_report",
    "ccp4mg_general": "ccp4i2.wrappers.ccp4mg_general.script.ccp4mg_general_report:ccp4mg_general_report",
    "chainsaw": "ccp4i2.wrappers.chainsaw.script.chainsaw_report:chainsaw_report",
    "chltofom": "ccp4i2.wrappers.chltofom.script.chltofom_report:chltofom_report",
    "cif2mtz": "ccp4i2.wrappers.cif2mtz.script.cif2mtz_report:cif2mtz_report",
    "clustalw": "ccp4i2.wrappers.clustalw.script.clustalw_report:clustalw_report",
    "cmapcoeff": "ccp4i2.wrappers.cmapcoeff.script.cmapcoeff_report:cmapcoeff_report",
    "comit": "ccp4i2.wrappers.comit.script.comit_report:comit_report",
    "coordinate_selector": "ccp4i2.wrappers.coordinate_selector.script.coordinate_selector_report:coordinate_selector_report",
    "coot1": "ccp4i2.wrappers.coot1.script.coot1_report:coot1_report",
    "coot_find_ligand": "ccp4i2.wrappers.coot_find_ligand.script.coot_find_ligand_report:coot_find_ligand_report",
    "coot_find_waters": "ccp4i2.wrappers.coot_find_waters.script.coot_find_waters_report:coot_find_waters_report",
    "coot_rebuild": "ccp4i2.wrappers.coot_rebuild.script.coot_rebuild_report:coot_rebuild_report",
    "coot_rsr_morph": "ccp4i2.wrappers.coot_rsr_morph.script.coot_rsr_morph_report:coot_rsr_morph_report",
    "coot_script_lines": "ccp4i2.wrappers.coot_script_lines.script.coot_script_lines_report:coot_script_lines_report",
    "cpatterson": "ccp4i2.wrappers.cpatterson.script.cpatterson_report:cpatterson_report",
    "cphasematch": "ccp4i2.wrappers.cphasematch.script.cphasematch_report:cphasematch_report",
    "crank2": "ccp4i2.pipelines.crank2.script.crank2_report:crank2_report",
    "crank2_comb_phdmmb": "ccp4i2.pipelines.crank2.wrappers.crank2_comb_phdmmb.script.crank2_comb_phdmmb_report:crank2_comb_phdmmb_report",
    "crank2_dmfull": "ccp4i2.pipelines.crank2.wrappers.crank2_dmfull.script.crank2_dmfull_report:crank2_dmfull_report",
    "crank2_faest": "ccp4i2.pipelines.crank2.wrappers.crank2_faest.script.crank2_faest_report:crank2_faest_report",
    "crank2_handdet": "ccp4i2.pipelines.crank2.wrappers.crank2_handdet.script.crank2_handdet_report:crank2_handdet_report",
    "crank2_mbref": "ccp4i2.pipelines.crank2.wrappers.crank2_mbref.script.crank2_mbref_report:crank2_mbref_report",
    "crank2_phas": "ccp4i2.pipelines.crank2.wrappers.crank2_phas.script.crank2_phas_report:crank2_phas_report",
    "crank2_phdmmb": "ccp4i2.pipelines.crank2.wrappers.crank2_phdmmb.script.crank2_phdmmb_report:crank2_phdmmb_report",
    "crank2_ref": "ccp4i2.pipelines.crank2.wrappers.crank2_ref.script.crank2_ref_report:crank2_ref_report",
    "crank2_refatompick": "ccp4i2.pipelines.crank2.wrappers.crank2_refatompick.script.crank2_refatompick_report:crank2_refatompick_report",
    "crank2_substrdet": "ccp4i2.pipelines.crank2.wrappers.crank2_substrdet.script.crank2_substrdet_report:crank2_substrdet_report",
    "csymmatch": "ccp4i2.pipelines.phaser_pipeline.script.phaser_pipeline_report:csymmatch_report",
    "ctruncate": "ccp4i2.wrappers.ctruncate.script.ctruncate_report:ctruncate_report",
    "density_calculator": "ccp4i2.wrappers.density_calculator.script.density_calculator_report:density_calculator_report",
    "dials_image": "ccp4i2.wrappers.dials_image.script.dials_image_report:dials_image_report",
    "dials_rlattice": "ccp4i2.wrappers.dials_rlattice.script.dials_rlattice_report:dials_rlattice_report",
    "dr_mr_modelbuild_pipeline": "ccp4i2.pipelines.dr_mr_modelbuild_pipeline.script.dr_mr_modelbuild_pipeline_report:dr_mr_modelbuild_pipeline_report",
    "dui": "ccp4i2.wrappers.dui.script.dui_report:dui_report",
    "editbfac": "ccp4i2.wrappers.editbfac.script.editbfac_report:editbfac_report",
    "edstats": "ccp4i2.wrappers.edstats.script.edstats_report:edstats_report",
    "fft": "ccp4i2.wrappers.fft.script.fft_report:fft_report",
    "findmyseq": "ccp4i2.wrappers.findmyseq.script.findmyseq_report:findmyseq_report",
    "freerflag": "ccp4i2.wrappers.freerflag.script.freerflag_report:freerflag_report",
    "gesamt": "ccp4i2.wrappers.gesamt.script.gesamt_report:gesamt_report",
    "i2Dimple": "ccp4i2.wrappers.i2Dimple.script.i2Dimple_report:i2Dimple_report",
    "imosflm": "ccp4i2.wrappers.imosflm.script.imosflm_report:imosflm_report",
    "import_files": "ccp4i2.wrappers2.import_files.script.import_files_report:import_files_report",
    "import_merged": "ccp4i2.pipelines.import_merged.script.import_merged_report:import_merged_report",
    "import_mosflm": "ccp4i2.wrappers.import_mosflm.script.import_mosflm_report:import_mosflm_report",
    "import_serial": "ccp4i2.wrappers.import_serial.script.import_serial_report:import_serial_report",
    "import_serial_pipe": "ccp4i2.pipelines.import_serial_pipe.script.import_serial_pipe_report:import_serial_pipe_report",
    "import_xia2": "ccp4i2.pipelines.import_xia2.script.import_xia2_report:import_xia2_report",
    "lorestr_i2": "ccp4i2.wrappers.lorestr_i2.script.lorestr_i2_report:lorestr_i2_report",
    "mergeMtz": "ccp4i2.wrappers.mergeMtz.script.mergeMtz_report:mergeMtz_report",
    "metalCoord": "ccp4i2.wrappers.metalCoord.script.metalCoord_report:metalCoord_report",
    "modelASUCheck": "ccp4i2.wrappers.modelASUCheck.script.modelASUCheck_report:modelASUCheck_report",
    "modelcraft": "ccp4i2.pipelines.phaser_ep.script.phaser_EP_report:modelcraft_report",
    "molrep_den": "ccp4i2.wrappers.molrep_den.script.molrep_den_report:molrep_den_report",
    "molrep_mr": "ccp4i2.wrappers.molrep_mr.script.molrep_mr_report:molrep_mr_report",
    "molrep_pipe": "ccp4i2.pipelines.molrep_pipe.script.molrep_pipe_report:molrep_pipe_report",
    "molrep_selfrot": "ccp4i2.wrappers.molrep_selfrot.script.molrep_selfrot_report:molrep_selfrot_report",
    "morda_i2": "ccp4i2.wrappers.morda_i2.script.morda_i2_report:morda_i2_report",
    "mosflm": "ccp4i2.wrappers.mosflm.script.mosflm_report:mosflm_report",
    "mrbump_basic": "ccp4i2.wrappers.mrbump_basic.script.mrbump_basic_report:mrbump_basic_report",
    "mrparse": "ccp4i2.wrappers.mrparse.script.mrparse_report:mrparse_report",
    "mrparse_simple": "ccp4i2.pipelines.dr_mr_modelbuild_pipeline.wrappers.mrparse_simple.script.mrparse_simple_report:mrparse_simple_report",
    "pairef": "ccp4i2.wrappers.pairef.script.pairef_report:pairef_report",
    "parrot": "ccp4i2.wrappers.parrot.script.parrot_report:parrot_report",
    "pdb_redo_api": "ccp4i2.wrappers.pdb_redo_api.script.pdb_redo_api_report:pdb_redo_api_report",
    "pdbset_ui": "ccp4i2.wrappers.pdbset_ui.script.pdbset_ui_report:pdbset_ui_report",
    "pdbview_edit": "ccp4i2.wrappers.pdbview_edit.script.pdbview_edit_report:pdbview_edit_report",
    "phaser_EP": "ccp4i2.pipelines.phaser_ep.script.phaser_EP_report:phaser_EP_report",
    "phaser_EP_AUTO": "ccp4i2.pipelines.phaser_pipeline.wrappers.phaser_EP_AUTO.script.phaser_EP_AUTO_report:phaser_EP_AUTO_report",
    "phaser_EP_LLG": "ccp4i2.pipelines.phaser_pipeline.wrappers.phaser_EP_LLG.script.phaser_EP_LLG_report:phaser_EP_LLG_report",
    "phaser_MR_AUTO": "ccp4i2.pipelines.phaser_rnp_pipeline.script.phaser_rnp_pipeline_report:phaser_MR_AUTO_report",
    "phaser_MR_FRF": "ccp4i2.pipelines.phaser_pipeline.wrappers.phaser_MR_FRF.script.phaser_MR_FRF_report:phaser_MR_FRF_report",
    "phaser_MR_FTF": "ccp4i2.pipelines.phaser_pipeline.wrappers.phaser_MR_FTF.script.phaser_MR_FTF_report:phaser_MR_FTF_report",
    "phaser_MR_PAK": "ccp4i2.pipelines.phaser_pipeline.wrappers.phaser_MR_PAK.script.phaser_MR_PAK_report:phaser_MR_PAK_report",
    "phaser_MR_RNP": "ccp4i2.pipelines.phaser_pipeline.wrappers.phaser_MR_RNP.script.phaser_MR_RNP_report:phaser_MR_RNP_report",
    "phaser_analysis": "ccp4i2.wrappers.phaser_analysis.script.phaser_analysis_report:phaser_analysis_report",
    "phaser_ensembler": "ccp4i2.wrappers.phaser_ensembler.script.phaser_ensembler_report:phaser_ensembler_report",
    "phaser_phil": "ccp4i2.wrappers.phaser_phil.script.phaser_phil_report:phaser_phil_report",
    "phaser_pipeline": "ccp4i2.pipelines.phaser_simple.script.phaser_simple_report:phaser_pipeline_report",
    "phaser_rnp_pipeline": "ccp4i2.pipelines.phaser_rnp_pipeline.script.phaser_rnp_pipeline_report:phaser_rnp_pipeline_report",
    "phaser_simple": "ccp4i2.pipelines.phaser_simple.script.phaser_simple_report:phaser_simple_report",
    "phaser_singleMR": "ccp4i2.wrappers.phaser_singleMR.script.phaser_singleMR_report:phaser_singleMR_report",
    "pisa_analyse": "ccp4i2.pipelines.pisapipe.wrappers.pisa_analyse.script.pisa_analyse_report:pisa_analyse_report",
    "pisa_xml": "ccp4i2.pipelines.pisapipe.script.pisapipe_report:pisa_xml_report",
    "pisapipe": "ccp4i2.pipelines.pisapipe.script.pisapipe_report:pisapipe_report",
    "pointless": "ccp4i2.pipelines.phaser_rnp_pipeline.script.phaser_rnp_pipeline_report:pointless_report",
    "pointless_reindexToMatch": "ccp4i2.wrappers.pointless_reindexToMatch.script.pointless_reindexToMatch_report:pointless_reindexToMatch_report",
    "privateer": "ccp4i2.wrappers.privateer.script.privateer_report:privateer_report",
    "prosmart": "ccp4i2.wrappers.prosmart.script.prosmart_report:prosmart_report",
    "prosmart_refmac": "ccp4i2.pipelines.prosmart_refmac.script.prosmart_refmac_report:prosmart_refmac_report",
    "pyphaser_mr": "ccp4i2.wrappers.pyphaser_mr.script.pyphaser_mr_report:pyphaser_mr_report",
    "qtpisa": "ccp4i2.wrappers.qtpisa.script.qtpisa_report:qtpisa_report",
    "refmac": "ccp4i2.pipelines.PrepareDeposit.script.PrepareDeposit_report:refmac_report",
    "scaleit": "ccp4i2.wrappers.scaleit.script.scaleit_report:scaleit_report",
    "sculptor": "ccp4i2.wrappers.sculptor.script.sculptor_report:sculptor_report",
    "servalcat": "ccp4i2.wrappers.servalcat.script.servalcat_report:servalcat_report",
    "servalcat_pipe": "ccp4i2.pipelines.servalcat_pipe.script.servalcat_pipe_report:servalcat_pipe_report",
    "sheetbend": "ccp4i2.pipelines.dr_mr_modelbuild_pipeline.script.dr_mr_modelbuild_pipeline_report:sheetbend_report",
    "shelx": "ccp4i2.pipelines.shelx.script.shelx_report:shelx_report",
    "shelxeMR": "ccp4i2.wrappers.shelxeMR.script.shelxeMR_report:shelxeMR_report",
    "slicendice": "ccp4i2.wrappers.slicendice.script.slicendice_report:slicendice_report",
    "splitMtz": "ccp4i2.wrappers.splitMtz.script.splitMtz_report:splitMtz_report",
    "validate_protein": "ccp4i2.wrappers.validate_protein.script.validate_protein_report:validate_protein_report",
    "xia2_dials": "ccp4i2.wrappers.xia2_dials.script.xia2_dials_report:xia2_dials_report",
    "xia2_multiplex": "ccp4i2.wrappers.xia2_multiplex.script.xia2_multiplex_report:xia2_multiplex_report",
    "xia2_ssx_reduce": "ccp4i2.wrappers.xia2_ssx_reduce.script.xia2_ssx_reduce_report:xia2_ssx_reduce_report",
    "xia2_xds": "ccp4i2.wrappers.xia2_xds.script.xia2_xds_report:xia2_xds_report",
    "zanuda": "ccp4i2.wrappers.zanuda.script.zanuda_report:zanuda_report",
}

_CACHE = {}


def get_report_class(task_name: str):
    """
    Get a report class by task name.

    The report is imported lazily on first access and cached.

    Args:
        task_name: Name of the task (e.g., "refmac", "pointless")

    Returns:
        Report class, or None if not found
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
        warnings.warn(f"Failed to import report {task_name}: {e}")
        return None
