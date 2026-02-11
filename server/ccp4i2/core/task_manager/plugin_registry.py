"""
This file provides lazy loading
of plugin classes using explicit import statements.
"""

from typing import Optional, Type, Dict, Any


def _get_plugin_class(plugin_name: str) -> Optional[Type]:
    """
    Get a plugin class by name using explicit imports.

    This function uses explicit import statements for each plugin,
    providing clear traceability and IDE support.
    """
    if plugin_name == 'AMPLE':
        from ccp4i2.wrappers.AMPLE.script.AMPLE import AMPLE
        return AMPLE
    if plugin_name == 'AUSPEX':
        from ccp4i2.wrappers.AUSPEX.script.auspex_wrapper import AUSPEX
        return AUSPEX
    if plugin_name == 'AcedrgLink':
        from ccp4i2.wrappers.AcedrgLink.script.AcedrgLink import AcedrgLink
        return AcedrgLink
    if plugin_name == 'AlternativeImportXIA2':
        from ccp4i2.wrappers.AlternativeImportXIA2.script.AlternativeImportXIA2 import AlternativeImportXIA2
        return AlternativeImportXIA2
    if plugin_name == 'Lidia':
        from ccp4i2.wrappers.Lidia.script.Lidia import Lidia
        return Lidia
    if plugin_name == 'LidiaAcedrgNew':
        from ccp4i2.pipelines.LidiaAcedrgNew.script.LidiaAcedrgNew import LidiaAcedrgNew
        return LidiaAcedrgNew
    if plugin_name == 'MakeLink':
        from ccp4i2.pipelines.MakeLink.script.MakeLink import MakeLink
        return MakeLink
    if plugin_name == 'MakeMonster':
        from ccp4i2.wrappers.MakeMonster.script.MakeMonster import MakeMonster
        return MakeMonster
    if plugin_name == 'MakeProjectsAndDoLigandPipeline':
        from ccp4i2.pipelines.MakeProjectsAndDoLigandPipeline.script.MakeProjectsAndDoLigandPipeline import MakeProjectsAndDoLigandPipeline
        return MakeProjectsAndDoLigandPipeline
    if plugin_name == 'Platonyzer':
        from ccp4i2.wrappers.Platonyzer.script.Platonyzer import Platonyzer
        return Platonyzer
    if plugin_name == 'PrepareDeposit':
        from ccp4i2.pipelines.PrepareDeposit.script.PrepareDeposit import PrepareDeposit
        return PrepareDeposit
    if plugin_name == 'ProvideAlignment':
        from ccp4i2.wrappers.ProvideAlignment.script.ProvideAlignment import ProvideAlignment
        return ProvideAlignment
    if plugin_name == 'ProvideAsuContents':
        from ccp4i2.wrappers.ProvideAsuContents.script.ProvideAsuContents import ProvideAsuContents
        return ProvideAsuContents
    if plugin_name == 'ProvideSequence':
        from ccp4i2.wrappers.ProvideSequence.script.ProvideSequence import ProvideSequence
        return ProvideSequence
    if plugin_name == 'ProvideTLS':
        from ccp4i2.wrappers.ProvideTLS.script.ProvideTLS import ProvideTLS
        return ProvideTLS
    if plugin_name == 'SIMBAD':
        from ccp4i2.wrappers.SIMBAD.script.SIMBAD import SIMBAD
        return SIMBAD
    if plugin_name == 'ShelxCD':
        from ccp4i2.wrappers.ShelxCDE.script.ShelxCD import ShelxCD
        return ShelxCD
    if plugin_name == 'SubstituteLigand':
        from ccp4i2.pipelines.SubstituteLigand.script.SubstituteLigand import SubstituteLigand
        return SubstituteLigand
    if plugin_name == 'SubtractNative':
        from ccp4i2.wrappers.SubtractNative.script.SubtractNative import SubtractNative
        return SubtractNative
    if plugin_name == 'TestObsConversions':
        from ccp4i2.wrappers.TestObsConversions.script.TestObsConversions import TestObsConversions
        return TestObsConversions
    if plugin_name == 'acedrg':
        from ccp4i2.wrappers.acedrg.script.acedrg import acedrg
        return acedrg
    if plugin_name == 'acedrgNew':
        from ccp4i2.wrappers.acedrgNew.script.acedrgNew import acedrgNew
        return acedrgNew
    if plugin_name == 'acorn':
        from ccp4i2.wrappers.acorn.script.acorn import acorn
        return acorn
    if plugin_name == 'add_fractional_coords':
        from ccp4i2.wrappers.add_fractional_coords.script.add_fractional_coords import add_fractional_coords
        return add_fractional_coords
    if plugin_name == 'adding_stats_to_mmcif_i2':
        from ccp4i2.wrappers.adding_stats_to_mmcif_i2.script.adding_stats_to_mmcif_i2 import adding_stats_to_mmcif_i2
        return adding_stats_to_mmcif_i2
    if plugin_name == 'aimless':
        from ccp4i2.wrappers.aimless.script.aimless import aimless
        return aimless
    if plugin_name == 'aimless_pipe':
        from ccp4i2.pipelines.aimless_pipe.script.aimless_pipe import aimless_pipe
        return aimless_pipe
    if plugin_name == 'arcimboldo':
        from ccp4i2.wrappers.arcimboldo.script.arcimboldo import arcimboldo
        return arcimboldo
    if plugin_name == 'arp_warp_classic':
        from ccp4i2.wrappers.arp_warp_classic.script.arp_warp_classic import arp_warp_classic
        return arp_warp_classic
    if plugin_name == 'baverage':
        from ccp4i2.wrappers.baverage.script.baverage import baverage
        return baverage
    if plugin_name == 'buster':
        from ccp4i2.wrappers.buster.script.buster import buster
        return buster
    if plugin_name == 'cad_copy_column':
        from ccp4i2.wrappers.cad_copy_column.script.cad_copy_column import cad_copy_column
        return cad_copy_column
    if plugin_name == 'ccp4mg_edit_model':
        from ccp4i2.wrappers.ccp4mg_edit_model.script.ccp4mg_edit_model import ccp4mg_edit_model
        return ccp4mg_edit_model
    if plugin_name == 'ccp4mg_edit_nomrbump':
        from ccp4i2.wrappers.ccp4mg_edit_nomrbump.script.ccp4mg_edit_nomrbump import ccp4mg_edit_nomrbump
        return ccp4mg_edit_nomrbump
    if plugin_name == 'ccp4mg_general':
        from ccp4i2.wrappers.ccp4mg_general.script.ccp4mg_general import ccp4mg_general
        return ccp4mg_general
    if plugin_name == 'chainsaw':
        from ccp4i2.wrappers.chainsaw.script.chainsaw import chainsaw
        return chainsaw
    if plugin_name == 'chltofom':
        from ccp4i2.wrappers.chltofom.script.chltofom import chltofom
        return chltofom
    if plugin_name == 'cif2mtz':
        from ccp4i2.wrappers.cif2mtz.script.cif2mtz import cif2mtz
        return cif2mtz
    if plugin_name == 'clustalw':
        from ccp4i2.wrappers.clustalw.script.clustalw import clustalw
        return clustalw
    if plugin_name == 'cmapcoeff':
        from ccp4i2.wrappers.cmapcoeff.script.cmapcoeff import cmapcoeff
        return cmapcoeff
    if plugin_name == 'comit':
        from ccp4i2.wrappers.comit.script.comit import comit
        return comit
    if plugin_name == 'convert2mtz':
        from ccp4i2.wrappers.convert2mtz.script.convert2mtz import convert2mtz
        return convert2mtz
    if plugin_name == 'coordinate_selector':
        from ccp4i2.wrappers.coordinate_selector.script.coordinate_selector import coordinate_selector
        return coordinate_selector
    if plugin_name == 'coot1':
        from ccp4i2.wrappers.coot1.script.coot1 import coot1
        return coot1
    if plugin_name == 'coot_find_ligand':
        from ccp4i2.wrappers.coot_find_ligand.script.coot_find_ligand import coot_find_ligand
        return coot_find_ligand
    if plugin_name == 'coot_find_waters':
        from ccp4i2.wrappers.coot_find_waters.script.coot_find_waters import coot_find_waters
        return coot_find_waters
    if plugin_name == 'coot_rebuild':
        from ccp4i2.wrappers.coot_rebuild.script.coot_rebuild import coot_rebuild
        return coot_rebuild
    if plugin_name == 'coot_rsr_morph':
        from ccp4i2.wrappers.coot_rsr_morph.script.coot_rsr_morph import coot_rsr_morph
        return coot_rsr_morph
    if plugin_name == 'coot_script_lines':
        from ccp4i2.wrappers.coot_script_lines.script.coot_script_lines import coot_script_lines
        return coot_script_lines
    if plugin_name == 'cpatterson':
        from ccp4i2.wrappers.cpatterson.script.cpatterson import cpatterson
        return cpatterson
    if plugin_name == 'cphasematch':
        from ccp4i2.wrappers.cphasematch.script.cphasematch import cphasematch
        return cphasematch
    if plugin_name == 'crank2':
        from ccp4i2.pipelines.crank2.script.crank2_script import crank2
        return crank2
    if plugin_name == 'crank2_comb_phdmmb':
        from ccp4i2.pipelines.crank2.wrappers.crank2_comb_phdmmb.script.crank2_comb_phdmmb import crank2_comb_phdmmb
        return crank2_comb_phdmmb
    if plugin_name == 'crank2_createfree':
        from ccp4i2.pipelines.crank2.wrappers.crank2_createfree.script.crank2_createfree import crank2_createfree
        return crank2_createfree
    if plugin_name == 'crank2_dmfull':
        from ccp4i2.pipelines.crank2.wrappers.crank2_dmfull.script.crank2_dmfull import crank2_dmfull
        return crank2_dmfull
    if plugin_name == 'crank2_faest':
        from ccp4i2.pipelines.crank2.wrappers.crank2_faest.script.crank2_faest import crank2_faest
        return crank2_faest
    if plugin_name == 'crank2_handdet':
        from ccp4i2.pipelines.crank2.wrappers.crank2_handdet.script.crank2_handdet import crank2_handdet
        return crank2_handdet
    if plugin_name == 'crank2_mbref':
        from ccp4i2.pipelines.crank2.wrappers.crank2_mbref.script.crank2_mbref import crank2_mbref
        return crank2_mbref
    if plugin_name == 'crank2_phas':
        from ccp4i2.pipelines.crank2.wrappers.crank2_phas.script.crank2_phas import crank2_phas
        return crank2_phas
    if plugin_name == 'crank2_phdmmb':
        from ccp4i2.pipelines.crank2.wrappers.crank2_phdmmb.script.crank2_phdmmb import crank2_phdmmb
        return crank2_phdmmb
    if plugin_name == 'crank2_ref':
        from ccp4i2.pipelines.crank2.wrappers.crank2_ref.script.crank2_ref import crank2_ref
        return crank2_ref
    if plugin_name == 'crank2_refatompick':
        from ccp4i2.pipelines.crank2.wrappers.crank2_refatompick.script.crank2_refatompick import crank2_refatompick
        return crank2_refatompick
    if plugin_name == 'crank2_substrdet':
        from ccp4i2.pipelines.crank2.wrappers.crank2_substrdet.script.crank2_substrdet import crank2_substrdet
        return crank2_substrdet
    if plugin_name == 'csymmatch':
        from ccp4i2.wrappers.csymmatch.script.csymmatch import csymmatch
        return csymmatch
    if plugin_name == 'ctruncate':
        from ccp4i2.wrappers.ctruncate.script.ctruncate import ctruncate
        return ctruncate
    if plugin_name == 'density_calculator':
        from ccp4i2.wrappers.density_calculator.script.density_calculator import density_calculator
        return density_calculator
    if plugin_name == 'dials_image':
        from ccp4i2.wrappers.dials_image.script.dials_image import dials_image
        return dials_image
    if plugin_name == 'dials_rlattice':
        from ccp4i2.wrappers.dials_rlattice.script.dials_rlattice import dials_rlattice
        return dials_rlattice
    if plugin_name == 'dr_mr_modelbuild_pipeline':
        from ccp4i2.pipelines.dr_mr_modelbuild_pipeline.script.dr_mr_modelbuild_pipeline import dr_mr_modelbuild_pipeline
        return dr_mr_modelbuild_pipeline
    if plugin_name == 'dui':
        from ccp4i2.wrappers.dui.script.dui import dui
        return dui
    if plugin_name == 'editbfac':
        from ccp4i2.wrappers.editbfac.script.editbfac import editbfac
        return editbfac
    if plugin_name == 'edstats':
        from ccp4i2.wrappers.edstats.script.edstats import edstats
        return edstats
    if plugin_name == 'fft':
        from ccp4i2.wrappers.fft.script.fft import fft
        return fft
    if plugin_name == 'findmyseq':
        from ccp4i2.wrappers.findmyseq.script.findmyseq import findmyseq
        return findmyseq
    if plugin_name == 'freerflag':
        from ccp4i2.wrappers.freerflag.script.freerflag import freerflag
        return freerflag
    if plugin_name == 'gesamt':
        from ccp4i2.wrappers.gesamt.script.gesamt import gesamt
        return gesamt
    if plugin_name == 'hklin2cif':
        from ccp4i2.pipelines.PrepareDeposit.wrappers.hklin2cif.script.hklin2cif import hklin2cif
        return hklin2cif
    if plugin_name == 'i2Dimple':
        from ccp4i2.wrappers.i2Dimple.script.i2Dimple import i2Dimple
        return i2Dimple
    if plugin_name == 'imosflm':
        from ccp4i2.wrappers.imosflm.script.imosflm import imosflm
        return imosflm
    if plugin_name == 'import_merged':
        from ccp4i2.pipelines.import_merged.script.import_merged import import_merged
        return import_merged
    if plugin_name == 'import_mosflm':
        from ccp4i2.wrappers.import_mosflm.script.import_mosflm import import_mosflm
        return import_mosflm
    if plugin_name == 'import_serial':
        from ccp4i2.wrappers.import_serial.script.import_serial import import_serial
        return import_serial
    if plugin_name == 'import_serial_pipe':
        from ccp4i2.pipelines.import_serial_pipe.script.import_serial_pipe import import_serial_pipe
        return import_serial_pipe
    if plugin_name == 'import_xia2':
        from ccp4i2.pipelines.import_xia2.script.import_xia2 import import_xia2
        return import_xia2
    if plugin_name == 'lorestr_i2':
        from ccp4i2.wrappers.lorestr_i2.script.lorestr_i2 import lorestr_i2
        return lorestr_i2
    if plugin_name == 'mergeMtz':
        from ccp4i2.wrappers.mergeMtz.script.mergeMtz import mergeMtz
        return mergeMtz
    if plugin_name == 'metalCoord':
        from ccp4i2.wrappers.metalCoord.script.metalCoord import metalCoord
        return metalCoord
    if plugin_name == 'modelASUCheck':
        from ccp4i2.wrappers.modelASUCheck.script.modelASUCheck import modelASUCheck
        return modelASUCheck
    if plugin_name == 'modelcraft':
        from ccp4i2.wrappers.modelcraft.script.modelcraft import modelcraft
        return modelcraft
    if plugin_name == 'molrep_den':
        from ccp4i2.wrappers.molrep_den.script.molrep_den import molrep_den
        return molrep_den
    if plugin_name == 'molrep_mr':
        from ccp4i2.wrappers.molrep_mr.script.molrep_mr import molrep_mr
        return molrep_mr
    if plugin_name == 'molrep_pipe':
        from ccp4i2.pipelines.molrep_pipe.script.molrep_pipe import molrep_pipe
        return molrep_pipe
    if plugin_name == 'molrep_selfrot':
        from ccp4i2.wrappers.molrep_selfrot.script.molrep_selfrot import molrep_selfrot
        return molrep_selfrot
    if plugin_name == 'morda_i2':
        from ccp4i2.wrappers.morda_i2.script.morda_i2 import morda_i2
        return morda_i2
    if plugin_name == 'mosflm':
        from ccp4i2.wrappers.mosflm.script.mosflm import mosflm
        return mosflm
    if plugin_name == 'mrbump_basic':
        from ccp4i2.wrappers.mrbump_basic.script.mrbump_basic import mrbump_basic
        return mrbump_basic
    if plugin_name == 'mrbump_model_prep':
        from ccp4i2.pipelines.dr_mr_modelbuild_pipeline.wrappers.mrbump_model_prep.script.mrbump_model_prep import mrbump_model_prep
        return mrbump_model_prep
    if plugin_name == 'mrparse':
        from ccp4i2.wrappers.mrparse.script.mrparse_wrapper import mrparse
        return mrparse
    if plugin_name == 'mrparse_simple':
        from ccp4i2.pipelines.dr_mr_modelbuild_pipeline.wrappers.mrparse_simple.script.mrparse_simple_wrapper import mrparse_simple
        return mrparse_simple
    if plugin_name == 'mtzheader':
        from ccp4i2.wrappers.mtzheader.script.mtzheader import mtzheader
        return mtzheader
    if plugin_name == 'mtzutils':
        from ccp4i2.wrappers.mtzutils.script.mtzutils import mtzutils
        return mtzutils
    if plugin_name == 'pairef':
        from ccp4i2.wrappers.pairef.script.pairef import pairef
        return pairef
    if plugin_name == 'parrot':
        from ccp4i2.wrappers.parrot.script.parrot import parrot
        return parrot
    if plugin_name == 'pdb_extract_wrapper':
        from ccp4i2.pipelines.PrepareDeposit.wrappers.pdb_extract_wrapper.script.pdb_extract_wrapper import pdb_extract_wrapper
        return pdb_extract_wrapper
    if plugin_name == 'pdb_redo_api':
        from ccp4i2.wrappers.pdb_redo_api.script.pdb_redo_api import pdb_redo_api
        return pdb_redo_api
    if plugin_name == 'pdbset_ui':
        from ccp4i2.wrappers.pdbset_ui.script.pdbset_ui import pdbset_ui
        return pdbset_ui
    if plugin_name == 'pdbview_edit':
        from ccp4i2.wrappers.pdbview_edit.script.pdbview_edit import pdbview_edit
        return pdbview_edit
    if plugin_name == 'phaser_EP':
        from ccp4i2.pipelines.phaser_ep.script.phaser_EP import phaser_EP
        return phaser_EP
    if plugin_name == 'phaser_EP_AUTO':
        from ccp4i2.pipelines.phaser_pipeline.wrappers.phaser_EP_AUTO.script.phaser_EP_AUTO import phaser_EP_AUTO
        return phaser_EP_AUTO
    if plugin_name == 'phaser_EP_LLG':
        from ccp4i2.pipelines.phaser_pipeline.wrappers.phaser_EP_LLG.script.phaser_EP_LLG import phaser_EP_LLG
        return phaser_EP_LLG
    if plugin_name == 'phaser_MR':
        from ccp4i2.pipelines.phaser_pipeline.wrappers.phaser_MR.script.phaser_MR import phaser_MR
        return phaser_MR
    if plugin_name == 'phaser_MR_AUTO':
        from ccp4i2.pipelines.phaser_pipeline.wrappers.phaser_MR_AUTO.script.phaser_MR_AUTO import phaser_MR_AUTO
        return phaser_MR_AUTO
    if plugin_name == 'phaser_MR_FRF':
        from ccp4i2.pipelines.phaser_pipeline.wrappers.phaser_MR_FRF.script.phaser_MR_FRF import phaser_MR_FRF
        return phaser_MR_FRF
    if plugin_name == 'phaser_MR_FTF':
        from ccp4i2.pipelines.phaser_pipeline.wrappers.phaser_MR_FTF.script.phaser_MR_FTF import phaser_MR_FTF
        return phaser_MR_FTF
    if plugin_name == 'phaser_MR_PAK':
        from ccp4i2.pipelines.phaser_pipeline.wrappers.phaser_MR_PAK.script.phaser_MR_PAK import phaser_MR_PAK
        return phaser_MR_PAK
    if plugin_name == 'phaser_MR_RNP':
        from ccp4i2.pipelines.phaser_pipeline.wrappers.phaser_MR_RNP.script.phaser_MR_RNP import phaser_MR_RNP
        return phaser_MR_RNP
    if plugin_name == 'phaser_analysis':
        from ccp4i2.wrappers.phaser_analysis.script.phaser_analysis import phaser_analysis
        return phaser_analysis
    if plugin_name == 'phaser_ensembler':
        from ccp4i2.wrappers.phaser_ensembler.script.phaser_ensembler import phaser_ensembler
        return phaser_ensembler
    if plugin_name == 'phaser_phil':
        from ccp4i2.wrappers.phaser_phil.script.phaser_phil import phaser_phil
        return phaser_phil
    if plugin_name == 'phaser_pipeline':
        from ccp4i2.pipelines.phaser_pipeline.script.phaser_pipeline import phaser_pipeline
        return phaser_pipeline
    if plugin_name == 'phaser_rnp_pipeline':
        from ccp4i2.pipelines.phaser_rnp_pipeline.script.phaser_rnp_pipeline import phaser_rnp_pipeline
        return phaser_rnp_pipeline
    if plugin_name == 'phaser_simple':
        from ccp4i2.pipelines.phaser_simple.script.phaser_simple import phaser_simple
        return phaser_simple
    if plugin_name == 'phaser_singleMR':
        from ccp4i2.wrappers.phaser_singleMR.script.phaser_singleMR import phaser_singleMR
        return phaser_singleMR
    if plugin_name == 'pisa_analyse':
        from ccp4i2.pipelines.pisapipe.wrappers.pisa_analyse.script.pisa_analyse import pisa_analyse
        return pisa_analyse
    if plugin_name == 'pisa_list':
        from ccp4i2.pipelines.pisapipe.wrappers.pisa_list.script.pisa_list import pisa_list
        return pisa_list
    if plugin_name == 'pisa_xml':
        from ccp4i2.pipelines.pisapipe.wrappers.pisa_xml.script.pisa_xml import pisa_xml
        return pisa_xml
    if plugin_name == 'pisapipe':
        from ccp4i2.pipelines.pisapipe.script.pisapipe import pisapipe
        return pisapipe
    if plugin_name == 'pointless':
        from ccp4i2.wrappers.pointless.script.pointless import pointless
        return pointless
    if plugin_name == 'pointless_reindexToMatch':
        from ccp4i2.wrappers.pointless_reindexToMatch.script.pointless_reindexToMatch import pointless_reindexToMatch
        return pointless_reindexToMatch
    if plugin_name == 'privateer':
        from ccp4i2.wrappers.privateer.script.privateer_wrapper import privateer
        return privateer
    if plugin_name == 'prosmart':
        from ccp4i2.wrappers.prosmart.script.prosmart import prosmart
        return prosmart
    if plugin_name == 'prosmart_refmac':
        from ccp4i2.pipelines.prosmart_refmac.script.prosmart_refmac import prosmart_refmac
        return prosmart_refmac
    if plugin_name == 'pyphaser_mr':
        from ccp4i2.wrappers.pyphaser_mr.script.pyphaser_mr import pyphaser_mr
        return pyphaser_mr
    if plugin_name == 'qtpisa':
        from ccp4i2.wrappers.qtpisa.script.qtpisa import qtpisa
        return qtpisa
    if plugin_name == 'refmac':
        from ccp4i2.wrappers.refmac.script.refmac import refmac
        return refmac
    if plugin_name == 'scaleit':
        from ccp4i2.wrappers.scaleit.script.scaleit import scaleit
        return scaleit
    if plugin_name == 'scalepack2mtz':
        from ccp4i2.wrappers.scalepack2mtz.script.scalepack2mtz import scalepack2mtz
        return scalepack2mtz
    if plugin_name == 'sculptor':
        from ccp4i2.wrappers.sculptor.script.sculptor import sculptor
        return sculptor
    if plugin_name == 'servalcat':
        from ccp4i2.wrappers.servalcat.script.servalcat import servalcat
        return servalcat
    if plugin_name == 'servalcat_pipe':
        from ccp4i2.pipelines.servalcat_pipe.script.servalcat_pipe import servalcat_pipe
        return servalcat_pipe
    if plugin_name == 'sheetbend':
        from ccp4i2.wrappers.sheetbend.script.sheetbend import sheetbend
        return sheetbend
    if plugin_name == 'shelx':
        from ccp4i2.pipelines.shelx.script.shelx_script import shelx
        return shelx
    if plugin_name == 'shelxeMR':
        from ccp4i2.wrappers.shelxeMR.script.shelxeMR import shelxeMR
        return shelxeMR
    if plugin_name == 'slicendice':
        from ccp4i2.wrappers.slicendice.script.slicendice import slicendice
        return slicendice
    if plugin_name == 'splitMtz':
        from ccp4i2.wrappers.splitMtz.script.splitMtz import splitMtz
        return splitMtz
    if plugin_name == 'tableone':
        from ccp4i2.pipelines.tableone.script.tableone import tableone
        return tableone
    if plugin_name == 'validate_protein':
        from ccp4i2.wrappers.validate_protein.script.validate_protein import validate_protein
        return validate_protein
    if plugin_name == 'x2mtz':
        from ccp4i2.wrappers.x2mtz.script.x2mtz import x2mtz
        return x2mtz
    if plugin_name == 'xia2_aimless':
        from ccp4i2.pipelines.import_xia2.wrappers.xia2_aimless.script.xia2_aimless import xia2_aimless
        return xia2_aimless
    if plugin_name == 'xia2_ctruncate':
        from ccp4i2.pipelines.import_xia2.wrappers.xia2_ctruncate.script.xia2_ctruncate import xia2_ctruncate
        return xia2_ctruncate
    if plugin_name == 'xia2_dials':
        from ccp4i2.wrappers.xia2_dials.script.xia2_dials import xia2_dials
        return xia2_dials
    if plugin_name == 'xia2_integration':
        from ccp4i2.pipelines.import_xia2.wrappers.xia2_integration.script.xia2_integration import xia2_integration
        return xia2_integration
    if plugin_name == 'xia2_multiplex':
        from ccp4i2.wrappers.xia2_multiplex.script.xia2_multiplex import xia2_multiplex
        return xia2_multiplex
    if plugin_name == 'xia2_pointless':
        from ccp4i2.pipelines.import_xia2.wrappers.xia2_pointless.script.xia2_pointless import xia2_pointless
        return xia2_pointless
    if plugin_name == 'xia2_run':
        from ccp4i2.pipelines.import_xia2.wrappers.xia2_run.script.xia2_run import xia2_run
        return xia2_run
    if plugin_name == 'xia2_ssx_reduce':
        from ccp4i2.wrappers.xia2_ssx_reduce.script.xia2_ssx_reduce import xia2_ssx_reduce
        return xia2_ssx_reduce
    if plugin_name == 'xia2_xds':
        from ccp4i2.wrappers.xia2_xds.script.xia2_xds import xia2_xds
        return xia2_xds
    if plugin_name == 'zanuda':
        from ccp4i2.wrappers.zanuda.script.zanuda import zanuda
        return zanuda
    return None


# Plugin names for fast lookup without loading metadata
PLUGIN_NAMES: set[str] = {
    'AMPLE',
    'AUSPEX',
    'AcedrgLink',
    'AlternativeImportXIA2',
    'Lidia',
    'LidiaAcedrgNew',
    'MakeLink',
    'MakeMonster',
    'MakeProjectsAndDoLigandPipeline',
    'Platonyzer',
    'PrepareDeposit',
    'ProvideAlignment',
    'ProvideAsuContents',
    'ProvideSequence',
    'ProvideTLS',
    'SIMBAD',
    'ShelxCD',
    'SubstituteLigand',
    'SubtractNative',
    'TestObsConversions',
    'acedrg',
    'acedrgNew',
    'acorn',
    'add_fractional_coords',
    'adding_stats_to_mmcif_i2',
    'aimless',
    'aimless_pipe',
    'arcimboldo',
    'arp_warp_classic',
    'baverage',
    'buster',
    'cad_copy_column',
    'ccp4mg_edit_model',
    'ccp4mg_edit_nomrbump',
    'ccp4mg_general',
    'chainsaw',
    'chltofom',
    'cif2mtz',
    'clustalw',
    'cmapcoeff',
    'comit',
    'convert2mtz',
    'coordinate_selector',
    'coot1',
    'coot_find_ligand',
    'coot_find_waters',
    'coot_rebuild',
    'coot_rsr_morph',
    'coot_script_lines',
    'cpatterson',
    'cphasematch',
    'crank2',
    'crank2_comb_phdmmb',
    'crank2_createfree',
    'crank2_dmfull',
    'crank2_faest',
    'crank2_handdet',
    'crank2_mbref',
    'crank2_phas',
    'crank2_phdmmb',
    'crank2_ref',
    'crank2_refatompick',
    'crank2_substrdet',
    'csymmatch',
    'ctruncate',
    'density_calculator',
    'dials_image',
    'dials_rlattice',
    'dr_mr_modelbuild_pipeline',
    'dui',
    'editbfac',
    'edstats',
    'fft',
    'findmyseq',
    'freerflag',
    'gesamt',
    'hklin2cif',
    'i2Dimple',
    'imosflm',
    'import_merged',
    'import_mosflm',
    'import_serial',
    'import_serial_pipe',
    'import_xia2',
    'lorestr_i2',
    'mergeMtz',
    'metalCoord',
    'modelASUCheck',
    'modelcraft',
    'molrep_den',
    'molrep_mr',
    'molrep_pipe',
    'molrep_selfrot',
    'morda_i2',
    'mosflm',
    'mrbump_basic',
    'mrbump_model_prep',
    'mrparse',
    'mrparse_simple',
    'mtzheader',
    'mtzutils',
    'pairef',
    'parrot',
    'pdb_extract_wrapper',
    'pdb_redo_api',
    'pdbset_ui',
    'pdbview_edit',
    'phaser_EP',
    'phaser_EP_AUTO',
    'phaser_EP_LLG',
    'phaser_MR',
    'phaser_MR_AUTO',
    'phaser_MR_FRF',
    'phaser_MR_FTF',
    'phaser_MR_PAK',
    'phaser_MR_RNP',
    'phaser_analysis',
    'phaser_ensembler',
    'phaser_phil',
    'phaser_pipeline',
    'phaser_rnp_pipeline',
    'phaser_simple',
    'phaser_singleMR',
    'pisa_analyse',
    'pisa_list',
    'pisa_xml',
    'pisapipe',
    'pointless',
    'pointless_reindexToMatch',
    'privateer',
    'prosmart',
    'prosmart_refmac',
    'pyphaser_mr',
    'qtpisa',
    'refmac',
    'scaleit',
    'scalepack2mtz',
    'sculptor',
    'servalcat',
    'servalcat_pipe',
    'sheetbend',
    'shelx',
    'shelxeMR',
    'slicendice',
    'splitMtz',
    'tableone',
    'validate_protein',
    'x2mtz',
    'xia2_aimless',
    'xia2_ctruncate',
    'xia2_dials',
    'xia2_integration',
    'xia2_multiplex',
    'xia2_pointless',
    'xia2_run',
    'xia2_ssx_reduce',
    'xia2_xds',
    'zanuda',
}


# Plugin metadata loaded lazily from JSON
_PLUGIN_METADATA: Optional[Dict[str, Dict[str, Any]]] = None


def _load_metadata() -> Dict[str, Dict[str, Any]]:
    """Load plugin metadata from JSON file."""
    global _PLUGIN_METADATA
    if _PLUGIN_METADATA is None:
        import json
        import os
        script_dir = os.path.dirname(os.path.abspath(__file__))
        json_path = os.path.join(script_dir, "plugin_lookup.json")
        try:
            with open(json_path, "r") as f:
                _PLUGIN_METADATA = json.load(f)
        except Exception:
            _PLUGIN_METADATA = {}
    return _PLUGIN_METADATA


class PluginRegistry:
    """Registry for lazy-loading plugin classes."""

    def __init__(self):
        self._cache: Dict[str, Type] = {}

    def get_plugin_class(self, task_name: str) -> Optional[Type]:
        """
        Get a plugin class by name.

        The plugin is imported lazily on first access and cached.

        Args:
            task_name: Name of the task/plugin

        Returns:
            Plugin class, or None if not found
        """
        cache_key = task_name
        if cache_key in self._cache:
            return self._cache[cache_key]

        if task_name not in PLUGIN_NAMES:
            return None

        try:
            plugin_class = _get_plugin_class(task_name)
            if plugin_class is not None:
                self._cache[cache_key] = plugin_class
            return plugin_class
        except Exception as e:
            import logging
            import traceback
            logger = logging.getLogger(f"ccp4i2:{__name__}")
            logger.error(f"Failed to import plugin {task_name}: {e}")
            logger.error(f"Traceback:\n{traceback.format_exc()}")
            return None


# Singleton instance
_registry = None


def get_registry() -> PluginRegistry:
    """Get the singleton plugin registry instance."""
    global _registry
    if _registry is None:
        _registry = PluginRegistry()
    return _registry
