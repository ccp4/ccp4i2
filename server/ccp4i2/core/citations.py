"""Which references a task cites, and which tasks cite nothing.

Kept here, free of any Django or CCP4 import, because two callers need it:
``lib.utils.jobs.bibliography`` (the Bibliography button) and
``report.metadata.ReferenceGroup`` (the per-report bibliography), and the
latter must stay importable without a configured Django.

A task belongs in NON_CITABLE when it genuinely has nothing to cite - i2
plumbing, format shims, glue.  It is *not* a place to silence a task whose
citation simply lives under another key: those belong in TASK_CITES.

Both live here rather than in ``bibliography`` so that the two consumers
share one map. They used to diverge -- the Bibliography button expanded
aliases and the per-report bibliography did not -- which left 52 of 143
report classes with an empty references section for programs that are
perfectly well cited elsewhere in the same application.
"""

from typing import Dict, List

NON_CITABLE = frozenset({
    # imports / providers — pure CCP4i2 data-plumbing
    "ImportAsuContent", "ImportCoordinate", "ImportDictionary", "ImportFreeR",
    "ImportMap", "ImportMapCoeffs", "ImportObs", "ImportPhases",
    "ImportSequence", "ProvideAlignment", "ProvideAsuContents",
    "ProvideSequence", "ProvideTLS",
    # format converters / column shims — no citable upstream program
    "coordinate_selector", "splitMtz", "mergeMtz", "cad_copy_column",
    "TestObsConversions", "adding_stats_to_mmcif_i2", "cif2mtz", "convert2mtz",
    "x2mtz", "scalepack2mtz", "mtzutils", "mtzheader", "hklin2cif",
    "add_fractional_coords", "editbfac", "pdbview_edit", "chltofom",
    "cmapcoeff", "coot_script_lines",
    # i2 wrapper/glue tasks with no distinct publication
    "MakeLink", "MakeMonster", "MakeProjectsAndDoLigandPipeline",
    "SubtractNative", "pdb_extract_wrapper", "PrepareDeposit", "Lidia",
    # i2 utility tasks with no distinct upstream program paper
    "density_calculator", "dm_multidomain", "findmyseq",
    # subjob-composing pipelines: their constituents surface via Job.children
    # burrowing, so the pipeline name itself needs no citation of its own.
    "prosmart_refmac", "import_merged",
})


TASK_CITES: Dict[str, List[str]] = {
    # --- crank2: an experimental-phasing pipeline that drives many programs
    # internally (see pipelines/crank2/crank2/programs/references/*.nbib). Every
    # crank2 step task cites the full toolchain it can invoke.
    "crank2": [
        "crank2", "shelxc", "shelxd", "shelxe", "refmac", "parrot",
        "buccaneer", "bp3", "afro", "crunch2", "prasa", "solomon", "multicomb",
    ],
    "crank2_comb_phdmmb": ["crank2", "parrot", "buccaneer", "solomon", "multicomb"],
    "crank2_createfree": ["crank2"],
    "crank2_dmfull": ["crank2", "parrot", "solomon"],
    "crank2_faest": ["crank2", "afro"],
    "crank2_handdet": ["crank2"],
    "crank2_mbref": ["crank2", "buccaneer", "refmac"],
    "crank2_phas": ["crank2", "bp3", "refmac"],
    "crank2_phdmmb": ["crank2", "parrot", "buccaneer"],
    "crank2_ref": ["crank2", "refmac"],
    "crank2_refatompick": ["crank2", "refmac"],
    "crank2_substrdet": ["crank2", "shelxd", "prasa", "crunch2"],
    # --- modelcraft: iterative autobuild driving buccaneer / nautilus (nucleic
    # acids) / sheetbend / parrot / refmac internally.
    "modelcraft": [
        "modelcraft", "buccaneer", "nautilus", "sheetbend", "parrot", "refmac",
    ],
    # --- MR pipelines that drive MR + refinement programs internally.
    "mrbump_basic": ["mrbump", "phaser", "molrep", "refmac"],
    "mrparse": ["mrparse", "phaser", "molrep"],
    "mrparse_simple": ["mrparse", "phaser", "molrep"],
    "slicendice": ["slicendice", "phaser", "refmac"],
    "arcimboldo": ["arcimboldo", "shelx", "phaser", "prosmart", "refmac"],
    # --- phaser family: all modes cite the core Phaser paper.
    "phaser_MR": ["phaser"],
    "phaser_MR_AUTO": ["phaser"],
    "phaser_mr_auto_phil": ["phaser"],
    "phaser_ep_auto_phil": ["phaser"],
    "phaser_mr_rnp_phil": ["phaser"],
    "phaser_ep_phil": ["phaser", "parrot", "shelxc", "shelxd", "modelcraft"],
    "phaser_pipeline_phil": ["phaser", "sheetbend", "refmac"],
    "phaser_rnp_pipeline_phil": ["phaser", "pointless", "sheetbend", "refmac"],
    "phaser_simple_phil": ["phaser", "sheetbend", "refmac"],
    "phaser_MR_RNP": ["phaser"],
    "phaser_MR_FRF": ["phaser"],
    "phaser_MR_FTF": ["phaser"],
    "phaser_MR_PAK": ["phaser"],
    "phaser_EP": ["phaser"],
    "phaser_EP_AUTO": ["phaser"],
    "phaser_EP_LLG": ["phaser"],
    "phaser_pipeline": ["phaser"],
    "phaser_rnp_pipeline": ["phaser"],
    "phaser_simple": ["phaser"],
    "phaser_singleMR": ["phaser"],
    "phaser_ensembler": ["phaser"],
    "phaser_analysis": ["phaser"],
    "phaser_phil": ["phaser"],
    "phasertng_picard": ["phasertng"],
    "phasertng_riker": ["phasertng"],
    # --- pisa family: all cite the PISA paper.
    "pisa": ["pisa"],
    "pisa_analyse": ["pisa"],
    "pisa_list": ["pisa"],
    "pisa_xml": ["pisa"],
    "qtpisa": ["pisa"],
    "pisapipe": ["pisa"],
    # --- xia2 family: all cite xia2 (+ DIALS where DIALS is the engine).
    "xia2_run": ["xia2"],
    "xia2_dials": ["xia2", "dials"],
    "xia2_xds": ["xia2"],
    "xia2_ssx_reduce": ["xia2", "dials"],
    "xia2_multiplex": ["xia2", "dials"],
    "xia2_integration": ["xia2", "dials"],
    "xia2_aimless": ["xia2", "aimless"],
    "xia2_ctruncate": ["xia2", "ctruncate"],
    "xia2_pointless": ["xia2", "pointless"],
    "AlternativeImportXIA2": ["xia2"],
    "import_xia2": ["xia2"],
    # --- dials image/lattice tools cite DIALS.
    "dials_image": ["dials"],
    "dials_rlattice": ["dials"],
    "dnatco": ["dnatco"],
    "dnatco_pipe": ["dnatco"],
    "dui": ["dials"],
    # --- mosflm variants.
    "imosflm": ["mosflm"],
    "import_mosflm": ["mosflm"],
    # --- shelx variants.
    "shelx": ["shelx"],
    "shelxeMR": ["shelxe"],
    "ShelxCD": ["shelxc", "shelxd"],
    # --- coot variants cite coot.
    "coot1": ["coot"],
    "coot_rsr_morph": ["coot"],
    "coot_script_lines": ["coot"],
    # --- acedrg variants.
    "acedrgNew": ["acedrg"],
    "AcedrgLink": ["acedrg"],
    # --- molrep variants.
    "molrep_den": ["molrep"],
    "molrep_selfrot": ["molrep"],
    # --- serial-crystallography import pipeline uses DIALS.
    "import_serial": ["dials"],
    "import_serial_pipe": ["dials"],
    # --- dimple difference-map wrapper.
    "i2Dimple": ["dimple"],
    # --- thin variants aliasing to a canonical program.
    "pointless_reindexToMatch": ["pointless"],
    "mrbump_model_prep": ["mrbump"],
    "tableone": ["phaser"],
}
