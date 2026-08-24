"""Tasks that have no citable upstream program.

Kept here, free of any Django or CCP4 import, because two callers need it:
``lib.utils.jobs.bibliography`` (the Bibliography button) and
``report.metadata.ReferenceGroup`` (the per-report bibliography), and the
latter must stay importable without a configured Django.

A task belongs here when it genuinely has nothing to cite - i2 plumbing,
format shims, glue.  It is *not* a place to silence a task whose citation
simply lives under another key: those resolve through the alias map in
``bibliography._TASK_CITATIONS`` instead.
"""

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
