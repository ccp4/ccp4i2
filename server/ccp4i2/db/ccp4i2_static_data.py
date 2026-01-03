USE_PERFORMANCE_CLASSES = True

UUIDTYPE = str

PRIVILEGE_NONE = 0
PRIVILEGE_READ = 1
PRIVILEGE_COMMENT = 2
PRIVILEGE_WRITE = 3
PRIVILEGE_DELETE = 4
PRIVILEGE_PROJECT = 5
PRIVILEGE_EXTEND = 6

JOB_STATUS_UNKNOWN = 0
JOB_STATUS_PENDING = 1
JOB_STATUS_QUEUED = 2
JOB_STATUS_RUNNING = 3
JOB_STATUS_INTERRUPTED = 4
JOB_STATUS_FAILED = 5
JOB_STATUS_FINISHED = 6
JOB_STATUS_REMOTE = 7
JOB_STATUS_FILE_HOLDER = 8
JOB_STATUS_TO_DELETE = 9
JOB_STATUS_UNSATISFACTORY = 10

USER_ROLE_UNKNOWN = 0
USER_ROLE_MANAGER = 1
USER_ROLE_OWNER = 2
USER_ROLE_USER = 3
USER_ROLE_REMOVED = 4

JOB_STATUS_TEXT = [
    "Unknown",
    "Pending",
    "Queued",
    "Running",
    "Interrupted",
    "Failed",
    "Finished",
    "Running remotely",
    "File holder",
    "To delete",
    "Unsatisfactory",
]

FINISHED_JOB_STATUS = ["Finished", "Interrupted", "To delete", "Unsatisfactory"]
JOB_EVALUATION_TEXT = ["Unknown", "Best", "Good", "Rejected"]
USER_AGENT_TEXT = ["Unknown", "CCP4i2", "CCP4mg", "Coot"]
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
MINIMTZFILETYPES = [10, 11, 12, 13]
FILE_ROLE_OUT = 0
FILE_ROLE_IN = 1
FILE_ROLE_IMPORT = 2

PATH_FLAG_JOB_DIR = 1
PATH_FLAG_IMPORT_DIR = 2

"""
Do not ever be tempted to renumber or delete anything here. The position of a type in this list is 
critically important.
"""
FILETYPELIST = [
    (0, "Unknown", "File type unknown"),
    (1, "application/CCP4-seq", "Model sequence"),
    (2, "chemical/x-pdb", "Model coordinates"),
    (3, "MultiPDB", "Multiple model coordinates"),
    (4, "application/CCP4-mtz", "Merged experimental data"),
    (5, "application/CCP4-mtz-unmerged", "Unmerged experimental data"),
    (
        6,
        "application/CCP4-unmerged-experimental",
        "Unmerged experimental data any format",
    ),
    (7, "application/CCP4-map", "Electron density map"),
    (8, "application/refmac-dictionary", "Refmac dictionary"),
    (9, "application/refmac-TLS", "Refmac TLS"),
    (10, "application/CCP4-mtz-freerflag", "FreeR flag"),
    (11, "application/CCP4-mtz-observed", "Observed intensities and structure factors"),
    (12, "application/CCP4-mtz-phases", "Phases"),
    (13, "application/CCP4-mtz-map", "Map coefficients"),
    (14, "Dummy", "Dummy"),
    (15, "application/CCP4-seqalign", "Sequence alignment"),
    (16, "application/CCP4-mtz-mini", "Experimental data object"),
    (17, "application/coot-script", "Coot script"),
    (18, "application/refmac-external-restraints", "Refmac external restraints"),
    (19, "application/CCP4-scene", "CCP4mg scene file"),
    (20, "application/CCP4-shelx-FA", "Shelx FA"),
    (21, "application/phaser-sol", "Phaser solutions"),
    (22, "chemical/x-mdl-molfile", "MDL Molfile"),
    (23, "application/iMosflm-xml", "iMosflm data"),
    (24, "application/CCP4-image", "Image file"),
    (25, "application/CCP4-generic-reflections", "Merged reflection data"),
    (26, "application/CCP4-HHPred-alignments", "HHPred sequence search results"),
    (27, "application/CCP4-Blast-alignments", "Blast sequence search results"),
    (28, "chemical/x-pdb-ensemble", "Ensemble model coordinates"),
    (29, "application/CCP4-asu-content", "Asu content"),
    (30, "application/dials-jfile", "Dials json data file"),
    (31, "application/dials-pfile", "Dials pickle data file"),
    (32, "application/phaser-rfile", "Phaser rotation solutions"),
    (33, "application/refmac-keywords", "Refmac5 keyword file"),
    (34, "application/x-pdf", "PDF File"),
    (35, "application/postscript", "Postscript file"),
    (36, "application/EBI-validation-xml", "Validation XML"),
    (37, "chemical/x-cif", "mmCif reflection data"),
]

KEYTYPELIST = [
    (0, "Unknown", "Key type unknown"),
    (1, "RFactor", "R Factor"),
    (2, "RFree", "Free R Factor"),
    (3, "completeness", "model completeness"),
    (4, "spaceGroup", "space group"),
    (5, "highResLimit", "high resolution limit"),
    (6, "rMeas", "Rmeas data consistency"),
    (7, "FOM", "figure of merit of phases"),
    (8, "CFOM", "correlation FOM"),
    (9, "Hand1Score", "Hand 1 score"),
    (10, "Hand2Score", "Hand 2 score"),
    (11, "CC", "correlation coefficient between Fo and Fc"),
    (12, "nAtoms", "number of atoms in model"),
    (13, "nResidues", "number of residues in model"),
    (14, "phaseError", "phase error"),
    (15, "weightedPhaseError", "weighted phase error"),
    (16, "reflectionCorrelation", "reflection correlation"),
    (17, "RMSxyz", "RMS displacement"),
    (18, "cutoff", "Pairef cutoff"),
    (19, "ccHalf", "correlation coefficient between two half datasets"),
]

FILEASSOCIATIONTYPELIST = [
    (0, "Unknown", "File association type unknown"),
    (1, "Observed-Free", "Observed data and FreeR set"),
]

FILEASSOCIATIONROLELIST = [
    (0, "Unknown", "File association role unknown", 0),
    (1, "Observed data", "Observed data", 1),
    (2, "FreeR set", "FreeR set", 1),
]
