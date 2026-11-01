"""The controlled vocabularies OneDep's form enforces but the PDBx dictionary
does not enumerate (or enumerates elsewhere).

What the DepUI rejected or dropped on the Ahir3252 deposition (2026-11-01):
a ``_software.name`` outside its software list is dropped silently, and a
combined classification ("OXIDOREDUCTASE/DE NOVO PROTEIN") is "Value not in
the enumeration list". The lists here are what deposited entries carry; the
wwPDB's own software list was not reachable to check, so a program the
depositor needs that is missing here is reported, not written.
"""

# _software.name as OneDep spells it. Keys are lower-case aliases.
SOFTWARE_NAMES = {
    "aimless": "Aimless", "pointless": "POINTLESS", "scala": "SCALA", "truncate": "TRUNCATE",
    "ctruncate": "CTRUNCATE", "dials": "DIALS", "xia2": "xia2", "xds": "XDS", "xscale": "XSCALE",
    "autoproc": "autoPROC", "staraniso": "STARANISO", "mosflm": "MOSFLM", "imosflm": "iMOSFLM",
    "hkl-2000": "HKL-2000", "hkl2000": "HKL-2000", "hkl-3000": "HKL-3000", "hkl3000": "HKL-3000",
    "molrep": "MOLREP", "phaser": "PHASER", "phasertng": "PHASER", "amore": "AMoRe",
    "shelxc": "SHELXC", "shelxd": "SHELXD", "shelxe": "SHELXE", "shelxl": "SHELXL",
    "crank2": "Crank2", "crank": "Crank2", "parrot": "Parrot", "buccaneer": "Buccaneer",
    "modelcraft": "ModelCraft", "arp/warp": "ARP/wARP", "arp_warp": "ARP/wARP",
    "coot": "Coot", "refmac": "REFMAC", "refmac5": "REFMAC", "refmacat": "REFMAC",
    "servalcat": "Servalcat", "phenix": "PHENIX", "phenix.refine": "PHENIX",
    "buster": "BUSTER", "buster-tnt": "BUSTER", "cns": "CNS", "ccp4": "CCP4",
    "pdb_extract": "PDB_EXTRACT", "molprobity": "MolProbity", "pisa": "PISA",
    "alphafold": "AlphaFold", "alphafold2": "AlphaFold", "alphafold3": "AlphaFold",
    "rosettafold": "RoseTTAFold", "rosetta": "Rosetta", "esmfold": "ESMFold",
    "swiss-model": "SWISS-MODEL", "modeller": "MODELLER",
    # Design tools are new to the archive's list; these spellings are the
    # programs' own and are written as given.
    "rfdiffusion": "RFdiffusion", "proteinmpnn": "ProteinMPNN", "chai-1": "Chai-1", "boltz": "Boltz",
}

# Names OneDep knows nothing of: written to the report, not the file.
SOFTWARE_NOT_WRITTEN = {"ccp4i2"}


def canonical_software_name(name):
    """OneDep's spelling of a program name, or None when the name is not in
    the list (CCP4i2 itself, a script, a typo)."""
    key = str(name or "").strip().lower()
    if not key or key in SOFTWARE_NOT_WRITTEN:
        return None
    return SOFTWARE_NAMES.get(key)


# _struct_keywords.pdbx_keywords: one term from the DepUI classification list.
CLASSIFICATIONS = (
    "ANTIBIOTIC", "ANTIFUNGAL PROTEIN", "ANTIMICROBIAL PROTEIN", "ANTITOXIN",
    "ANTITUMOR PROTEIN", "ANTIVIRAL PROTEIN", "APOPTOSIS", "BIOSYNTHETIC PROTEIN",
    "BLOOD CLOTTING", "CELL ADHESION", "CELL CYCLE", "CELL INVASION", "CHAPERONE",
    "CIRCADIAN CLOCK PROTEIN", "CONTRACTILE PROTEIN", "CYTOKINE", "CYTOSOLIC PROTEIN",
    "DE NOVO PROTEIN", "DNA", "DNA BINDING PROTEIN", "DNA/RNA", "ELECTRON TRANSPORT",
    "ENDOCYTOSIS", "EXOCYTOSIS", "FLAVOPROTEIN", "FLUORESCENT PROTEIN", "GENE REGULATION",
    "HORMONE", "HYDROLASE", "IMMUNE SYSTEM", "ISOMERASE", "LIGASE", "LIPID BINDING PROTEIN",
    "LIPID TRANSPORT", "LUMINESCENT PROTEIN", "LYASE", "MEMBRANE PROTEIN",
    "METAL BINDING PROTEIN", "METAL TRANSPORT", "MOTOR PROTEIN", "NEUROPEPTIDE",
    "NUCLEAR PROTEIN", "ONCOPROTEIN", "OXIDOREDUCTASE", "OXYGEN STORAGE", "OXYGEN TRANSPORT",
    "PEPTIDE BINDING PROTEIN", "PHOTOSYNTHESIS", "PLANT PROTEIN", "PROTEIN BINDING",
    "PROTEIN FIBRIL", "PROTEIN TRANSPORT", "PROTON TRANSPORT", "RECOMBINATION",
    "REPLICATION", "RIBOSOME", "RNA", "RNA BINDING PROTEIN", "SIGNALING PROTEIN",
    "SPLICING", "STRUCTURAL GENOMICS", "STRUCTURAL PROTEIN", "SUGAR BINDING PROTEIN",
    "SURFACTANT PROTEIN", "TOXIN", "TRANSCRIPTION", "TRANSFERASE", "TRANSLATION",
    "TRANSLOCASE", "TRANSPORT PROTEIN", "UNKNOWN FUNCTION", "VIRAL PROTEIN", "VIRUS",
    "VIRUS LIKE PARTICLE",
)

# _pdbx_struct_assembly_auth_evidence.experimental_support (dictionary enumeration).
ASSEMBLY_EVIDENCE = (
    "none", "gel filtration", "light scattering", "homology", "equilibrium centrifugation",
    "mass spectrometry", "SAXS", "immunoprecipitation", "isothermal titration calorimetry",
    "cross-linking", "native gel electrophoresis", "microscopy", "assay for oligomerization",
    "scanning transmission electron microscopy", "surface plasmon resonance",
    "fluorescence resonance energy transfer", "NMR relaxation study", "electron microscopy",
)

# _pdbx_initial_refinement_model (dictionary enumerations).
STARTING_MODEL_TYPES = ("experimental model", "in silico model", "integrative model", "other")
STARTING_MODEL_SOURCES = ("PDB", "AlphaFold", "Robetta", "RoseTTAFold", "ModelArchive",
                          "SwissModel", "Modeller", "ITasser", "PDB-Dev", "PHYRE", "Other")
