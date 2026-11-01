"""The repeatable records of a deposition's metadata form: authors, related
entries, software, and -- the one with structure -- each crystallised
molecule with the segments its construct is made of.

A segment says where a stretch of the construct came from. Natural and
engineered segments reference UniProt; a designed scaffold the PDB already
holds references that entry; a de novo design references nothing and is a
first-class case, not a UniProt entry with the accession left blank. A
chemically synthesised peptide is recorded as such (``_pdbx_entity_src_syn``)
rather than as an expressed gene.

Resolvable by the def.xml class-name lookup because ``prepare_deposition``
declares this module in its ``Task.dataTypes`` (core/tasks.py).
"""
from ccp4i2.core.base_object.class_metadata import content
from ccp4i2.core.CCP4Data import CData

SEGMENT_KINDS = ["natural", "engineered", "designed", "scaffold", "synthetic_peptide",
                 "tag", "linker"]
SEGMENT_KIND_LABELS = [
    "Natural protein (UniProt)",
    "Engineered variant of a natural protein (UniProt + differences)",
    "De novo design (no database reference)",
    "Designed scaffold with a PDB entry",
    "Chemically synthesised peptide",
    "Expression or purification tag",
    "Linker",
]

SOFTWARE_CLASSES = ["data collection", "data reduction", "data scaling", "phasing",
                    "model building", "refinement", "data processing", "other"]


class CDepositionAuthor(CData):
    """One author: surname and initials as the PDB writes them, and an ORCID."""

    class Meta:
        contents_order = ["NAME", "ORCID"]
        qualifiers = {"allowUndefined": True, "guiLabel": "Author"}

    NAME = content("CString", allowUndefined=False, guiLabel="Name",
                   toolTip="Surname, initials: Noble, M.E.M.")
    ORCID = content("CString", allowUndefined=True, guiLabel="ORCID",
                    toolTip="0000-0002-3595-9807 (leave blank if none)")


class CRelatedEntry(CData):
    """A related database entry (the same protein by another method, the
    scaffold alone, an EM map)."""

    class Meta:
        contents_order = ["DB_NAME", "DB_ID", "CONTENT_TYPE", "DETAILS"]
        qualifiers = {"allowUndefined": True, "guiLabel": "Related entry"}

    DB_NAME = content("CString", enumerators=["PDB", "EMDB", "BMRB"], onlyEnumerators=True,
                      default="PDB", guiLabel="Database")
    DB_ID = content("CString", allowUndefined=False, guiLabel="Entry", toolTip="8PVE, EMD-17963")
    CONTENT_TYPE = content(
        "CString", onlyEnumerators=True, default="unspecified", guiLabel="Relationship",
        enumerators=["unspecified", "re-refinement", "split", "associated EM volume",
                     "complete structure", "other"])
    DETAILS = content("CString", allowUndefined=True, guiLabel="Details",
                      toolTip="Same engineered AHIR by cryoEM; the DARPin alone")


class CSourceSegment(CData):
    """A stretch of a construct and where it came from."""

    class Meta:
        contents_order = ["BEG", "END", "KIND", "ORGANISM", "TAXID", "STRAIN", "GENE",
                          "UNIPROT", "PDB_ID", "ENGINEERED", "HOST", "HOST_TAXID", "DETAILS"]
        qualifiers = {"allowUndefined": True, "guiLabel": "Segment"}

    BEG = content("CInt", allowUndefined=False, min=1, guiLabel="From",
                  toolTip="First residue of the segment, in the construct numbering")
    END = content("CInt", allowUndefined=False, min=1, guiLabel="To",
                  toolTip="Last residue of the segment, in the construct numbering")
    KIND = content("CString", onlyEnumerators=True, default="natural", guiLabel="Origin",
                   enumerators=SEGMENT_KINDS, menuText=SEGMENT_KIND_LABELS)
    ORGANISM = content("CString", allowUndefined=True, guiLabel="Organism",
                       toolTip="Scientific name of the source organism (natural, engineered)")
    TAXID = content("CInt", allowUndefined=True, guiLabel="Taxonomy id",
                    toolTip="NCBI taxonomy id of the source organism")
    STRAIN = content("CString", allowUndefined=True, guiLabel="Strain")
    GENE = content("CString", allowUndefined=True, guiLabel="Gene")
    UNIPROT = content("CString", allowUndefined=True, guiLabel="UniProt",
                      toolTip="Accession of the natural protein (Q5SJ03); the sequence is fetched "
                              "and the differences to the construct recorded")
    PDB_ID = content("CString", allowUndefined=True, guiLabel="PDB entry",
                     toolTip="For a designed scaffold the PDB already holds (2JAB)")
    ENGINEERED = content("CString", allowUndefined=True, guiLabel="Mutations",
                         toolTip="Intended substitutions, comma separated: L33H, Q141P")
    HOST = content("CString", allowUndefined=True, default="Escherichia coli",
                   guiLabel="Expression host")
    HOST_TAXID = content("CInt", allowUndefined=True, default=562, guiLabel="Host taxonomy id")
    DETAILS = content("CString", allowUndefined=True, guiLabel="Details")


class CEntityMetadata(CData):
    """One crystallised molecule (an AU-contents sequence): what to call it,
    what it does, and the segments its construct is made of."""

    class Meta:
        contents_order = ["SEQUENCE_NAME", "DESCRIPTION", "COMMON_NAME", "EC", "MUTATION",
                          "DETAILS", "SEGMENTS"]
        qualifiers = {"allowUndefined": True, "guiLabel": "Molecule"}

    SEQUENCE_NAME = content("CString", allowUndefined=True, guiLabel="AU-contents sequence",
                            toolTip="Which sequence of the AU contents this describes")
    DESCRIPTION = content("CString", allowUndefined=False, guiLabel="Name",
                          toolTip="The molecule's name as the PDB should show it: "
                                  "Ketol-acid reductoisomerase (NADP(+))")
    COMMON_NAME = content("CString", allowUndefined=True, guiLabel="Other names",
                          toolTip="Comma separated synonyms")
    EC = content("CString", allowUndefined=True, guiLabel="EC number")
    MUTATION = content("CString", allowUndefined=True, guiLabel="Mutations",
                       toolTip="Summary for the entity record: L33H")
    DETAILS = content("CString", allowUndefined=True, guiLabel="Details",
                      toolTip="How a design was made, what a chimera is made of")
    # A nested list: the class-name lookup content() uses only knows core/
    # classes, so the item class is given to the CList as its subItem
    # qualifier, the way the def.xml handler would.
    SEGMENTS = content("CList", listMinLength=0, guiLabel="Segments",
                       subItem={"class": CSourceSegment, "qualifiers": {}})


class CSoftwareEntry(CData):
    """A program used, beyond those the lineage already records."""

    class Meta:
        contents_order = ["NAME", "VERSION", "CLASSIFICATION"]
        qualifiers = {"allowUndefined": True, "guiLabel": "Program"}

    NAME = content("CString", allowUndefined=False, guiLabel="Program",
                   toolTip="RFdiffusion, ProteinMPNN, AlphaFold2, Phaser, Coot ...")
    VERSION = content("CString", allowUndefined=True, guiLabel="Version")
    CLASSIFICATION = content("CString", onlyEnumerators=True, default="other",
                             enumerators=SOFTWARE_CLASSES, guiLabel="Used for")
