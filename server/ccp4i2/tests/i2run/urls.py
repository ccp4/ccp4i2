_PDBE = "https://www.ebi.ac.uk/pdbe"
_RCSB = "https://files.rcsb.org"
_REDO = "https://pdb-redo.eu/db"
_UNIPROT = "https://rest.uniprot.org"
_EMDB = "https://ftp.ebi.ac.uk/pub/databases/emdb/structures"


def emdb_map(code: str):
    """Primary (post-processed) map for an EMDB entry, e.g. code='12042'."""
    return f"{_EMDB}/EMD-{code}/map/emd_{code}.map.gz"


def emdb_half_map(code: str, n: int):
    """Unfiltered half map ``n`` (1 or 2) for an EMDB entry."""
    return f"{_EMDB}/EMD-{code}/other/emd_{code}_half_map_{n}.map.gz"


def pdbe_fasta(code: str):
    return f"{_PDBE}/api/v2/pdb/entry/{code}/fasta"


def pdbe_mmcif(code: str):
    return f"{_PDBE}/entry-files/download/{code}.cif"


def pdbe_pdb(code: str):
    return f"{_PDBE}/entry-files/download/pdb{code}.ent"


def pdbe_sfcif(code: str):
    return f"{_PDBE}/entry-files/download/r{code}sf.ent"


def rcsb_ligand_cif(code: str):
    return f"{_RCSB}/ligands/download/{code}.cif"


def rcsb_ligand_sdf(code: str):
    return f"{_RCSB}/ligands/download/{code}_ideal.sdf"


def rcsb_mmcif(code: str):
    return f"{_RCSB}/download/{code}.cif"


def rcsb_pdb(code: str):
    return f"{_RCSB}/download/{code}.pdb"


def redo_cif(code: str):
    return f"{_REDO}/{code}/{code}_final.cif"


def redo_mtz(code: str):
    return f"{_REDO}/{code}/{code}_final.mtz"


def uniprot_fasta(entry: str):
    return f"{_UNIPROT}/uniprotkb/{entry}.fasta"
