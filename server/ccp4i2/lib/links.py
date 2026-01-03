_PDBE = "https://www.ebi.ac.uk/pdbe"
_RCSB = "https://files.rcsb.org/download"
_REDO = "https://pdb-redo.eu/db"


def pdbe_eds_map(code: str):
    return f"{_PDBE}/entry-files/{code}.ccp4"


def pdbe_fasta(code: str):
    return f"{_PDBE}/entry/pdb/{code}/fasta"


def pdbe_mmcif(code: str):
    return f"{_PDBE}/entry-files/download/{code}.cif"


def pdbe_mmcif_updated(code: str):
    return f"{_PDBE}/entry-files/download/{code}_updated.cif"


def pdbe_pdb(code: str):
    return f"{_PDBE}/entry-files/download/pdb{code}.ent"


def pdbe_pdb_gz(code: str):
    base = "https://ftp.ebi.ac.uk/pub/databases/rcsb/pdb-remediated/data"
    return f"{base}/structures/divided/pdb/{code[1:3]}/pdb{code}.ent.gz"

def pdbe_sfcif(code: str):
    return f"{_PDBE}/entry-files/download/r{code}sf.ent"


def rcsb_fasta(code: str):
    return f"https://www.rcsb.org/fasta/entry/{code}"


def rcsb_mmcif(code: str):
    return f"{_RCSB}/{code}.cif"


def rcsb_mmcif_gz(code: str):
    return f"{_RCSB}/{code}.cif.gz"


def rcsb_pdb(code: str):
    return f"{_RCSB}/{code}.pdb"


def rcsb_pdb_gz(code: str):
    return f"{_RCSB}/{code}.pdb.gz"


def rcsb_sfcif(code: str):
    return f"{_RCSB}/{code}-sf.cif"


def rcsb_sfcif_gz(code: str):
    return f"{_RCSB}/{code}-sf.cif.gz"


def redo_cif(code: str):
    return f"{_REDO}/{code}/{code}_final.cif"


def redo_mtz(code: str):
    return f"{_REDO}/{code}/{code}_final.mtz"


def redo_pdb(code: str):
    return f"{_REDO}/{code}/{code}_final.pdb"
