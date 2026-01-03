_PDBE = "https://www.ebi.ac.uk/pdbe"
_RCSB = "https://files.rcsb.org"
_REDO = "https://pdb-redo.eu/db"


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
