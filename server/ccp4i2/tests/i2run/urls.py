# Copyright (C) 2025-2026 Newcastle University
# Copyright (C) 2025-2026 University of York
#
# This file is part of CCP4i2.
#
# CCP4i2 is free software: you can redistribute it and/or modify it
# under the terms of the GNU Lesser General Public License version 3,
# modified in accordance with the provisions of the license to address
# the requirements of UK law.
#
# See https://www.ccp4.ac.uk/ccp4license.php for details.
_PDBE = "https://www.ebi.ac.uk/pdbe"
_RCSB = "https://files.rcsb.org"
_REDO = "https://pdb-redo.eu/db"
_UNIPROT = "https://rest.uniprot.org"


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
