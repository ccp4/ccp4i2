import gemmi
from pathlib import Path


def detect_file_type(file_path):
    path = Path(file_path)
    if not path.exists():
        return "File does not exist"

    # Try detecting an MTZ file
    try:
        mtz = gemmi.read_mtz_file(str(path))
        if mtz:  # If it successfully reads, it's an MTZ file
            return "MTZ file"
    except Exception:
        pass

    # Try detecting an mmCIF file and its subtype
    try:
        cif_doc = gemmi.cif.read(str(path))
        if cif_doc:  # If it successfully reads, it's an mmCIF file
            # Check for coordinate CIF (structure)
            for block in cif_doc:
                if block.find_mmcif_category("_atom_site"):
                    return "mmCIF coordinate file"
                if block.find_mmcif_category("_refln"):
                    return "mmCIF reflection file"
                if block.find_mmcif_category("_chem_comp"):
                    return "mmCIF ligand/geometry file"
            return "mmCIF file (unknown subtype)"
    except Exception:
        pass

    # Try detecting a PDB file
    try:
        pdb_structure = gemmi.read_structure(str(path))
        if pdb_structure:  # If it successfully reads, it's a PDB file
            return "PDB file"
    except Exception:
        pass

    # Try detecting a sequence file (FASTA, PIR, etc.) using BioPython
    # Import lazily to avoid numpy dependency at module load time
    try:
        from Bio import SeqIO
        with path.open("r") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                return "FASTA file"
    except Exception:
        pass

    try:
        from Bio import SeqIO
        with path.open("r") as handle:
            for record in SeqIO.parse(handle, "pir"):
                return "PIR file"
    except Exception:
        pass

    return "Unknown file type"


if __name__ == "__main__":
    file_path = input("Enter the file path: ")
    file_type = detect_file_type(file_path)
    print(f"Detected file type: {file_type}")
