import gemmi
from .utils import demoData, i2run


def test_gamma_model():
    """Test adding fractional coordinates to gamma model."""
    args = ["add_fractional_coords"]
    args += ["--XYZIN", demoData("gamma", "gamma_model.pdb")]
    with i2run(args) as job:
        # The wrapper writes xyzout.cif. (This test looked for XYZOUT.cif,
        # which only a case-insensitive filesystem finds.)
        xyzout = job / "xyzout.cif"
        assert xyzout.exists(), f"No xyzout.cif: {list(job.iterdir())}"
        st = gemmi.read_structure(str(xyzout))
        assert len(st[0]) > 0, "Output structure has no chains"
        # From a PDB file, residues must still get their label_seq_id (they
        # came out as '.'), and every atom its fractional coordinates.
        atoms = gemmi.cif.read(str(xyzout))[0].get_mmcif_category("_atom_site.")
        assert all(key in atoms for key in ("fract_x", "fract_y", "fract_z"))
        polymer = [s for s, g in zip(atoms["label_seq_id"], atoms["group_PDB"])
                   if g == "ATOM"]
        assert polymer and all(s not in (None, False, ".") for s in polymer)
