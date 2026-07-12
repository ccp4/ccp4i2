from gemmi import read_ccp4_map, read_mtz_file
from .utils import i2run


def test_8xfm(cif8xfm):
    args = ["density_calculator", "--XYZIN", cif8xfm]
    with i2run(args) as job:
        read_ccp4_map(str(job / "mapout.map"))
        read_mtz_file(str(job / "FPHIOUT.mtz"))
