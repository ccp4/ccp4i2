import xml.etree.ElementTree as ET
import gemmi
from .urls import pdbe_mmcif
from .utils import download, i2run


def _check_output(job):
    max_allowed_rwork = 0.20
    max_allowed_rfree = 0.25
    n_protocols = 10

    gemmi.read_structure(str(job / "XYZOUT.pdb"), format=gemmi.CoorFormat.Pdb)
    for name in ["DIFFPHI", "FPHI"]:
        gemmi.read_mtz_file(str(job / f"{name}OUT.mtz"))

    xml = ET.parse(job / "program.xml")
    root = xml.getroot()

    # Check lorestr protocol structure
    protocols_node = root.find(".//Protocols")
    assert protocols_node is not None, "Missing Protocols section in program.xml"
    protocols = [
        p
        for p in list(protocols_node)
        if p.tag.startswith("P") and p.tag[1:].isdigit()
    ]
    assert len(protocols) == n_protocols, f"Expected {n_protocols} protocols, got {len(protocols)}"

    for protocol in protocols:
        status = protocol.find("Status")
        assert status is not None and status.text == "Finished", "Protocol has not finished successfully"
        rfree = protocol.find("Rfree")
        assert rfree is not None and rfree.text is not None, "Rfree value missing in protocol"
        assert float(rfree.text) > 0.0, "Rfree value invalid in protocol"

    # Best protocol should have acceptable Rfree
    best_protocol_elem = root.find(".//BestProtocol")
    assert best_protocol_elem is not None and best_protocol_elem.text is not None, "Missing BestProtocol"
    best_protocol_num = int(best_protocol_elem.text)
    best_protocol = root.find(f".//Protocols/P{best_protocol_num}")
    assert best_protocol is not None, f"Best protocol P{best_protocol_num} not found"
    best_rwork_elem = best_protocol.find("Rfact")
    assert best_rwork_elem is not None and best_rwork_elem.text is not None, "Missing Rfact in best protocol"
    best_rwork = float(best_rwork_elem.text)
    assert best_rwork < max_allowed_rwork, f"Best Rwork {best_rwork} exceeds threshold {max_allowed_rwork}"
    best_rfree_elem = best_protocol.find("Rfree")
    assert best_rfree_elem is not None and best_rfree_elem.text is not None, "Missing Rfree in best protocol"
    best_rfree = float(best_rfree_elem.text)
    assert best_rfree < max_allowed_rfree, f"Best Rfree {best_rfree} exceeds threshold {max_allowed_rfree}"

    """
    # Check that validation metrics exist in best protocol
    assert best_protocol.find("ramaFav") is not None, "Ramachandran data missing"
    assert best_protocol.find("molprobPercentile") is not None, "Molprobity data missing"
    assert best_protocol.find("clashPercentile") is not None, "Clash percentile data missing"
    """


def test_8xfm(cif8xfm, mtz8xfm):
    with download(pdbe_mmcif("2ac3")) as cif2AC3:
        args = ["lorestr_i2"]
        args += ["--XYZIN", cif8xfm]
        args += ["--F_SIGF", f"fullPath={mtz8xfm}", "columnLabels=/*/*/[FP,SIGFP]"]
        args += ["--FREERFLAG", f"fullPath={mtz8xfm}", "columnLabels=/*/*/[FREE]"]
        args += ["--REFERENCE_MODEL", f"{cif2AC3}"]
        with i2run(args) as job:
            # assert hasLongLigandName(job / "CIFFILE.pdb")
            _check_output(job)