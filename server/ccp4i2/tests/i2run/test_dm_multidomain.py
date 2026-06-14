import xml.etree.ElementTree as ET

import gemmi

from .utils import demoData, i2run


def test_dm_multidomain():
    """Multi-domain NCS averaging on the AHIR driver case (6 NCS copies, T sym).

    The DARPIN / NADP-knob / helical-shell domains follow different NCS
    transformations, so each gets its own operator set and mask. Observed
    intensities are French-Wilsoned to F by makeHklin0; phases come from the
    model (no experimental phases for this case).
    """
    args = ["dm_multidomain"]
    args += ["--F_SIGF", demoData("ahir", "ahir_observed.mtz")]
    args += ["--ABCD", demoData("ahir", "ahir_phases.mtz")]
    args += ["--FREERFLAG", demoData("ahir", "ahir_freer.mtz")]
    args += ["--XYZIN", demoData("ahir", "ahir_model.cif")]
    args += ["--DOMAINS", "chainId=A", "firstRes=340", "lastRes=485", "mode=average"]
    args += ["--DOMAINS", "chainId=A", "firstRes=140", "lastRes=339", "mode=refine"]
    args += ["--DOMAINS", "chainId=A", "firstRes=13", "lastRes=139", "mode=exclude"]
    with i2run(args) as job:
        for name in ["ABCDOUT", "FPHIOUT"]:
            gemmi.read_mtz_file(str(job / f"{name}.mtz"))
        ET.parse(job / "program.xml")
