"""Add to the MDM2 project (scenario_refine.py, scenario_tools.py) the runs
the search-model preparation pages are illustrated from: MDM2 treated as
unknown, solved with its homologue MDMX.

PDB entry 3dab holds four copies of the MDMX N-terminal domain (57%
identical to MDM2) with p53 peptides. ClustalW aligns the MDM2 and MDMX
sequences; Chainsaw and Sculptor prune MDMX to the alignment; the Phaser
ensembler superposes an MDM2 and an MDMX structure into one ensemble; and
Phaser places the Chainsaw model in the MDM2 data.

    env CCP4I2_HOME=/tmp/docs-home ccp4-python ../docs/user/tools/scenario_prep.py

Run scenario_refine.py and scenario_tools.py first. Each preparation task's
last job gets an unrun clone.
"""
import gemmi

from scenario_common import clone_last, fetch, i2run as _i2run, inputs_dir
from scenario_refine import DEMO, PROJECT


def i2run(*args):
    _i2run(PROJECT, *args)


def main():
    mdmx = fetch("https://www.ebi.ac.uk/pdbe/entry-files/download/pdb3dab.ent",
                 "3dab.pdb")
    chain_a = gemmi.read_structure(str(mdmx))[0]["A"].get_polymer()
    fasta = inputs_dir() / "mdmx.fasta"
    fasta.write_text(">MDMX 3dab chain A\n"
                     + gemmi.one_letter_code([r.name for r in chain_a]) + "\n")

    i2run("clustalw", "--SEQUENCELISTORALIGNMENT", "SEQUENCELIST",
          "--SEQIN", f"fullPath={DEMO / '4hg7.seq'}",
          "--SEQIN", f"fullPath={fasta}")
    alignment = "fileOut=clustalw[-1].ALIGNMENTOUT"
    i2run("chainsaw", "--XYZIN", f"fullPath={mdmx}", "selection/text=A/",
          "--ALIGNIN", alignment)
    i2run("sculptor", "--XYZIN", f"fullPath={mdmx}", "selection/text=A/",
          "--ALIGNMENTORSEQUENCEIN", "ALIGNMENT", "--ALIGNIN", alignment)
    i2run("phaser_ensembler",
          "--XYZIN_LIST", f"fullPath={DEMO / '4qo4.pdb'}",
          "--XYZIN_LIST", f"fullPath={mdmx}")

    # The point of it: the pruned homologue solves the structure.
    i2run("phaser_simple_phil",
          "--F_SIGF", "fileOut=aimless_pipe[-1].HKLOUT[0]",
          "--FREERFLAG", "fileOut=aimless_pipe[-1].FREEROUT",
          "--COMP_BY", "ASU",
          "--ASUFILE", "fileOut=ProvideAsuContents[-1].ASUCONTENTFILE",
          "--XYZIN", "fileOut=chainsaw[-1].XYZOUT",
          "--SEARCHSEQUENCEIDENTITY", "0.57")

    for task in ("clustalw", "chainsaw", "sculptor", "phaser_ensembler"):
        clone_last(PROJECT, task)


if __name__ == "__main__":
    main()
