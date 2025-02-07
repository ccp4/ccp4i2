#!/usr/bin/python3
''' monitor_differences.py is written to analyze differences
between refinement input and output.

Usage: python3 monitor_differences.py file1.mmcif file2.mmcif output.csv

The report is written into CSV file if an output file name is given.
'''

__author__ = "Petr Kolenko, Martin Maly"
__license__ = "Creative Commons Attribution 4.0 International License"
__email__ = "kolenpe1@cvut.cz"
__version__ = "0.5"

from argparse import ArgumentParser
import gemmi
gemmi.set_leak_warnings(False)
import pandas as pd


def makeAddressStr(cra):
    address = cra.chain.name + "/" + cra.residue.name + " "
    address += str(cra.residue.seqid.num)
    if cra.residue.seqid.icode.strip():
        address += str(cra.residue.seqid.icode)
    address += "/"
    address += cra.atom.name
    if cra.atom.has_altloc():
        address += "."
        address += cra.atom.altloc
    return address


def search(st1Cras, st2Cras, output, minCoordDev, minAdpDev):
    records = []
    for cra1 in st1Cras:
        for j, cra2 in enumerate(st2Cras):
            if (
                cra1.atom.name == cra2.atom.name
                and cra1.atom.altloc == cra2.atom.altloc
                and cra1.residue.name == cra2.residue.name
                and cra1.residue.seqid == cra2.residue.seqid
                and cra1.chain.name == cra2.chain.name
            ):
                coordDev = cra1.atom.pos.dist(cra2.atom.pos)
                adpDev = cra2.atom.b_iso - cra1.atom.b_iso
                if coordDev >= minCoordDev or abs(adpDev) >= minAdpDev:
                    record = {
                        "AtomAddress": makeAddressStr(cra1),
                        "CoordDev": round(coordDev, 2),
                        "ADPDev": round(adpDev, 2),
                    }
                    records.append(record)
                del st2Cras[j]
                break
    df = pd.DataFrame.from_records(records)
    if output:
        df.to_csv(output, index=False)
    return df


def main(file1, file2, output=None, minCoordDev=0, minAdpDev=0):
    print("File 1 is " + file1 + ".")
    print("File 2 is " + file2 + ".")
    st1 = gemmi.read_structure(file1)
    st2 = gemmi.read_structure(file2)
    st1Cras = list(st1[0].all())
    st2Cras = list(st2[0].all())
    if len(st1Cras) != len(st2Cras):
        print(
            "Number of atoms in the 1st file does not match the number"
            " of atoms in the second file."
        )
    return search(st1Cras, st2Cras, output, minCoordDev, minAdpDev)


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("file1")
    parser.add_argument("file2")
    parser.add_argument("output")
    parser.add_argument("--minCoordDev", type=float, default=0)
    parser.add_argument("--minAdpDev", type=float, default=0)
    args = parser.parse_args()

    main(args.file1, args.file2, args.output, args.minCoordDev, args.minAdpDev)
