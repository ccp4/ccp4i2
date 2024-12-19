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
import csv
import io
from math import sqrt
import gemmi


def calculate_shift(pos1, pos2):
    return sqrt((pos1.x-pos2.x)**2 + (pos1.y-pos2.y)**2 + (pos1.z-pos2.z)**2)


def make_address_str(cra):
    address = cra.chain.name + "/" +  cra.residue.name + " "
    address += str(cra.residue.seqid.num)
    if cra.residue.seqid.icode.strip():
        address += str(cra.residue.seqid.icode)
    address += "/"
    address += cra.atom.name
    if cra.atom.has_altloc():
        address += "."
        address += cra.atom.altloc
    return address


def search(lookup1, lookup2, output, minCoordDev, minADPDev):
    data = []
    for entry1 in lookup1:
        for j, entry2 in enumerate(lookup2):
            if (
                entry1.atom.name == entry2.atom.name
                and entry1.atom.altloc == entry2.atom.altloc
                and entry1.residue.name == entry2.residue.name
                and entry1.residue.seqid.num == entry2.residue.seqid.num
                and entry1.residue.seqid.icode == entry2.residue.seqid.icode
                and entry1.chain.name == entry2.chain.name
            ):
                coordDev = calculate_shift(entry1.atom.pos, entry2.atom.pos)
                ADPDev = entry2.atom.b_iso - entry1.atom.b_iso
                if coordDev >= minCoordDev or abs(ADPDev) >= minADPDev:
                    address = make_address_str(entry1)
                    entry_dict = {
                        "AtomAddress": address,
                        "CoordDev": round(coordDev, 2),
                        "ADPDev": round(ADPDev, 2),
                    }
                    data.append(entry_dict)
                del lookup2[j]
                break
    # JSON
    # import json
    # if output:
    #     with open(output + ".json", "w") as file:
    #         data_json = json.dumps(data, indent=4)
    #         file.writelines(data_json)
    # CSV
    csv_io = io.StringIO()
    fieldnames = ["AtomAddress", "CoordDev", "ADPDev"]
    writer = csv.DictWriter(csv_io, fieldnames=fieldnames, quoting=csv.QUOTE_NONNUMERIC)
    writer.writeheader()
    for entry in data:
        writer.writerow(entry)
    csv_string = csv_io.getvalue()
    if output:
        with open(output, "w") as file:
            file.write(csv_string)
    return csv_string


def main(file1, file2, output=None, minCoordDev=0, minADPDev=0):
    st1 = gemmi.read_structure(file1)
    lookup1 = [entry for entry in st1[0].all()]
    print("File 1 is " + file1 + ".")
    st2 = gemmi.read_structure(file2)
    lookup2 = [entry for entry in st2[0].all()]
    print("File 2 is " + file2 + ".")
    if len(lookup1) != len(lookup2):
        print("Number of atoms in the 1st file does not match the number"
            + " of atoms in the second file.")
    csv_string = search(lookup1, lookup2, output, minCoordDev, minADPDev)
    return csv_string


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("file1")
    parser.add_argument("file2")
    parser.add_argument("output", nargs="?")
    parser.add_argument("--minCoordDev", type=float, default=0)
    parser.add_argument("--minADPDev", type=float, default=0)
    args = parser.parse_args()

    main(args.file1, args.file2, args.output, args.minCoordDev, args.minADPDev)
