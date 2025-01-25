"""
     python/ui/selection_protocols.py: CCP4MG Molecular Graphics Program
     Copyright (C) 2001-2008 University of York, CCLRC
     Copyright (C) 2009 University of York

     This library is free software: you can redistribute it and/or
     modify it under the terms of the GNU Lesser General Public License
     version 3, modified in accordance with the provisions of the 
     license to address the requirements of UK law.
 
     You should have received a copy of the modified GNU Lesser General 
     Public License along with this library.  If not, copies may be 
     downloaded from http://www.ccp4.ac.uk/ccp4license.php
 
     This program is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU Lesser General Public License for more details.
"""

solvent = ["HOH", "H2O", "WAT", "SOL", "TIP", "DOD", "D2O", "DUM"]

dummy = ["DUM"]

solute = [
    "SUL",
    "SO4",
    "NO3",
    "MOH",
    "EOH",
    "GOL",
    "EDO",
    "PO4",
    "CL",
    "ATP",
    "PHO",
]

saccharide = ["GLC", "MAN", "NAG", "RAM", "RIB", "XYS"]

amino_acid = [
    "GLY",
    "ALA",
    "VAL",
    "PRO",
    "SER",
    "THR",
    "LEU",
    "ILE",
    "CYS",
    "ASP",
    "GLU",
    "ASN",
    "GLN",
    "ARG",
    "LYS",
    "MET",
    "MSE",
    "HIS",
    "PHE",
    "TYR",
    "TRP",
    "UNK",
]

nucleic_acid = [
    "DA",
    "DC",
    "DG",
    "DT",
    "Ad",
    "Cd",
    "Gd",
    "Td",
    "A",
    "G",
    "I",
    "U",
    "T",
    "C",
    "ADE",
    "CYT",
    "GUA",
    "INO",
    "THY",
    "URA",
    "YG",
    "PSU",
]

nucleic_acid_code = {
    "GLY": "G",
    "DA": "A",
    "DC": "C",
    "DG": "G",
    "DT": "T",
    "Ad": "A",
    "Cd": "C",
    "Gd": "G",
    "Td": "T",
    "I": "I",
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "U": "U",
    "ADE": "A",
    "CYT": "C",
    "GUA": "G",
    "INO": "I",
    "THY": "T",
    "URA": "U",
    "YG": "W",
    "PSU": "Q",
}

amino_acid_code = {
    "GLY": "G",
    "PRO": "P",
    "ALA": "A",
    "VAL": "V",
    "LEU": "L",
    "ILE": "I",
    "MET": "M",
    "MSE": "M",
    "CYS": "C",
    "PHE": "F",
    "TYR": "Y",
    "TRP": "W",
    "HIS": "H",
    "LYS": "K",
    "GLN": "Q",
    "ASN": "N",
    "ARG": "R",
    "GLU": "E",
    "ASP": "D",
    "SER": "S",
    "THR": "T",
    "YPO": "Y",
    "TPO": "T",
}
