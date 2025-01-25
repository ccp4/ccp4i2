"""
    CCP4ClipperUtils.py: CCP4i2 Project
    Copyright (C) 2015 The University of York

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

    Jon Agirre Oct 2015 - created the file and started developing little
                          crystallographic utilities with clipper-python

    A word on conventions
    ---------------------

    In order to keep the operations traceable by i2, non-trivial functions must output:

        * plain text log string (not a file)
        * XML tree in the form of a etree.Element
        * output files (path) if any
"""


def is_aminoacid(residue_name):
    return residue_name in {
        'UNK', 'ALA', 'GLY', 'VAL', 'LEU', 'ILE', 'PRO',
        'PHE', 'TYR', 'TRP', 'SER', 'THR', 'CYS', 'MET',
        'ASN', 'GLN', 'LYS', 'ARG', 'HIS', 'ASP', 'GLU',
    }
