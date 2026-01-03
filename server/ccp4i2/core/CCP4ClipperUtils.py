"""
    Copyright (C) 2015 The University of York

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
