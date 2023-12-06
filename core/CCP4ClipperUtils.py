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
    """

"""
    Jon Agirre Oct 2015 - created the file and started developing little
                          crystallographic utilities with clipper-python

    A word on conventions
    ---------------------

    In order to keep the operations traceable by i2, non-trivial functions must output:

        * plain text log string (not a file)
        * XML tree in the form of a etree.Element
        * output files (path) if any

    """

import clipper

def read_pdb ( pdbin = "undefined" ) :

    log_string = "\n\t--------- clipper-python snippet: read_pdb ---------\n\n"
    log_string += "\tpdbin: %s\n\n" % pdbin

    from lxml import etree
    xml_root = etree.Element('read_pdb')

    f = clipper.MMDBfile()
    try:
        import ccp4mg
        import mmdb2
        f.SetFlag(mmdb2.MMDBF_IgnoreBlankLines| mmdb2.MMDBF_IgnoreDuplSeqNum | mmdb2.MMDBF_IgnoreNonCoorPDBErrors | mmdb2.MMDBF_IgnoreRemarks)
    except:
        pass
    f.read_file ( pdbin )

    log_string += "\tFile opened and read\n"

    mmol = clipper.MiniMol ()
    f.import_minimol ( mmol )

    log_string += "\tMiniMol imported\n"

    log_string += "\n\t---------------------------------------------------\n\n"

    return log_string, xml_root, mmol



def is_aminoacid ( residue_name = "None" ) :
    residue_names={ 'UNK', 'ALA', 'GLY', 'VAL', 'LEU', 'ILE', 'PRO', 'PHE', 'TYR', 'TRP', 'SER', 'THR', 'CYS', 'MET', 'ASN', 'GLN', 'LYS', 'ARG', 'HIS', 'ASP', 'GLU' }

    if residue_name == "None" :
        return False

    if residue_name in residue_names :
        return True
    else:
        return False



def is_mainchain ( atom_name = "None" ) :

    atom_names = { 'C', 'O', 'CA', 'N', 'C  A', 'O  A', 'CA A', 'N  A' }

    if atom_name in atom_names :
        return True
    else :
        return False
