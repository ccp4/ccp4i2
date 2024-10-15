#!/usr/bin/python3
''' monitor_refinement_differences.py is written to analyze differences
between refinement input and output.

Usage: python3 monitor_refinement_differences.py INPUT OUTPUT

The report is written into the JSON file (default output.json).
'''

__author__ = "Petr Kolenko"
__license__ = "Creative Commons Attribution 4.0 International License"
__email__ = "kolenpe1@cvut.cz"
__version__ = "0.3"

#https://patorjk.com/software/taag/#p=display&f=Stick%20Letters&t=%20write_smthg



import sys
import os
import argparse
import textwrap
import re
from datetime import datetime
from math import sqrt


# Checking for file types
def is_cif_file(file_name):
    """Check if the given file has a .cif or .mmcif extension."""
    return file_name.lower().endswith('.cif') or file_name.lower().endswith('.mmcif')

def is_pdb_file(file_name):
    """Check if the given file has a .pdb extension."""
    return file_name.lower().endswith('.pdb')



#    __   ___       __       __   __   __  
#   |__) |__   /\  |  \     |__) |  \ |__) 
#   |  \ |___ /~~\ |__/ ___ |    |__/ |__) 

def read_pdb(filename):
    with open(filename, "r") as file:
        lines = file.readlines()

    # initiation of lists
    atomSiteId = []
    atomSiteLabelAtomId = []
    atomSiteLabelAltId = []
    atomSiteLabelCompId = []
    atomSitePDBInsCode = []
    atomSiteAuthAsymId = []
    atomSiteAuthSeqId = []
    atomSiteCartnX = []
    atomSiteCartnY = []
    atomSiteCartnZ = []
    atomSiteOccupancy = []
    atomSiteBIsoOrEquiv = []

    for i in range(len(lines)):
        if lines[i][0:4] == "ATOM" or lines[i][0:6] == "HETATM":
            atomSiteId.append(lines[i][6:11])
            atomSiteLabelAtomId.append(lines[i][13:16])
            atomSiteLabelAltId.append(lines[i][16:17])
            atomSiteLabelCompId.append(lines[i][17:20])
            atomSiteAuthAsymId.append(lines[i][21:22])
            atomSiteAuthSeqId.append(lines[i][22:26])
            atomSitePDBInsCode.append(lines[i][26:27])
            atomSiteCartnX.append(float(lines[i][30:38]))
            atomSiteCartnY.append(float(lines[i][38:46]))
            atomSiteCartnZ.append(float(lines[i][46:54]))
            atomSiteOccupancy.append(float(lines[i][54:60]))
            atomSiteBIsoOrEquiv.append(float(lines[i][60:66]))

    return atomSiteId, atomSiteLabelAtomId, atomSiteLabelAltId, \
        atomSiteLabelCompId,atomSitePDBInsCode, atomSiteAuthAsymId, \
        atomSiteAuthSeqId, atomSiteCartnX, atomSiteCartnY, atomSiteCartnZ, \
        atomSiteOccupancy, atomSiteBIsoOrEquiv



#    __   ___       __       __     ___ 
#   |__) |__   /\  |  \     /  ` | |__  
#   |  \ |___ /~~\ |__/ ___ \__, | |    
                                       
# Reading cif files
def read_cif(filename):

    # simple reading lines
    with open(filename, "r") as file:
        lines = file.readlines()

    # initiation of lists
    counter = -1
    atomSiteId = []             # cislo atomu
    atomSiteLabelAtomId = []    # znacka atomu, napr. CA pro Calfa
    atomSiteLabelAltId = []     # alternativa
    atomSiteLabelCompId = []    # tripismenovy kod pro reziduum
    atomSitePDBInsCode = []
    atomSiteAuthAsymId = []     # asi retezec?
    atomSiteAuthSeqId = []      # poradi v sekvenci
    atomSiteCartnX = []         # souradnice X
    atomSiteCartnY = []         # souradnice Y
    atomSiteCartnZ = []         # souradnice Z
    atomSiteOccupancy = []      # okupance
    atomSiteBIsoOrEquiv = []    # ADP
    
    # analysis of positions within the CIF file
    for i in range(len(lines)):
        line = lines[i].split()
        if lines[i][0:4] == "loop":
            counter = -1
        if "atom_site.id" in lines[i] and len(line) == 1:
            positionSiteId = counter
        if "atom_site.label_atom_id" in lines[i] and len(line) ==1:
            positionSiteLabelAtomId = counter
        if "atom_site.label_alt_id" in lines[i] and len(line) ==1:
            positionSiteLabelAltId = counter
        if "atom_site.label_comp_id" in lines[i] and len(line) ==1:
            positionSiteLabelCompId = counter
        if "atom_site.pdbx_PDB_ins_code" in lines[i] and len(line) ==1:
            positionPDB_ins_code = counter
        if "atom_site.auth_asym_id" in lines[i] and len(line) ==1:
            positionSiteAuthAsymId = counter
        if "atom_site.auth_seq_id" in lines[i] and len(line) ==1:
            positionSiteAuthSeqId = counter
        if "atom_site.Cartn_x" in lines[i] and len(line) ==1:
            positionSiteCartnX = counter
        if "atom_site.Cartn_y" in lines[i] and len(line) ==1:
            positionSiteCartnY = counter
        if "atom_site.Cartn_z" in lines[i] and len(line) ==1:
            positionSiteCartnZ = counter
        if "atom_site.occupancy" in lines[i] and len(line) ==1:
            positionSiteOccupancy = counter
        if "atom_site.B_iso_or_equiv" in lines[i] and len(line) ==1:
            positionSiteBIsoOrEquiv = counter
        counter = counter + 1

        
    # Filling of lists with values
    for i in range(len(lines)):
        # looking for loop
        if "loop" in lines[i] or lines[i] == []:
            scan = False
        # looking for loop and X coordinates
        if "atom_site.Cartn_x" in lines[i]:
            scan = True


        line = lines[i].split()
        # generation of arrays
        if len(line) > positionSiteBIsoOrEquiv and scan:
            atomSiteId.append(line[positionSiteId])
            atomSiteLabelAtomId.append(line[positionSiteLabelAtomId])
            atomSiteLabelAltId.append(line[positionSiteLabelAltId])
            atomSiteLabelCompId.append(line[positionSiteLabelCompId])
            atomSitePDBInsCode.append(line[positionPDB_ins_code])
            atomSiteAuthAsymId.append(line[positionSiteAuthAsymId])
            atomSiteAuthSeqId.append(line[positionSiteAuthSeqId])
            atomSiteCartnX.append(float(line[positionSiteCartnX]))
            atomSiteCartnY.append(float(line[positionSiteCartnY]))
            atomSiteCartnZ.append(float(line[positionSiteCartnZ]))
            atomSiteOccupancy.append(float(line[positionSiteOccupancy]))
            atomSiteBIsoOrEquiv.append(float(line[positionSiteBIsoOrEquiv]))

    return atomSiteId, atomSiteLabelAtomId, atomSiteLabelAltId, \
        atomSiteLabelCompId,atomSitePDBInsCode, atomSiteAuthAsymId, \
        atomSiteAuthSeqId, atomSiteCartnX, atomSiteCartnY, atomSiteCartnZ, \
        atomSiteOccupancy, atomSiteBIsoOrEquiv



def calculate_shift(a, b, c, d, e, f):
    return sqrt((a-b)**2 + (c-d)**2 + (e-f)**2)



#         __    ___  ___           __   __       
#   |  | |__) |  |  |__         | /__` /  \ |\ | 
#   |/\| |  \ |  |  |___ ___ \__/ .__/ \__/ | \| 
                                                

def write_json(atomSiteLabelAtomId, atomSiteLabelCompId, atomSiteAuthSeqId, \
        atomSitePDBInsCode, atomSiteLabelAltId, atomSiteAuthAsymId, \
        atomSiteCartnX, atomSiteCartnY, atomSiteCartnZ, \
        atomSiteCartnX2, atomSiteCartnY2, atomSiteCartnZ2, \
        atomSiteBIsoOrEquiv, atomSiteBIsoOrEquiv2, \
        output):
    data_json = []
    data_json.append("[")
    for i in range(len(atomSiteCartnX)):
        data_json.append("\n    {\n")
        data_json.append(8 * " " + "\"Atom_Name\": \"" 
            + atomSiteLabelAtomId[i].strip() + "\",\n")
        data_json.append(8 * " " + "\"Residue_Name\": \"" 
            + atomSiteLabelCompId[i].strip() + "\",\n")
        if (atomSitePDBInsCode[i] != " ") and (atomSitePDBInsCode[i] != ".") \
            and (atomSitePDBInsCode[i] != "?"):
            data_json.append(8 * " " + "\"Residue_Number\": \"" 
                + atomSiteAuthSeqId[i].strip() + "." 
                + atomSitePDBInsCode[i] + "\",\n")
        else:
            data_json.append(8 * " " + "\"Residue_Number\": \"" 
                + atomSiteAuthSeqId[i].strip() + "\",\n")
        data_json.append(8 * " " + "\"Alternative\": \"" 
            + atomSiteLabelAltId[i].strip() + "\",\n")
        data_json.append(8 * " " + "\"Chain_Id\": \"" 
            + atomSiteAuthAsymId[i].strip() + "\",\n")
        data_json.append(8 * " " + "\"Coordinate_Deviation\": " 
            + str(calculate_shift(atomSiteCartnX[i], atomSiteCartnX2[i], \
                atomSiteCartnY[i], atomSiteCartnY2[i], atomSiteCartnZ[i],\
                atomSiteCartnZ2[i])) + ",\n")
        data_json.append(8 * " " + "\"ADP_Deviation\": " 
            + str(atomSiteBIsoOrEquiv2[i] - atomSiteBIsoOrEquiv[i]) + "\n")
        data_json.append("    }")
        data_json.append(",")
    data_json.pop()
    data_json.append("\n]")

    # writing the log file
    with open(output, "w") as file:
        file.writelines(data_json)

    # printing the log file
    f = open(output, 'r')
    content = f.read()
    print(content)
    f.close()



#    __   ___       __   __                 __    ___  ___           __   __       
#   /__` |__   /\  |__) /  ` |__|     |  | |__) |  |  |__         | /__` /  \ |\ | 
#   .__/ |___ /~~\ |  \ \__, |  | ___ |/\| |  \ |  |  |___ ___ \__/ .__/ \__/ | \| 
                                                                                  
def search_write_json(atomSiteLabelAtomId, atomSiteLabelAtomId2, \
        atomSiteLabelCompId, atomSiteLabelCompId2, \
        atomSiteAuthSeqId, atomSiteAuthSeqId2,\
        atomSitePDBInsCode, atomSitePDBInsCode2, \
        atomSiteLabelAltId, atomSiteLabelAltId2, \
        atomSiteAuthAsymId, atomSiteAuthAsymId2, \
        atomSiteCartnX, atomSiteCartnY, atomSiteCartnZ, \
        atomSiteCartnX2, atomSiteCartnY2, atomSiteCartnZ2, \
        atomSiteBIsoOrEquiv, atomSiteBIsoOrEquiv2, \
        output, minCoordDev, minADPDev):
    data_json = []
    data_json.append("[")

    for i in range(len(atomSiteLabelAtomId)):
        if atomSiteLabelAltId[i].strip() == "." or atomSiteLabelAltId[i].strip() == "?":
            atomSiteLabelAltId[i] = ""
        if atomSitePDBInsCode[i].strip() == "." or atomSitePDBInsCode[i].strip() == "?":
            atomSitePDBInsCode[i] = ""
        for j in range(len(atomSiteLabelAtomId2)):
            if atomSiteLabelAltId2[j].strip() == "." or atomSiteLabelAltId2[j].strip() == "?":
                atomSiteLabelAltId2[j] = ""
            if atomSitePDBInsCode2[j].strip() == "." or atomSitePDBInsCode2[j].strip() == "?":
                atomSitePDBInsCode2[j] = ""
            if (atomSiteLabelAtomId[i].strip() == atomSiteLabelAtomId2[j].strip()) \
                and (atomSiteLabelAltId[i].strip() == atomSiteLabelAltId2[j].strip()) \
                and (atomSiteLabelCompId[i].strip() == atomSiteLabelCompId2[j].strip()) \
                and (atomSiteAuthAsymId[i].strip() == atomSiteAuthAsymId2[j].strip()) \
                and (atomSiteAuthSeqId[i].strip() == atomSiteAuthSeqId2[j].strip())\
                and (atomSitePDBInsCode[i].strip() == atomSitePDBInsCode2[j].strip()):
                    coordDev = calculate_shift(atomSiteCartnX[i], \
                        atomSiteCartnX2[j], atomSiteCartnY[i], atomSiteCartnY2[j], \
                        atomSiteCartnZ[i], atomSiteCartnZ2[j])
                    ADPDev = atomSiteBIsoOrEquiv2[j] - atomSiteBIsoOrEquiv[i]
                    if coordDev >= minCoordDev or abs(ADPDev) >= minADPDev:
                        # writing procedure
                        data_json.append("\n    {\n")
                        # data_json.append(8 * " " + "\"Atom_Name\": \"" 
                        #     + atomSiteLabelAtomId[i].strip() + "\",\n")
                        # data_json.append(8 * " " + "\"Residue_Name\": \"" 
                        #     + atomSiteLabelCompId[i].strip() + "\",\n")
                        # if (atomSitePDBInsCode[i] != " ") \
                        #     and (atomSitePDBInsCode[i] != ".") \
                        #     and (atomSitePDBInsCode[i] != "?"):
                        #     data_json.append(8 * " " + "\"Residue_Number\": \"" 
                        #         + atomSiteAuthSeqId[i].strip() + "." 
                        #         + atomSitePDBInsCode[i].strip() + "\",\n")
                        # else:
                        #     data_json.append(8 * " " + "\"Residue_Number\": \"" 
                        #         + atomSiteAuthSeqId[i].strip() + "\",\n")
                        # data_json.append(8 * " " + "\"Alternative\": \"" 
                        #     + atomSiteLabelAltId[i].strip() + "\",\n")
                        # data_json.append(8 * " " + "\"Chain_Id\": \"" 
                        #     + atomSiteAuthAsymId[i].strip() + "\",\n")
                        address = atomSiteAuthAsymId[i].strip()
                        address += "/"
                        address += atomSiteLabelCompId[i].strip()
                        address += " "
                        address += atomSiteAuthSeqId[i].strip()
                        if atomSitePDBInsCode[i].strip():
                            address += atomSitePDBInsCode[i].strip()
                        address += "/"
                        address += atomSiteLabelAtomId[i].strip()
                        if atomSiteLabelAltId[i].strip():
                            address += "."
                            address += atomSiteLabelAltId[i].strip()
                        data_json.append(8 * " " + "\"AtomAddress\": \"" + address + "\",\n")
                        data_json.append(8 * " " + "\"CoordDev\": " + str(round(coordDev, 2)) + ",\n")
                        data_json.append(8 * " " + "\"ADPDev\": " + str(round(ADPDev, 2)) + "\n")
                        data_json.append("    }")
                        data_json.append(",")
                    break
    data_json.pop()
    data_json.append("\n]")

    # writing the log file
    with open(output, "w") as file:
        file.writelines(data_json)

    # printing the log file
    f = open(output, 'r')
    content = f.read()
    # print(content)
    f.close()








def main(argv):
    # check for command line integrity
    if len(argv) < 2 or len(argv) >= 6:
        print("Usage: python3 monitor_refinement_differences.py" 
            +" <file1> <file2> [output - optional] [minimal coordination deviation] [minimal ADP deviation]")
        sys.exit(1)
    
    # Assigns the arguments
    file1 = argv[0]
    file2 = argv[1]

    # Check if the optional argument is provided
    if len(argv) >= 3:
        output = argv[2]
    else:
        output = "report.json"

    if len(argv) >= 4:
        try:
            minCoordDev = float(argv[3])
        except:
            minCoordDev = 0
    else:
        minCoordDev = 0
    if len(argv) >= 5:
        try:
            minADPDev = float(argv[4])
        except:
            minADPDev = 0
    else:
        minADPDev = 0
        

    # Checking for file formats and reading values
    if is_cif_file(file1):
        atomSiteId, \
        atomSiteLabelAtomId, \
        atomSiteLabelAltId, \
        atomSiteLabelCompId, \
        atomSitePDBInsCode, \
        atomSiteAuthAsymId, \
        atomSiteAuthSeqId, \
        atomSiteCartnX, \
        atomSiteCartnY, \
        atomSiteCartnZ, \
        atomSiteOccupancy, \
        atomSiteBIsoOrEquiv = read_cif(file1)
        print("File 1 is " + file1 + ".")
    elif is_pdb_file(file1):
        atomSiteId, \
        atomSiteLabelAtomId, \
        atomSiteLabelAltId, \
        atomSiteLabelCompId, \
        atomSitePDBInsCode, \
        atomSiteAuthAsymId, \
        atomSiteAuthSeqId, \
        atomSiteCartnX, \
        atomSiteCartnY, \
        atomSiteCartnZ, \
        atomSiteOccupancy, \
        atomSiteBIsoOrEquiv = read_pdb(file1)
        print("File 1 is " + file1 + ".")
    else:
        print("File 1 does not have the right format.")
        sys.exit(1)
    if is_cif_file(file2):
        atomSiteId2, \
        atomSiteLabelAtomId2, \
        atomSiteLabelAltId2, \
        atomSiteLabelCompId2, \
        atomSitePDBInsCode2, \
        atomSiteAuthAsymId2, \
        atomSiteAuthSeqId2, \
        atomSiteCartnX2, \
        atomSiteCartnY2, \
        atomSiteCartnZ2, \
        atomSiteOccupancy2, \
        atomSiteBIsoOrEquiv2 = read_cif(file2)
        print("File 2 is " + file2 + ".")
    elif is_pdb_file(file2):
        atomSiteId2, \
        atomSiteLabelAtomId2, \
        atomSiteLabelAltId2, \
        atomSiteLabelCompId2, \
        atomSitePDBInsCode2, \
        atomSiteAuthAsymId2, \
        atomSiteAuthSeqId2, \
        atomSiteCartnX2, \
        atomSiteCartnY2, \
        atomSiteCartnZ2, \
        atomSiteOccupancy2, \
        atomSiteBIsoOrEquiv2 = read_pdb(file2)
        print("File 2 is " + file2 + ".")
    else:
        print("File 2 does not have the right format.")
        sys.exit(1)

    # Writing the JSON file
    if len(atomSiteLabelAtomId) == len(atomSiteLabelAtomId2):
        pass
        #write_json(atomSiteLabelAtomId, atomSiteLabelCompId, atomSiteAuthSeqId,\
        #    atomSitePDBInsCode, atomSiteLabelAltId, atomSiteAuthAsymId, \
        #    atomSiteCartnX, atomSiteCartnY, atomSiteCartnZ, \
        #    atomSiteCartnX2, atomSiteCartnY2, atomSiteCartnZ2, \
        #    atomSiteBIsoOrEquiv, atomSiteBIsoOrEquiv2, \
        #    output)
    else:
        print("Number of atoms in the 1st file does not match the number"
            + " of atoms in the second file.")
    search_write_json(atomSiteLabelAtomId, atomSiteLabelAtomId2, \
        atomSiteLabelCompId, atomSiteLabelCompId2, \
        atomSiteAuthSeqId, atomSiteAuthSeqId2,\
        atomSitePDBInsCode, atomSitePDBInsCode2,\
        atomSiteLabelAltId, atomSiteLabelAltId2,\
        atomSiteAuthAsymId, atomSiteAuthAsymId2,\
        atomSiteCartnX, atomSiteCartnY, atomSiteCartnZ, \
        atomSiteCartnX2, atomSiteCartnY2, atomSiteCartnZ2, \
        atomSiteBIsoOrEquiv, atomSiteBIsoOrEquiv2, \
        output, minCoordDev, minADPDev)
            


if __name__ == "__main__":
    main(sys.argv[1:])

