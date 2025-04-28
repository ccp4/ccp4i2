import os
import sys

from iotbx import file_reader

from ....qtgui import CCP4TaskWidget


class CTaskimport_serial(CCP4TaskWidget.CTaskWidget):

    TASKNAME = 'import_serial'
    TASKVERSION = 1.0
    TASKTITLE = 'Import Serial'
    SHORTTASKTITLE = "Import Serial Data"
    DESCRIPTION = 'Import merged data from CrystFEL'
    WHATNEXT = []

#################################################################
#
# THIS CTask SHOULD NOT BE USED
# 
# THIS CODE IS NOT LONGER UPDATED
#
# IMPORT_SERIAL SHOULD BE RUN FROM pipelines/import_serial_pipe
#
#################################################################

    def __init__(self,parent):
        CCP4TaskWidget.CTaskWidget.__init__(self,parent)

    #@QtCore.Slot()
    #def updateFromCifFile(self):
    #    #print 'CTaskCif2mtz.updateFromCifFile',self.container.inputData.HKLIN.fileContent
    #    self.container.inputData.SPACEGROUPCELL.cell.set(self.container.inputData.HKLIN.fileContent.cell)
    #    self.container.inputData.SPACEGROUPCELL.spaceGroup.set(self.container.inputData.HKLIN.fileContent.spaceGroup)

    def drawContents(self):
        #self.setProgramHelpFile('cif2mtz')
        # TO DO: help
        #self.createLine(['advice', 'Define Resolution Shells (semi-automated)'])
        self.openFolder(folderFunction='inputData', title='Input Data', followFrom=False)
        self.createLine(['tip', 'Merged data file (I)', 'widget', 'HKLIN'])
        self.connectDataChanged('HKLIN', self.updateFromHklinFile)
        # TO DO: display .hkl filen in the window
        self.createLine(["subtitle", "Half data sets for calculation of statistics"])
        self.openSubFrame(frame=[True])
        self.createLine(['tip', 'Merged data file (I)', 'widget', 'HKLIN1'])
        self.createLine(['tip', 'Merged data file (I)', 'widget', 'HKLIN2'])
        self.createLine(['label', 'Number of resolution bins', 'widget', 'N_BINS'])
        self.closeSubFrame()

        self.createLine(["subtitle", "Symmetry (required for CrystFEL)"])
        self.openSubFrame(frame=[True])
        self.createLine(['label', 'Load symmetry from ', 'widget', '-guiMode', 'radio', 'SYMMETRY_SOURCE'])
        self.createLine(['tip', 'Reference structure or data', 'widget', 'REFERENCEFILE'], toggle=['SYMMETRY_SOURCE', 'open', 'reference'])
        self.connectDataChanged('REFERENCEFILE', self.updateFromReferenceFile)
        self.createLine(['tip', 'Cell file from CrystFEL', 'widget', 'CELLFILE'], toggle=['SYMMETRY_SOURCE', 'open', 'cellfile'])
        self.connectDataChanged('CELLFILE', self.updateFromCellFile)
        self.createLine(['tip', 'Stream file from CrystFEL', 'widget', 'STREAMFILE'], toggle=['SYMMETRY_SOURCE', 'open', 'streamfile'])
        self.connectDataChanged('STREAMFILE', self.updateFromCellFile)
        self.createLine(['label', 'Space group', 'widget', 'SPACEGROUP'])
        self.createLine(['label', 'Unit cell', 'widget', 'CELL', 'label', '(6 parameters divided by spaces)'])
        # TO DO: check that were given 6 floats
        self.closeSubFrame()

        self.createLine(["subtitle", "Resolution"])
        self.openSubFrame(frame=[True])
        self.createLine(['label', 'Low resolution cutoff', 'widget', 'D_MAX'])
        self.createLine(['label', 'High resolution cutoff', 'widget', 'D_MIN'])
        # TO DO: use resolution parameters not only for statistics but really cut the data
        self.closeSubFrame()

    def updateFromCellFile(self):
        def get_cell_cellfile(cellfile):
            with open(cellfile, 'r') as f:
                lines = f.readlines()
            cell = [None, None, None, None, None, None]
            cell_string = None
            for line in lines:
                try:
                    if line.split()[0] == "al" and line.split()[-1] == "deg":
                        cell[3] = float(line.split()[-2])
                    elif line.split()[0] == "be" and line.split()[-1] == "deg":
                        cell[4] = float(line.split()[-2])
                    elif line.split()[0] == "ga" and line.split()[-1] == "deg":
                        cell[5] = float(line.split()[-2])
                    elif line.split()[0] == "a" and line.split()[-1] == "A":
                        cell[0] = float(line.split()[-2])
                    elif line.split()[0] == "b" and line.split()[-1] == "A":
                        cell[1] = float(line.split()[-2])
                    elif line.split()[0] == "c" and line.split()[-1] == "A" and not "centering" in line:
                        cell[2] = float(line.split()[-2])
                except:
                    continue
            print("")
            if (cell[0] and cell[1] and cell[2] and cell[3] and cell[4] and cell[5]):
                cell_string = " ".join(map(str, cell))
                print(f"Unit cell parameters found in file {cellfile}:")
                print(cell_string)
            else:
                sys.stderr.write(
                    f"WARNING: Unit cell parameters could not be parsed from "
                    f"the file {cellfile}.\n"
                    f"Attempt to find the unit cell parameters found "
                    f"in this file: " + " ".join(map(str, cell)) + "\n")
                cell = None
            return cell, cell_string
        if os.path.isfile(self.container.inputData.CELLFILE.fullPath.__str__()):
            cell, cell_string = get_cell_cellfile(self.container.inputData.CELLFILE.fullPath.__str__())
            if cell_string:
                self.container.inputParameters.CELL.set(cell_string)
        else:
            # TO DO - ERROR
            pass


    def get_cs_reference(self, filename):
        file = file_reader.any_file(filename)
        if file.file_type == "hkl" or file.file_type == "pdb":
            try:
                cs = file.crystal_symmetry()
                spacegroup = file.crystal_symmetry().space_group().info()
                cell = file.crystal_symmetry().unit_cell().parameters()
                cell_string = " ".join(map(str, cell))
                print("")
                print(f"Symmetry from the reference file {filename}:")
                print(str(cs))
            except NotImplementedError:
                sys.stderr.write(
                    f"WARNING: Symmetry could not be found in the provided "
                    f"reference file {filename}.\n")
                return None, None, None
            return cs, spacegroup, cell_string


    def updateFromHklinFile(self):
        if os.path.isfile(self.container.inputData.HKLIN.fullPath.__str__()):
            cs, spacegroup, cell_string = self.get_cs_reference(self.container.inputData.HKLIN.fullPath.__str__())
            if spacegroup and cell_string:
                self.container.inputParameters.SPACEGROUP.set(spacegroup)
                self.container.inputParameters.CELL.set(cell_string)
        else:
            # TO DO - ERROR
            pass


    def updateFromReferenceFile(self):
        if os.path.isfile(self.container.inputData.REFERENCEFILE.fullPath.__str__()):
            cs, spacegroup, cell_string = self.get_cs_reference(self.container.inputData.REFERENCEFILE.fullPath.__str__())
            if spacegroup and cell_string:
                self.container.inputParameters.SPACEGROUP.set(spacegroup)
                self.container.inputParameters.CELL.set(cell_string)
        else:
            # TO DO - ERROR
            pass


    def updateFromStreamFile(self):
        def get_cell_streamfile(streamfile):
            cell = [None, None, None, None, None, None]
            cell_string = None
            if os.path.isfile(streamfile + "_tmp"):
                os.remove(streamfile + "_tmp")
            with open(streamfile, "r") as file1:
                with open(streamfile + "_tmp", "a+") as file2:
                    for line in file1:
                        if "Cell parameters " in line and len(line.split()) == 10:
                            file2.write(line)
            if os.path.isfile(streamfile + "_tmp"):
                cell_df = pd.read_csv(
                    streamfile + "_tmp", header=None, index_col=False, sep=r'\s+',
                    names=("none1", "none2", "a", "b", "c", "none3", "alpha", "beta", "gamma", "none4"))
                cell_df = cell_df.drop(columns=["none1", "none2", "none3", "none4"])
                cell_df = cell_df.astype(float)
                cell[0] = round(cell_df.mean()["a"] * 10, 2)
                cell[1] = round(cell_df.mean()["b"] * 10, 2)
                cell[2] = round(cell_df.mean()["c"] * 10, 2)
                cell[3] = round(cell_df.mean()["alpha"], 2)
                cell[4] = round(cell_df.mean()["beta"], 2)
                cell[5] = round(cell_df.mean()["gamma"], 2)
                cell_string = " ".join(map(str, cell))
                print("")
                print(f"Unit cell parameters fit using file {streamfile}:")
                print(cell_string)
                os.remove(streamfile + "_tmp")
            else:
                sys.stderr.write(
                    f"WARNING: Unit cell parameters could not be fitted from "
                    f"the file {streamfile}.\n")
                cell = None
            return cell, cell_string
        if os.path.isfile(self.container.inputData.STREAMFILE.fullPath.__str__()):
            cell, cell_string = get_cs_reference(self.container.inputData.STREAMFILE.fullPath.__str__())
            if cell_string:
                self.container.inputParameters.CELL.set(cell_string)
        else:
            # TO DO - ERROR
            pass
