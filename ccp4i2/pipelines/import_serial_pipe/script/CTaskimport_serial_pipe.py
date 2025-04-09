import os
import sys

from iotbx import mtz, file_reader, reflection_file_reader
import pandas as pd

from ....qtgui import CCP4TaskWidget


class CTaskimport_serial_pipe(CCP4TaskWidget.CTaskWidget):

    TASKNAME = 'import_serial_pipe'
    TASKVERSION = 1.0
    TASKMODULE =['data_entry']
    TASKTITLE = 'Import Serial Pipeline'
    SHORTTASKTITLE = "Import Serial Data"
    DESCRIPTION = 'Import merged data from CrystFEL'
    WHATNEXT = []


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

        self.createLine(["subtitle", "Symmetry (required for CrystFEL) and wavelength"])
        self.openSubFrame(frame=[True])
        self.createLine(['label', 'Load symmetry from ', 'widget', '-guiMode', 'radio', 'SYMMETRY_SOURCE'])
        self.createLine(['tip', 'Reference structure or data', 'widget', 'REFERENCEFILE'], toggle=['SYMMETRY_SOURCE', 'open', 'reference'])
        self.connectDataChanged('REFERENCEFILE', self.updateFromReferenceFile)
        self.createLine(['tip', 'Cell file from CrystFEL', 'widget', 'CELLFILE'], toggle=['SYMMETRY_SOURCE', 'open', 'cellfile'])
        self.connectDataChanged('CELLFILE', self.updateFromCellFile)
        self.createLine(['tip', 'Stream file from CrystFEL', 'widget', 'STREAMFILE'], toggle=['SYMMETRY_SOURCE', 'open', 'streamfile'])
        self.connectDataChanged('STREAMFILE', self.updateFromStreamFile)
        self.createLine(['label', 'Space group', 'widget', 'SPACEGROUP'])
        self.createLine(['label', 'Unit cell', 'widget', 'CELL', 'label', '(6 parameters divided by spaces)'])
        # TO DO: check that were given 6 floats
        self.createLine(['label', 'Wavelength (A)', 'widget', 'WAVELENGTH'])
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
                cell = list(file.crystal_symmetry().unit_cell().parameters())
                for i in range(len(cell)):
                    cell[i] = round(cell[i], 2)
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
            if reflection_file_reader.any_reflection_file(self.container.inputData.HKLIN.fullPath.__str__()).file_type() == 'ccp4_mtz':
                cs, spacegroup, cell_string = self.get_cs_reference(self.container.inputData.HKLIN.fullPath.__str__())
                if spacegroup and cell_string:
                    self.container.inputParameters.SPACEGROUP.set(spacegroup)
                    self.container.inputParameters.CELL.set(cell_string)
            else:
                pass
        else:
            # TO DO - ERROR
            pass


    def updateFromReferenceFile(self):
        def get_wavelength_reference(ref):
            if reflection_file_reader.any_reflection_file(ref).file_type() == 'ccp4_mtz':
                try:
                    mtz_object = mtz.object(file_name=ref)
                    crystal = mtz.crystal(mtz_object=mtz_object, i_crystal=1)
                    dataset = mtz.dataset(mtz_crystal=crystal, i_dataset=0)
                    wavelength = round(float(dataset.wavelength()), 5)
                    print("")
                    print(f"Wavelength found in {ref}:")
                    print(str(wavelength))
                except:
                    sys.stderr.write(
                        f"WARNING: Wavelength could not found in "
                        f"the file {ref}.\n")
                    wavelength = 0
            else:
                wavelength = 0
            return wavelength
        if os.path.isfile(self.container.inputData.REFERENCEFILE.fullPath.__str__()):
            cs, spacegroup, cell_string = self.get_cs_reference(self.container.inputData.REFERENCEFILE.fullPath.__str__())
            if spacegroup and cell_string:
                self.container.inputParameters.SPACEGROUP.set(spacegroup)
                self.container.inputParameters.CELL.set(cell_string)
            wavelength = get_wavelength_reference(self.container.inputData.REFERENCEFILE.fullPath.__str__())
            if wavelength:
                self.container.inputParameters.WAVELENGTH.set(wavelength)
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
        def get_wavelength_streamfile(streamfile):
            energy_eV = None
            wavelength = None
            if os.path.isfile(streamfile + "_tmp"):
                os.remove(streamfile + "_tmp")
            with open(streamfile, "r") as file1:
                with open(streamfile + "_tmp", "a+") as file2:
                    for line in file1:
                        if "photon_energy_eV" in line and len(line.split()) == 3:
                            file2.write(line)
            if os.path.isfile(streamfile + "_tmp"):
                energy_eV_df = pd.read_csv(
                    streamfile + "_tmp", header=None, index_col=False, sep=r'\s+',
                    names=("none1", "none2", "energy_eV"))
                energy_eV_df = energy_eV_df.drop(columns=["none1", "none2"])
                energy_eV_df = energy_eV_df.astype(float)
                energy_eV = energy_eV_df.median()["energy_eV"]
                wavelength = 12398.425 / energy_eV
                wavelength = round(wavelength, 5)
                print("")
                print(f"Wavelength median using file {streamfile}:")
                print(str(wavelength))
                os.remove(streamfile + "_tmp")
            else:
                sys.stderr.write(
                    f"WARNING: Wavelength could not be fitted from "
                    f"the file {streamfile}.\n")
                wavelength = 0
            return wavelength
        if os.path.isfile(self.container.inputData.STREAMFILE.fullPath.__str__()):
            cell, cell_string = get_cell_streamfile(self.container.inputData.STREAMFILE.fullPath.__str__())
            if cell_string:
                self.container.inputParameters.CELL.set(cell_string)
            wavelength = get_wavelength_streamfile(self.container.inputData.STREAMFILE.fullPath.__str__())
            if wavelength:
                self.container.inputParameters.WAVELENGTH.set(wavelength)
        else:
            # TO DO - ERROR
            pass
