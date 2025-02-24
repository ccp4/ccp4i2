"""
    density_calculator.py: CCP4 GUI Project

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

import os
import gemmi
from core.CCP4ErrorHandling import SEVERITY_WARNING
from core.CCP4PluginScript import CPluginScript
from core.CCP4XtalData import CMapCoeffsDataFile


class density_calculator(CPluginScript):
    TASKMODULE = "expt_data_utility"
    TASKTITLE = "Density Calculator"
    TASKNAME = "density_calculator"
    TASKVERSION = 0.1
    MAINTAINER = "paul.bond@york.ac.uk"
    RUNEXTERNALPROCESS = False
    WHATNEXT = ["coot_rebuild"]

    def startProcess(self, command=None, handler=None, **kw):
        xyzin = os.path.join(self.getWorkDirectory(), "xyzin.xyz")
        hklout = os.path.join(self.getWorkDirectory(), "hklout.mtz")
        mapout = os.path.join(self.getWorkDirectory(), "mapout.map")

        # Write selected atoms and read the structure
        self.container.inputData.XYZIN.getSelectedAtomsPdbFile(xyzin)
        structure = gemmi.read_structure(xyzin, format=gemmi.CoorFormat.Detect)

        # Calculate the density map
        params = self.container.controlParameters
        if params.FORM_FACTOR == "xray":
            dencalc = gemmi.DensityCalculatorX()
        elif params.FORM_FACTOR == "electron":
            dencalc = gemmi.DensityCalculatorE()
        elif params.FORM_FACTOR == "neutron":
            dencalc = gemmi.DensityCalculatorN()
        dencalc.d_min = params.D_MIN
        dencalc.rate = params.RATE
        if params.BLUR_MODE == "refmac":
            dencalc.set_refmac_compatible_blur(structure[0])
        elif params.BLUR_MODE == "custom":
            dencalc.blur = params.BLUR
        dencalc.cutoff = params.CUTOFF
        use_mott_bethe = params.FORM_FACTOR == "xray" and params.MOTT_BETHE
        if use_mott_bethe:
            dencalc.addends.subtract_z()
        dencalc.set_grid_cell_and_spacegroup(structure)
        dencalc.put_model_density_on_grid(structure[0])

        # Write CCP4 map file
        ccp4 = gemmi.Ccp4Map()
        ccp4.grid = dencalc.grid
        ccp4.update_ccp4_header()
        ccp4.write_ccp4_map(mapout)

        # Write map coeffients to MTZ file
        sf_grid = gemmi.transform_map_to_f_phi(dencalc.grid)
        unblur = dencalc.blur if params.UNBLUR else 0
        asu_data = sf_grid.prepare_asu_data(
            dmin=dencalc.d_min, mott_bethe=use_mott_bethe, unblur=unblur
        )
        mtz = gemmi.Mtz(with_base=True)
        mtz.spacegroup = sf_grid.spacegroup
        mtz.set_cell_for_all(sf_grid.unit_cell)
        mtz.add_dataset("Calculated")
        mtz.add_column("FC", "F")
        mtz.add_column("PHIC", "P")
        mtz.set_data(asu_data)
        mtz.write_to_file(hklout)

        # Set output objects
        annotation_map = "Calculated map"
        annotation_hkl = "Calculated map coefficients"
        if params.BLUR_MODE != "none":
            annotation_map += " (blurred)"
            if not params.UNBLUR:
                annotation_hkl += " (blurred)"
        output_data = self.container.outputData
        output_data.MAPOUT.setFullPath(mapout)
        output_data.MAPOUT.annotation.set(annotation_map)
        output_data.FPHIOUT.annotation.set(annotation_hkl)
        output_data.FPHIOUT.subType.set(CMapCoeffsDataFile.SUBTYPE_NORMAL)
        output_data.FPHIOUT.contentFlag.set(CMapCoeffsDataFile.CONTENT_FLAG_FPHI)
        files = ["FPHIOUT"]
        columns = ["FC,PHIC"]
        error = self.splitHklout(files, columns, hklout)
        if error.maxSeverity() > SEVERITY_WARNING:
            return CPluginScript.FAILED

        return CPluginScript.SUCCEEDED
