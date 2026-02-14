import clipper

from ccp4i2.core.CCP4PluginScript import CPluginScript


class SubtractNative(CPluginScript):
    TASKNAME = "SubtractNative"

    def startProcess(self):
        mtz_file = clipper.CCP4MTZfile()
        hkl_info = clipper.HKL_info()
        coefficientsPath = str(self.container.inputData.MAPIN.fullPath)
        mtz_file.open_read(coefficientsPath)
        mtz_file.import_hkl_info(hkl_info)
        sg = hkl_info.spacegroup()
        cell = hkl_info.cell()
        fphidata = clipper.HKL_data_F_phi_float(hkl_info)
        mtz_file.import_hkl_data(fphidata, str("/*/*/[F,PHI]"))
        mtz_file.close_read()
        # Clipper will sample the output map according to Fourier theory and hte nominal resolution
        # for visualisation, it is generally nicer to make things a bit more finely sampled
        fudgedResolution = hkl_info.resolution()
        fudgedResolution.init((2.0 / 3.0) * hkl_info.resolution().limit())
        mygrid = clipper.Grid_sampling(sg, cell, fudgedResolution)
        mymap = clipper.Xmap_float(sg, cell, mygrid)
        mymap.fft_from(fphidata)
        f = clipper.MMDBfile()
        f.read_file(str(self.container.inputData.XYZIN.fullPath))
        mmol = clipper.MiniMol()
        f.import_minimol(mmol)
        atoms = mmol.atom_list()
        for atom in atoms:
            fraction = float(self.container.controlParameters.FRACTION)
            modifiedOccupancy = atom.occupancy() * fraction
            try:
                atom.set_occupancy(modifiedOccupancy)
            except Exception as err:
                print(err)

        xmap = clipper.Xmap_float(sg, cell, mygrid)
        try:
            EDCalculator = clipper.EDcalc_aniso_float(2.5)
            EDCalculator(xmap, atoms)
        except Exception as err:
            print(err.message)

        mymap -= xmap

        mapout = clipper.CCP4MAPfile()
        mapoutPath = str(self.container.outputData.MAPOUT.fullPath)
        mapout.open_write(mapoutPath)
        mapout.export_xmap_float(mymap)
        mapout.close_write()
        return CPluginScript.SUCCEEDED
