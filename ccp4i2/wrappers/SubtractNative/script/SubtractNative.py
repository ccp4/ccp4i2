import clipper

from ....core.CCP4MgImports import mmdb2 as mmdb, mmut
from ....core.CCP4PluginScript import CPluginScript


class SubtractNative(CPluginScript):
    TASKNAME = 'SubtractNative'   # Task name - should be same as class name and match pluginTitle in the .def.xml file
    TASKVERSION= 0.1               # Version of this plugin
    MAINTAINER = 'Martin Noble'
    ERROR_CODES = { 201 : {'description' : 'Failed to analyse output files' },
                    202 : {'description' : 'Failed applying selection ot PDB file' }
                    }
    PURGESEARCHLIST = [ [ 'hklin.mtz' , 0 ],
                       ['log_mtzjoin.txt', 0]
                       ]
    RUNEXTERNALPROCESS = False
    TASKCOMMAND="command"
    
    def __init__(self, *args, **kws):
        super(SubtractNative, self).__init__(*args, **kws)

    def processInputFiles(self):
        return CPluginScript.SUCCEEDED

    def makeCommandAndScript(self,**kw):
        return CPluginScript.SUCCEEDED
    
    def startProcess(self, command, **kw):
        print(1)
        print(2)
        mtz_file = clipper.CCP4MTZfile()
        print(3)
        hkl_info = clipper.HKL_info()
        coefficientsPath = str(self.container.inputData.MAPIN.fullPath)
        mtz_file.open_read (coefficientsPath)
        mtz_file.import_hkl_info ( hkl_info )
        sg, cell = hkl_info.spacegroup(), hkl_info.cell()
        fphidata = clipper.HKL_data_F_phi_float(hkl_info)
        mtz_file.import_hkl_data( fphidata, str("/*/*/[F,PHI]") );
        mtz_file.close_read()
        print(4)
        #Clipper will sample the output map according to Fourier theory and hte nominal resolution
        #for visualisation, it is generally nicer to make things a bit more finely sampled
        fudgedResolution = hkl_info.resolution()
        fudgedResolution.init((2./3.)*hkl_info.resolution().limit())
        mygrid=clipper.Grid_sampling ( hkl_info.spacegroup(), hkl_info.cell(), fudgedResolution )
        mymap = clipper.Xmap_float(hkl_info.spacegroup(), hkl_info.cell(), mygrid )
        print(dir(mymap))
        mymap.fft_from(fphidata)
        print(5)
        f=clipper.MMDBfile()
        f.read_file (str(self.container.inputData.XYZIN.fullPath))
        mmol=clipper.MiniMol()
        f.import_minimol (mmol)
        print(mmol)
        UsingMMDB=False
        if UsingMMDB:
            manager = mmdb.Manager()
            manager.ReadCoorFile(str(self.container.inputData.XYZIN.fullPath))
            hndl = manager.NewSelection()
            manager.SelectAtoms( hndl, 0, 0, mmdb.SKEY_NEW );
            print(8)
            selAtoms = mmdb.newPPCAtom()
            nSelAtoms = mmut.intp()
            print(nSelAtoms.value())
            print(9)
            selAtoms = mmut.GetAtomSelIndex(manager,hndl,nSelAtoms)
            print(typselAtoms)
            atoms = clipper.MMDBAtom_list( selAtoms, nSelAtoms.value() );
        atoms = mmol.atom_list()
        print(len(atoms))
        for iAtom in range(len(atoms)):
            atom = atoms[iAtom]
            modifiedOccupancy = atom.occupancy()*self.container.controlParameters.FRACTION
            print(atom.occupancy(),modifiedOccupancy)
            try:
                atom.set_occupancy(float(modifiedOccupancy))
            except Exception as err:
                print(err)
        print(atoms)
        print(hkl_info.resolution())
        print(hkl_info.cell())
        print(hkl_info.spacegroup())
        
        xmap = clipper.Xmap_float(hkl_info.spacegroup(), hkl_info.cell(), mygrid)
        print(xmap)
        print([a for a in dir(clipper) if 'ED' in a])
        try:
            EDCalculator = clipper.EDcalc_aniso_float(2.5)
            a = EDCalculator(xmap, atoms)
            print(a)
        except Exception as err:
            print(err.message)
        print(9)

        mymap -= xmap
        print(10)
        print(5)

        mapout = clipper.CCP4MAPfile()
        mapoutPath = str(self.container.outputData.MAPOUT.fullPath)
        print(6)
        mapout.open_write( mapoutPath )
        mapout.export_xmap_float( mymap )
        mapout.close_write()
        print(7)

        return CPluginScript.SUCCEEDED
    
    def processOutputFiles(self):
        with open(str(self.makeFileName('PROGRAMXML')),"w") as xmlFile:
            xmlFile.write('<SubtractNative></SubtractNative>')
        return CPluginScript.SUCCEEDED
