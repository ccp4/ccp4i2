from ccp4i2.core.CCP4PluginScript import CPluginScript
from ccp4i2.core.mgimports import mmdb2 as mmdb
from ccp4i2.core.mgimports import mmut


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
        import clipper
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
        '''
        print([a for a in dir(clipper) if 'FFT' in a.upper()])
        try:
            mymap.fft_to_float(fphidata)
        except Exception as err:
            print(dir(mymap))
            print("Error FFTing to {}".format(str(err)))
                                              
        print(5.5)


        mtzout = clipper.CCP4MTZfile()
        mtzoutPath = str(self.container.outputData.MAPCOEFFSOUT.fullPath)
        mtzout.open_write( mtzoutPath )
        mtzout.export_xmap( fphidata )
        mtzout.close_write()

              /* mmdb part of the calculation */

  CMMDBManager mmdb;
  mmdb.ReadPDBASCII( "input.pdb" );           // read pdb file
  int hndl = mmdb.NewSelection();             // make selection handle
  mmdb.SelectAtoms( hndl, 0, 0, SKEY_NEW );   // select all atoms
  PPCAtom psel;
  int nsel;
  mmdb.GetSelIndex( hndl, psel, nsel );       // get the selection

  /* Clipper part of the calculation */

  clipper::HKL_info hkls;       // make reflection lists for result
  /* *********************************************************** */
  /* NOTE: we need to initialise the reflection list 'hkls' here */
  /* *********************************************************** */
  clipper::HKL_data<clipper::data32::F_phi> fphi(hkls);  // and data list
  clipper::MMDBAtom_list atoms( psel, nsel );            // make atom list
  clipper::SFcalc_aniso_fft<float>( fphi, atoms );       // and do SF calc'''
        
        '''
      // prepare target map
  const HKL_info& hkls = fphidata.base_hkl_info();
  const Grid_sampling grid( hkls.spacegroup(), hkls.cell(), hkls.resolution() );
  Xmap<float> xmap( hkls.spacegroup(), hkls.cell(), grid );

  // work out how big a box we need to calc density over for each atom
  Grid_range gd( hkls.cell(), grid, 3.0 );
  Xmap<float>::Map_reference_coord i0, iu, iv, iw;
  // loop over atoms
  for ( int i = 0; i < atoms.size(); i++ ) if ( !atoms[i].is_null() ) {
    AtomShapeFn sf( atoms[i] );  // get atom shape fn
    Coord_frac uvw = atoms[i].coord_orth().coord_frac( hkls.cell() );
    Coord_grid g0 = uvw.coord_grid( grid ) + gd.min();
    Coord_grid g1 = uvw.coord_grid( grid ) + gd.max();
    i0 = Xmap<float>::Map_reference_coord( xmap, g0 );
    // sum all map contributions from this atoms
    for ( iu = i0; iu.coord().u() <= g1.u(); iu.next_u() )
      for ( iv = iu; iv.coord().v() <= g1.v(); iv.next_v() )
        for ( iw = iv; iw.coord().w() <= g1.w(); iw.next_w() )
          xmap[iw] += sf.rho( iw.coord_orth() );
  }

'''
        return CPluginScript.SUCCEEDED
    
    def processOutputFiles(self):
        with open(str(self.makeFileName('PROGRAMXML')),"w") as xmlFile:
            xmlFile.write('<SubtractNative></SubtractNative>')
        return CPluginScript.SUCCEEDED

