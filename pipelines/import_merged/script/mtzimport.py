from lxml import etree
import gemmi
import numpy

from pipelines.import_merged.script.importutils import addXMLelement, ReflectionDataTypes


class ImportMTZ():
    # Import MTZ into I2 objects
    # This is an alternative to using cmtzsplit if we want to
    # select a resolution range, which cmtzsplit cannot do
    # cf ConvertCIF
    def __init__(self, inputfilename,
                 outfile, freerfile,
                 obsColLabels, contentFlag, freeRcolumnLabel,
                 resorange=None):

        print("\n>** ",inputfilename, "\n",
                 outfile, "\n", freerfile,"\n",
                 obsColLabels, contentFlag, freeRcolumnLabel,
                 resorange)

        self.status = False  # True if OK

        # Column names as a list
        self.inputcolumns = obsColLabels
        self.outputcolumns = []
        self.freercolumn = None
        if freerfile is not None:
            self.freercolumn = freeRcolumnLabel

        mtz = gemmi.read_mtz_file(inputfilename)

        # Input resolution range
        da = mtz.make_d_array()
        resorangein = (max(da), min(da))
        if resorange is not None:
            # apply resolution cutoffs first
            mtz = self.setMtzResolutionLimits(mtz, resorange)

        # Use 1st column label to get dataset
        col1 = mtz.column_with_label(obsColLabels[0])
        print ('col1', col1, type(obsColLabels[0]), obsColLabels[0])
        datasetID = col1.dataset_id

        chosenDataset = None
        for dataset in mtz.datasets:
            if dataset.id == datasetID:
                chosenDataset = dataset
                break

        contentType = ReflectionDataTypes.CONTENT_TYPES[contentFlag]
        print(contentFlag, contentType)

        # Start to construct the MTZ file for main data
        mtzout = gemmi.Mtz(with_base=True)  # with h,k,l
        mtzspecs = ReflectionDataTypes.REFLECTION_DATA[contentType]
        print(mtzspecs)
        self.start_mtzout(mtz, mtzout, mtz.title, chosenDataset, mtzspecs)

        # Get column data
        coldata = []
        ncolumns = len(obsColLabels) + 3  # including H K L
        coldata.append(mtz.column_with_label('H'))
        coldata.append(mtz.column_with_label('K'))
        coldata.append(mtz.column_with_label('L'))
        nrefs = []
        for colobs in obsColLabels:
            print(colobs)
            coldata.append(mtz.column_with_label(colobs))
            nrefs.append(len(coldata[-1]))

        print("nrefs", nrefs)
        nrefunique = nrefs[0]
        data = numpy.zeros((nrefunique, ncolumns), dtype=numpy.float)

        for i in range(ncolumns):
            data[:,i] = coldata[i]

        mtzout.set_data(data)
        mtzout.write_to_file(outfile)

        print("File ", outfile, " written, number of reflections",
              len(data))

        #  FreeR data
        if freerfile is not None and freeRcolumnLabel is not None:
            coldata = []
            ncolumns = 4  # including H K L
            coldata.append(mtz.column_with_label('H'))
            coldata.append(mtz.column_with_label('K'))
            coldata.append(mtz.column_with_label('L'))
            coldata.append(mtz.column_with_label(freeRcolumnLabel))

            # Start to construct the MTZ file for FreeR
            mtzout = gemmi.Mtz(with_base=True)  # with h,k,l
            mtzspecs = ReflectionDataTypes.REFLECTION_DATA['FreeR flag']
            print(mtzspecs)
            self.start_mtzout(mtz, mtzout, mtz.title, None, mtzspecs)

            data = numpy.zeros((nrefunique, ncolumns), dtype=numpy.float)

            for i in range(ncolumns):
                data[:,i] = coldata[i]

            mtzout.set_data(data)
            mtzout.write_to_file(freerfile)

            print("File ", freerfile, " written, number of reflections",
              len(data))

        self.makeXMLreport(nrefunique, nrefunique, contentType,
                           resorangein, resorange)
                           
        self.status = True  # True if OK

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def start_mtzout(self, mtz, mtzout, title, chosenDataset, specs):
        mtzout.title = title
        mtzout.spacegroup = mtz.spacegroup
        mtzout.set_cell_for_all(mtz.cell)
        if chosenDataset is not None:  #  None for FreeR
            mtzout.add_dataset(chosenDataset.dataset_name)
            dts = mtzout.datasets[1]
            dts.project_name = chosenDataset.project_name
            dts.crystal_name = chosenDataset.crystal_name
            dts.dataset_name = chosenDataset.dataset_name
            dts.wavelength =  chosenDataset.wavelength
            mtzout.datasets[1] = dts

        mtzout.history = ['From mtzimport']
        
        # Add columns
        if specs is not None:
            for spec in specs:
                #  MTZ column label, column type
                self.addMTZcolumnLabelType(mtzout, spec)

    def setMtzResolutionLimits(self, mtz, resorange):
        dmin = resorange[1]
        if resorange[0] > 0.0:
            dmax = resorange[0]
            mtz.set_data(mtz.array[mtz.make_d_array() >= dmin])
            mtz.set_data(mtz.array[mtz.make_d_array() <= dmax])
        else:
            mtz.set_data(mtz.array[mtz.make_d_array() >= dmin])
        return mtz

    def addMTZcolumnLabelType(self, mtz, spec):
        #  MTZ column label, column type, fix Freer type
        clabel = spec.split()[1]
        ctype = spec.split()[2]
        if ctype == 's':
            ctype = 'I'
        print("add column", clabel, ctype)
        mtz.add_column(clabel, ctype)
        self.outputcolumns.append(clabel)

    def makeXMLreport(self, nrefunique, nreffreer, contentType,
                      resorangein, resorangeout):
        
        self.importXML = etree.Element('X2MTZ')
        addXMLelement(self.importXML, 'inputcolumnames',
                      ','.join(self.inputcolumns))
        addXMLelement(self.importXML, 'outputcolumnames',
                      ','.join(self.outputcolumns))
        
        addXMLelement(self.importXML, 'nrefoutput', str(nrefunique))
        if self.freercolumn is not None:
            addXMLelement(self.importXML, 'freercolumnname',
                          self.freercolumn)
            addXMLelement(self.importXML, 'freeRwritten', str(nreffreer))

        addXMLelement(self.importXML, 'datatype', contentType)
        addXMLelement(self.importXML, 'message', 'mtz import')
        # input resolution range
        rrngxml = self.resorangeXML(resorangein)
        self.importXML.append(rrngxml)
        # Output resolution range
        if resorangeout is not None:
            rrngxml = self.resorangeXML(resorangeout, 'cutresolution')
            self.importXML.append(rrngxml)

    def getstatus(self):
        return self.status

    def getXML(self):
        return self.importXML

    def resorangeXML(self, resorange, id=None):
        # XML version of resolution range, a tuple of (dmax, dmin)
        resoxml = etree.Element('ResolutionRange')
        if id is not None:
            resoxml.set('id', id)
        if resorange[0] > 0.0:
            addXMLelement(resoxml, 'min', f"{resorange[0]:.3f}")
        if resorange[1] > 0.0:
            addXMLelement(resoxml, 'max', f"{resorange[1]:.3f}")
        return resoxml
