# mtzimport.py

import sys
import math
import numpy
import gemmi

from lxml import etree

try:
    from importutils import *
except:
    from pipelines.import_merged.script.importutils import *

class ImportMTZ():
    # Import MTZ into I2 objects
    # This is an alternative to using cmtzsplit if we want to
    # select a resolution range, which cmtzsplit cannot do
    # cf ConvertCIF
    #
    # Optionally exclude reflections for which all selected columns
    #   are missing == NaN
    # That may involve changing the resolution limits

    # inputfilename    input MTZ file name
    # outfile          output filename for selected data
    # freerfile        output filename for FreeR data if present
    # obsColLabels     column labels selected from input file
    # contentFlag      used to select I+-, F+-, Imean, Fmean = 1,2,3,4
    # freeRcolumnLabel column label for FreeR, if present
    # resorange=None   optional resolution cutoff
    # excludemissing=False
    #           if True, omit reflections with all selected columns missing

    def __init__(self, inputfilename,
                 outfile, freerfile,
                 obsColLabels, contentFlag, freeRcolumnLabel,
                 resorange=None, excludemissing=False):

        print("\n>** Input file: ",inputfilename, "\nOutput files: ",
              outfile, "\n", freerfile,"\n",
              "Column labels: ",obsColLabels,
              " Content flag: ", contentFlag,
              "\nResolution range: ",resorange)
        if excludemissing:
            print("Exclude reflections with all selected columns missing")
        
        self.status = False  # True if OK

        # Column names as a list
        self.inputcolumns = obsColLabels
        self.outputcolumns = []
        mtz = gemmi.read_mtz_file(inputfilename)

        self.freercolumn = None
        if freerfile is not None:
            if mtz.rfree_column() is not None:
                self.freercolumn = freeRcolumnLabel
                
        if self.freercolumn is None:
            print("No FreeR in file")
        else:
            print("\nFreerLabel: ", self.freercolumn)

        # Input resolution range
        self.d_array = mtz.make_d_array()
        resorangein = (max(self.d_array), min(self.d_array))
        if resorange is not None:
            # apply resolution cutoffs first
            mtz = self.setMtzResolutionLimits(mtz, resorange)
        self.maxResIncluded = resorangein[1]
            
        # Use 1st column label to get dataset
        col1 = mtz.column_with_label(obsColLabels[0])
        print ('col1', col1, type(obsColLabels[0]), obsColLabels[0])
        datasetID = col1.dataset_id

        chosenDatasetIndex = -1
        chosenDataset = None
        for i in range(len(mtz.datasets)):
            if mtz.datasets[i].id == datasetID:
                chosenDatasetIndex = i
                chosenDataset = mtz.datasets[i]
                break
            
        # from importutils.py
        contentType = ReflectionDataTypes.CONTENT_TYPES[contentFlag]
        print("ContentFlag", contentFlag, ", contentType:", contentType)

        # Start to construct the MTZ file for main data
        mtzout = gemmi.Mtz(with_base=True)  # with h,k,l
        mtzspecs = ReflectionDataTypes.REFLECTION_DATA[contentType]
        #print(mtzspecs)
        self.start_mtzout(mtz, mtzout, mtz.title, chosenDataset, mtzspecs)

        # Get column data
        coldata = []
        ncolumns = len(obsColLabels) + 3  # including H K L
        coldata.append(mtz.column_with_label('H'))
        coldata.append(mtz.column_with_label('K'))
        coldata.append(mtz.column_with_label('L'))
        nrefs = []
        for colobs in obsColLabels:
            #print(colobs)
            coldata.append(mtz.column_with_label(colobs))
            nrefs.append(len(coldata[-1]))

        nrefunique = nrefs[0]
        self.nrefinput = nrefunique
        data = numpy.zeros((nrefunique, ncolumns), dtype=float)

        for i in range(ncolumns):
            data[:,i] = coldata[i]

        mtzout.set_data(data)  # Columns selected
        
        self.wanted = None
        self.nrefremoved = 0  # set in stripMissing
        # Always call if only to count missing reflections
        mtzout = self.stripMissing(data, mtzout, excludemissing)
        if excludemissing:
            print("Number of reflections removed because all missing:",
                  self.nrefremoved)
            
        nrefout = len(mtzout.array)
        mtzout.write_to_file(outfile)

        print("\nFile ", outfile, " written, number of reflections",
              nrefout)

        da = mtzout.make_d_array()
        resorangeout = (max(da), min(da))

        #  FreeR data
        if freerfile is not None and freeRcolumnLabel is not None \
           and self.freercolumn is not None:
            coldata = []
            ncolumns = 4  # including H K L
            coldata.append(mtz.column_with_label('H'))
            coldata.append(mtz.column_with_label('K'))
            coldata.append(mtz.column_with_label('L'))
            coldata.append(mtz.column_with_label(freeRcolumnLabel))

            # Start to construct the MTZ file for FreeR
            mtzout = gemmi.Mtz(with_base=True)  # with h,k,l
            mtzspecs = ReflectionDataTypes.REFLECTION_DATA['FreeR flag']
            #print(mtzspecs)
            
            self.start_mtzout(mtz, mtzout, mtz.title, None, mtzspecs)

            data = numpy.zeros((nrefunique, ncolumns), dtype=float)

            for i in range(ncolumns):
                data[:,i] = coldata[i]

            mtzout.set_data(data)

            if excludemissing:
                mtzout = self.stripMissing(data, mtzout, excludemissing)

            mtzout.write_to_file(freerfile)

            print("\nFile ", freerfile, " written, number of reflections",
                  nrefout)

        
        self.makeXMLreport(nrefout, nrefout, contentType, resorange,
                           resorangein, resorangeout, excludemissing)
                           
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

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def setMtzResolutionLimits(self, mtz, resorange):
        dmin = resorange[1]
        if resorange[0] > 0.0:
            dmax = resorange[0]
            mtz.set_data(mtz.array[mtz.make_d_array() >= dmin])
            mtz.set_data(mtz.array[mtz.make_d_array() <= dmax])
        else:
            mtz.set_data(mtz.array[mtz.make_d_array() >= dmin])
        return mtz

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def addMTZcolumnLabelType(self, mtz, spec):
        #  MTZ column label, column type, fix Freer type
        clabel = spec.split()[1]
        ctype = spec.split()[2]
        if ctype == 's':
            ctype = 'I'
        print("add column", clabel, ctype)
        mtz.add_column(clabel, ctype)
        self.outputcolumns.append(clabel)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def stripMissing(self, data, mtzout, excludemissing):
        # Count reflections where all selected columns are NaN
        # If excludemissing True, then exclude them
        # Sets self.nrefremoved  = number removed

        nrefremoved = 0
        if self.wanted is None:
            self.maxResIncluded = 10000.0
            # for each reflection, are the elements NaN (missing)?
            # Compile self.wanted list of included reflections
            testnan = numpy.isnan(data[:,3:])
            # compile flags for each reflection, = False if all selected
            #  columns are NaN
            nref = len(data)
            self.wanted = numpy.zeros(nref, dtype=bool)

            for j in range(nref):
                if testnan[j,:].all():
                    # all columns = NaN, reject
                    self.wanted[j] = False
                    nrefremoved += 1
                else:
                    self.wanted[j] = True
                    self.maxResIncluded = min(self.maxResIncluded, self.d_array[j])
                    
            self.nrefremoved = nrefremoved

        if excludemissing:
            mtzout = mtzout.filtered(self.wanted)  #  remove null reflections

        return mtzout
        
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def makeXMLreport(self, nrefout, nreffreer, contentType, resorange,
                      resorangein, resorangeout, excludemissing):
        
        self.importXML = etree.Element('X2MTZ')
        addXMLelement(self.importXML, 'inputcolumnames',
                      ','.join(self.inputcolumns))
        addXMLelement(self.importXML, 'outputcolumnames',
                      ','.join(self.outputcolumns))
        
        addXMLelement(self.importXML, 'nrefinput', str(self.nrefinput))
        addXMLelement(self.importXML, 'nrefoutput', str(nrefout))

        if self.nrefremoved > 0:   # flagged or removed
            if excludemissing:
                addXMLelement(self.importXML, 'nrefexcluded',
                              str(self.nrefremoved))
            else:
                addXMLelement(self.importXML, 'nrefflagged',
                              str(self.nrefremoved))                

        # Maximum resolution for reflections not flagged or excluded
        addXMLelement(self.importXML, 'maxResolutionAccepted',
                      "{:.3f}".format(self.maxResIncluded))
                
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
            if resorangeout != resorangein:
                rrngxml = self.resorangeXML(resorangeout, 'cutresolution')
                self.importXML.append(rrngxml)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def getstatus(self):
        return self.status

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def getXML(self):
        return self.importXML

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def resorangeXML(self, resorange, id=None):
        # XML version of resolution range, a tuple of (dmax, dmin)
        resoxml = etree.Element('ResolutionRange')
        if id is not None:
            resoxml.set('id', id)
        if resorange[0] > 0.0:
            addXMLelement(resoxml, 'min', "{:.3f}".format(resorange[0]))
        if resorange[1] > 0.0:
            addXMLelement(resoxml, 'max', "{:.3f}".format(resorange[1]))
        return resoxml

