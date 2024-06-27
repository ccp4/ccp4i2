#  aimless_cifstats

# Extract mmcif statistics from <AIMLESS> block in xml
#  Usage:
#    class CifStatsExtractFromXML:
#      def __init__(self, aimlessxml, outputfile, blkid=None):
#        aimlessxml   XML block <AIMLESS_PIPE>
#                     Most statistics come from <AIMLESS> block,
#                     except Wilson B estimate from <CTRUNCATE>
#        outputfile   output filename for mmcif statistics
#        blkid        entry.id, default None

# Typically Aimless is run with just one dataset, but if there are multiple
# datasets, each set of statistics is written to a separate data block,
# with the blockid blkid constructed as blkid-datasetname, taken from
# the last part of the project/crystal/dataset name

import os,sys
import xml.etree.ElementTree as ET
import math

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
def getItem(xblock, item, gettext=True):
    # Extract item from xblock
    #  if gettext True, return the 1st value, else return the block
    found = xblock.findall(item)
    if len(found) == 0: return None
    if gettext:
        return found[0].text

    return found  # a list

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
class CifStatsExtractFromXML:

    def __init__(self, aimlessxml, outputfile, blkid=None):
        # aimlessxml   XML block <AIMLESS_PIPE>
        # outputfile   output filename for mmcif statistics
        # blkid        entry.id, default None

        print("CifStatsExtractFromXML: ", aimlessxml, outputfile, blkid)
        
        axml = getItem(aimlessxml, 'AIMLESS', False)
        if axml is None: return
        self.aimlessxml = axml[0]
        ctruncates = getItem(aimlessxml, 'CTRUNCATES', False)
        self.wilsonB = []
        if ctruncates is not None:
            # One for each dataset: just for Wilson B
            truncatexmlList = ctruncates[0].findall('CTRUNCATE')
            for trnc in truncatexmlList:
                self.wilsonB.append(getItem(trnc, 'DataStatistics/WilsonB'))

        if blkid is None:
            blkid = 'xxxx'
        self.blkid = blkid
        outputfile = outputfile
        self.output = CifStatsOutput(outputfile)

        resultname = 'Result'
        resultblock = self.aimlessxml.findall(resultname)
        if len(resultblock) == 0: self.fail(resultname)
        resultblock = resultblock[0]
        # for each dataset, usually one
        self.datasetResults = resultblock.findall('Dataset')

        self.ndatasets = int(getItem(self.aimlessxml,
                                 'ReflectionData/NumberDatasets'))
        if self.ndatasets > 1:
            # Multiple datasets give multiple data blocks
            self.multipleDatasets()

        else:
            #  Normal case, one dataset
            self.output.addLine('data_'+self.blkid)
            self.output.set_pair('_entry.id', self.blkid, True)
            
            self.addSpaceGroup()
            self.addCell()
            self.addWavelength()
            self.addText()

            # Summary table (Table 1)
            self.parseResultTable()
            # Statistics by shell
            self.addShellStatistics()

        # All datasets
        self.output.print()

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def parseResultTable(self, dtsidx=0):
        # Overall summary from Result block
        
        # Translation table from Aimless summary results to mmCIF
        xtocTable = {
            'ResolutionLow' : 'd_resolution_low',
            'ResolutionHigh' : 'd_resolution_high',
            'RmergeOverall' : 'pdbx_Rmerge_I_all',
            'RmeasOverall' : 'pdbx_Rrim_I_all',
            'RpimOverall' : 'pdbx_Rpim_I_all',
            'NumberObservations' : 'pdbx_number_measured_all',
            'NumberReflections' : 'number_obs',
            'MeanIoverSD' : 'pdbx_netI_over_sigmaI',
            'Completeness' : 'percent_possible_obs',
            'Multiplicity' : 'pdbx_redundancy',
            'CChalf' : 'pdbx_CC_half',
            'AnomalousCompleteness' : 'pdbx_percent_possible_anomalous',
            'AnomalousMultiplicity' : 'pdbx_redundancy_anomalous',
            'AnomalousCChalf' : 'pdbx_CC_half_anomalous'
            # 'xxx' : ' pdbx_absDiff_over_sigma_anomalous'  # not in Aimless
            }
    
        datasetblock = self.datasetResults[dtsidx]
        self.lowres = datasetblock.findall('ResolutionLow/Overall')[0].text
        self.highres = datasetblock.findall('ResolutionHigh/Overall')[0].text

        self.output.addLine('\n')
        for key  in xtocTable:
            self.addResultItem(datasetblock, key, xtocTable[key])

        # Wilson B
        if len(self.wilsonB) > 0:
            self.output.set_pair('_reflns.B_iso_Wilson_estimate',
                                 self.checkValue(self.wilsonB[dtsidx]))

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def addResultItem(self, datasetblock, xmltag, mmciftag):
        value = getItem(datasetblock, xmltag+'/Overall')
        self.output.set_pair("_reflns."+mmciftag, self.checkValue(value))

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def addSpaceGroup(self, dtsidx=0):
        # Space group name
        self.output.set_pair('_symmetry.entry_id', self.blkid, True)
        sgname = getItem(self.datasetResults[dtsidx], 'SpacegroupName')
        self.output.set_pair('_symmetry.space_group_name_H-M',
                             "'"+sgname+"'")

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def addCell(self, dtsidx=0):
        # unit cell
        cellblock = self.datasetResults[dtsidx].findall('cell')[0]
        a = getItem(cellblock, 'a')
        b = getItem(cellblock, 'b')
        c = getItem(cellblock, 'c')
        alpha = getItem(cellblock, 'alpha')
        beta  = getItem(cellblock, 'beta')
        gamma = getItem(cellblock, 'gamma')

        self.output.set_pair('_cell.entry_id', self.blkid, True)

        self.output.set_pair('_cell.length_a', a)
        self.output.set_pair('_cell.length_b', b)
        self.output.set_pair('_cell.length_c', c)
        self.output.set_pair('_cell.angle_alpha', alpha)
        self.output.set_pair('_cell.angle_beta ', beta)
        self.output.set_pair('_cell.angle_gamma', gamma)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def addWavelength(self, dtsidx=0):
        dtsblock = getItem(self.aimlessxml,
                      'ReflectionData/Dataset', False)
        wvl = getItem(dtsblock[dtsidx], 'Wavelength')
        if wvl is not None:
            self.output.set_pair('_diffrn_radiation_wavelength.id',
                                 self.blkid, True)
            self.output.set_pair('_diffrn_radiation_wavelength.wavelength',
                                 wvl)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def addShellStatistics(self, dtsidx=0):
        graphTables = {}
        # CChalf
        CChalf = self.aimlessxml.findall("CCP4Table[@id='Graph-CChalf']")[dtsidx]
        graphTables['CChalf'] = CifStatsGraphTable(CChalf, 'CChalf')

        resstats = self.aimlessxml.findall("CCP4Table[@id='Graph-StatsVsResolution']")[dtsidx]
        graphTables['ResStats'] = CifStatsGraphTable(resstats, 'ResStats')

        resstatsall = self.aimlessxml.findall("CCP4Table[@id='Graph-StatsAllVsResolution']")[dtsidx]
        graphTables['ResStatsAll'] = CifStatsGraphTable(resstatsall, 'ResStatsAll')

        completeness = self.aimlessxml.findall("CCP4Table[@id='Graph-CompletenessVsResolution']")[dtsidx]
        graphTables['Completeness'] = CifStatsGraphTable(completeness, 'Completeness')

        # Extract and compare resolution values, and make list of ranges
        # Returns list of shell boundaries
        lowreslimits, highreslimits = \
                      self.makeResolutionShells(graphTables)

#        for gtable in graphTables.values():
#         gtable.dump()

        # Translation table for converion to mmcif
        # Items:
        #  mmcif label; graphTable containing the data; column label
        # Special values:
        #   'ordinal'   serial number for row, from 1
        #   'lowres'    lowreslimits
        #   'highres'   highreslimits
       
        shellLookupTable = [\
            ['pdbx_ordinal', 'ordinal', ''],
            ['d_res_low', 'lowres', ''],
            ['d_res_high', 'highres', ''],
            ['number_measured_all', 'Completeness', 'Nmeas'],
            ['number_unique_all', 'Completeness', 'Nref'],
            ['percent_possible_all', 'Completeness', '%poss'],
            ['Rmerge_I_all', 'ResStatsAll', 'RmrgOv'], 
            ['pdbx_Rrim_I_all', 'ResStatsAll', 'RmeasOv'],
            ['pdbx_Rpim_I_all', 'ResStatsAll', 'RpimOv'],
            ['pdbx_redundancy', 'Completeness', 'Mlplct'],
            ['pdbx_netI_over_sigmaI_all', 'ResStats', 'Mn(I/sd)'],
            ['pdbx_CC_half', 'CChalf', 'CC1/2'],
            ['pdbx_percent_possible_anomalous', 'Completeness', 'AnoCmp'],
            ['pdbx_redundancy_anomalous', 'Completeness', 'AnoMlt'],
            ['pdbx_CC_half_anomalous', 'CChalf', 'CCanom']\
            ]


        # self.nvalues is length of each column
        colvalues = []   # List of wanted columns
        for shellLookup in shellLookupTable:
            cifcode = shellLookup[0]
            table = shellLookup[1]
            colname = shellLookup[2]

            # Specials
            if shellLookup[1] == 'ordinal':
                # ordinal
                nums = list(range(1,self.nvalues+1))
                nums = [str(i) for i in nums]
                colvalues.append(nums)
            elif shellLookup[1] == 'lowres':
                colvalues.append(lowreslimits)
            elif shellLookup[1] == 'highres':
                colvalues.append(highreslimits)
            else:
                gt = graphTables[shellLookup[1]]
                coldata = gt.getColumn(shellLookup[2])
                colvalues.append(coldata)

        #        for cd in colvalues:
        #            print('\n', cd)

        tags = []
        for lookup in shellLookupTable:
            tags.append('_reflns_shell.'+lookup[0])

    
        # write loop headers
        self.output.addLine('\nloop_')

        for tag in tags:
            self.output.addLine(tag)
            
        ncolumns = len(colvalues)
        nvalues = len(colvalues[0])

        for j in range(nvalues):  # rows
            s = ''
            for i in range(ncolumns):  # columns
                val = self.checkValue(colvalues[i][j])
                s += ' ' + '{: <5}'.format(val)
            self.output.addLine(s)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def makeResolutionShells(self, graphTables):
        # Each graph table has resolution shells: compare them

        ntables = len(graphTables)
        rdvalues = []
        nvalues = []
        
        for graphtable in graphTables.values():
            cols = graphtable.getColumns()
            names = cols.getNames()
            rdvalues.append(cols.getColumn('1/d^2'))
            nvalues.append(len(rdvalues[-1]))

        if len(nvalues) != len(rdvalues):
            print('Unequal Dmid lengths:',len(nvalues), len(rdvalues))
            exit(1)
        self.nvalues = nvalues[0]

        # half width of ranges in 1/d^2
        halfwidth = 0.5*(float(rdvalues[0][1])-float(rdvalues[0][0]))
        nranges = nvalues[0]

        lowreslimits = [0.0]*nranges  # lower limits of each range
        highreslimits = [0.0]*nranges
        
        for i in range(1,nranges):
            dstar = float(rdvalues[0][i]) - halfwidth
            lowreslimits[i] = str('{:.3f}'.format(1.0/math.sqrt(dstar)))
        for i in range(nranges-1):
            highreslimits[i] = lowreslimits[i+1]

        # low resolution limit
        lowreslimits[0] = '{: <5}'.format(self.lowres.strip())
        # high resolution limit
        highreslimits[-1] = '{: <5}'.format(self.highres.strip())
        
        return lowreslimits, highreslimits

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def addText(self, ordinal='1'):
        #  Write an explanatory text block
        self.output.set_pair('_reflns.pdbx_ordinal', ordinal, True)
        text = []
        text.append('_reflns.details')
        text.append(';')
        text.append('PDB mmCIF categories _reflns.*_obs and _reflns_shell.*_obs are not used, as they seem to be based on a long-discredited idea that weak reflections should be discarded: this was never a Good Idea, and should not be used now')
        text.append(';')

        for line in text:
            self.output.addLine(line)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def getDatasetNames(self):
        reflectionblock = getItem(self.aimlessxml, 'ReflectionData', False)[0]
        datasetinfolist = reflectionblock.findall('Dataset')
        dsetnames = []
        for dset in datasetinfolist:
            dname = dset.attrib['name'].split('/')[-1]
            dsetnames.append(dname)
        return dsetnames

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def multipleDatasets(self):
        # Multiple datasets, get list of dataset names
        dsetnames = self.getDatasetNames()
        blkid = self.blkid # save base name

        for idts in range(self.ndatasets):   # Loop datasets
            self.blkid = blkid+'-'+dsetnames[idts]

            self.output.addLine('\ndata_'+self.blkid)
            
            self.output.set_pair('_entry.id', self.blkid, True)
            self.addSpaceGroup(idts)
            self.addCell(idts)
            self.addWavelength(idts)
            self.addText(str(idts+1))

            # Summary table (Table 1)
            self.parseResultTable(idts)
            # Statistics by shell
            self.addShellStatistics(idts)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def checkValue(self, svalue):
        #  Replace any "-" fields with "?"
        if svalue == '-':
            return '?'
        return svalue

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def fail(self, elementname):
        print('FATAL: Element', elementname, ' not found')
        exit(1)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
class CifStatsGraphTable:
    def __init__(self, graphblock, label):
        #  construct from CCP4Table
        self.label = label
        self.title = graphblock.attrib['title']
        self.headers = getItem(graphblock, 'headers')
        data = getItem(graphblock, 'data')
        self.columns = CifStatsColumns(self.headers, data)

        #print('\n'+self.label)
        names = self.columns.getNames()
        for name in names:
            col = self.columns.getColumn(name)
            #print('Name:', name)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def getTitle(self):
        return self.title

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def getLabel(self):
        return self.label

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def getColumns(self):
        return self.columns
        
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def getColumn(self, name):
        return self.columns.getColumn(name)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def dump(self):
        print('\n',self.title)
        s = ' '.join(self.headers.split())
        print(s)
        names = self.columns.getNames()
        cols = []
        for name in names:
            cols.append(self.columns.getColumn(name))

        for i in range(len(cols[0])):
            s = ''
            for j in range(len(cols)):
                s += ' '+cols[j][i]
            print(s)
        

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
class CifStatsColumns:
    def __init__(self, headers, table):
        self.colnames = headers.split()
        ncols = len(self.colnames)
        self.coldata = {}
        for name in self.colnames:
            self.coldata[name] = []

        lines = table.splitlines()
        for line in lines:
            fields = line.split()
            if len(fields) > 0:
                if len(fields) != ncols:
                    print("Wrong length", len(fields), ncols)
                for i, d in enumerate(fields):
                    self.coldata[self.colnames[i]].append(fields[i])

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def getNames(self):
        # Return list of column names
        return self.colnames
                    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def getColumn(self, name=None):
        if name in self.coldata:
            return  self.coldata[name]
        return self.coldata[self.getNames()[0]]

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
class CifStatsOutput:
    # Accumulate all mmcif output as list of lines (strings)
    def __init__(self, outputfile=None):
        self.outputfile = outputfile
        self.output = []

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def addLine(self, line):
        self.output.append(line)
        
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def set_pair(self, tag, value, newline=False):
        s = tag+'  '+value
        if newline: s = '\n'+s
        self.addLine(s)
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def print(self):
        if self.outputfile is not None:
            ofl = open(self.outputfile, 'w')
            for line in self.output:
                print(line, file=ofl)
        else:
            for line in self.output:
                print(line)


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#  For stand-alone running and testing

if __name__ == "__main__":
    args = sys.argv
    # Parse arguments from string list
    import argparse
    parser = argparse.ArgumentParser(args)
    if len(args) <= 1:
        print("Usage: python aimlessstatscif.py xmlfilename [options]")
        print("Options: -o output filename; -id  block name for output; -h help")
        exit(1)
            
    parser.add_argument('filename')
    parser.add_argument('-o', '--outputfile')
    parser.add_argument('-id', '--id')

    values = parser.parse_args()
    
    filename = values.filename
    blkid = values.id
    outputfile = values.outputfile


    # XML block is probably always AIMLESS
    aimlessname = 'AIMLESS_PIPE'

    tf = ET.parse(filename)
    #print(filename, tf.getroot().tag)

    if  tf.getroot().tag == aimlessname:
        aimlessxml = tf.getroot()
    else:
        aimlessxml = tf.findall(aimlessname)
    if type(aimlessxml) == list:
        aimlessxml = aimlessxml[0]

    CifStatsExtractFromXML = CifStatsExtractFromXML(aimlessxml, outputfile, blkid)
    
