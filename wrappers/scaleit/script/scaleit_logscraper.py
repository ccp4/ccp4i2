import os, sys
import lxml
#from lxml import etree
from xml.etree import ElementTree as ET

from pathlib import Path

CCP4I2_ROOT = os.environ.get("CCP4I2_ROOT",
                             str(Path(lxml.__file__).parents[1] / "ccp4i2"))
sys.path.append(str(CCP4I2_ROOT))

from core.CCP4PluginScript import CPluginScript
from core import CCP4ErrorHandling
from core import CCP4Utils

class makeGraphs:
    def __init__(self, logfilename):
        self.logfilename = logfilename
        self.scaleitgraphs = ET.Element('SCALEITGRAPHS')
        # Extract graphs from log file
        self.scrapeSmartieGraphs(self.scaleitgraphs)

    def getXML(self):
        return self.scaleitgraphs

    # - - - - - - - - -  - - - - - - - - -  - - - - - - - - - 
    def scrapeSmartieGraphs(self, smartieNode):
        from core import CCP4Utils
        smartiePath = os.path.join(CCP4Utils.getCCP4I2Dir(),'smartie')
        sys.path.append(smartiePath)
        import smartie
        
        logfile = smartie.parselog(self.logfilename)
        for smartieTable in logfile.tables():
            if smartieTable.ngraphs() > 0:
                tableelement = \
                          self.xmlForSmartieTable(smartieTable, smartieNode)
        
        return
    
    # - - - - - - - - -  - - - - - - - - -  - - - - - - - - - 
    def xmlForSmartieTable(self, table, parent):
        from pimple import MGQTmatplotlib
        tableetree = MGQTmatplotlib.CCP4LogToEtree(table.rawtable())
        parent.append(tableetree)
        return tableetree


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
class scaleitLogScraper:
    # cf refmacLogscraper, but this is simpler
    def __init__(self):
        self.datasetsinfile = None
        self.nderivatives = 0
        self.datasetsused = {}
        self.nativedataset = None
        self.anomused = []
        self.scales = None
        self.differenceLimits = []  # for each derivative
        self.resolutionTotals = []
        self.resolutionHeaders = []
        self.resolutionTotalLines = []
        self.normalProbabilityHeaders = None
        self.resolutionMax = None
        
        # 1st column label for each used derivative
        # an index into self.datasetsused
        self.derivativeIndex = []
        self.derivativeNames = []  # dnames
        
        # blockHandlers recognise triggerText and then call parseFunction
        self.blockHandlers = [{'triggerText':'Data line---',
                               'parseFunction':self.parseControlData},
                              {'triggerText':'OPENED INPUT MTZ FILE',
                               'parseFunction':self.parseMTZfile},
                              {'triggerText':'Derivative:   ',
                               'parseFunction':self.parseDerivativeNames},
                              {'triggerText':'APPLICATION OF SCALES ',
                               'parseFunction':self.parseScales},
                              {'triggerText':'Differences greater ',
                               'parseFunction':self.parseDifferences},
                              {'triggerText':'Normal probability ',
                               'parseFunction':self.parseNormalProbability},
                              {'triggerText':'Analysis v resolution ',
                               'parseFunction':self.parseResolutionAnalysis}
                              ]


    # - - - - - - - - -  - - - - - - - - -  - - - - - - - - - 
    def process(self, logfilename):
        with open(logfilename,'r') as logFile:
            filetext = logFile.read()

        self.lines = filetext.splitlines()
        self.nlines = len(self.lines)
        self.iline = 0
        while self.iline < self.nlines:
            line = self.lines[self.iline]
            for blockHandler in self.blockHandlers:
                if blockHandler['triggerText'] in line:
                    blockHandler['parseFunction'](line)
            self.iline += 1

        # Add the graphs
        self.graphs = makeGraphs(logfilename)

    # - - - - - - - - -  - - - - - - - - -  - - - - - - - - - 
    def fileInfo(self):
        # Return information about input files to Scaleit
        nativeDname = self.datasetsused[self.nativedataset][1]
        derivativeDnames = []
        for i in range(self.nderivatives):
            lab = self.derivativeIndex[i]
            derivativeDnames.append(self.datasetsused[lab][1])
        return nativeDname, derivativeDnames

        # - - - - - - - - -  - - - - - - - - -  - - - - - - - - - 
    def makeXML(self):
        # Returns an XML block SCALEITLOG
        #print(">>>logscraper makeXML")

        xmlroot = ET.Element('SCALEITLOG')

        if self.resolutionMax is not None:
            addElement(xmlroot, 'ResolutionMax', self.resolutionMax)

        # "native" is the first dataset to which others are compared
        native = self.datasetsused[self.nativedataset]  # [ID, dname]
        addElement(xmlroot, 'NativeDname', native[1])
        addElement(xmlroot, 'NativeColName', self.nativedataset)

        addElement(xmlroot, 'Nderivatives', str(self.nderivatives))
        for i in range(self.nderivatives):
            lab = self.derivativeIndex[i]
            e = ET.Element('Derivative')
            addElement(e, 'Name', self.datasetsused[lab][1])
            addElement(e, 'ColName', lab)
            
            m = 'false'
            if self.anomused[i]:
                # Anomalous
                m = "true"
            addElement(e, 'Anomalous', m)
            # Scales, B
            scales = self.scales[i]
            addElement(e,'ScaleFactor', scales[2])
            if len(scales) > 4:
                # Anisotropic B
                m = 'AnisoB'
            else:
                m = 'IsoB'
            addElement(e, m, ' '.join(scales[3:]))

            # Analysis v. resolution
            er = ET.Element('ResolutionAnalysis')
            addElement(er, 'Headers', ' '.join(self.resolutionHeaders[i]))
            addElement(er, 'Totals', ' '.join(self.resolutionTotalLines[i]))
            for key, value in self.resolutionTotals[i].items():
                addElement(er, self.validTag(key), value)
            e.append(er)
            xmlroot.append(e)
        # End loop derivatives

        # Large differences
        xmlroot.append(self.largeDifferenceList.makeXML())

        # Normal probability analysis if present
        if self.normalProbabilityHeaders is not None:
            e = ET.Element('NormalProbability')
            addElement(e, 'Headers', ' '.join(self.normalProbabilityHeaders))
            addElement(e, 'Centric', ' '.join(self.normalProbCentric))
            addElement(e, 'Acentric', ' '.join(self.normalProbAcentric))
            xmlroot.append(e)

        # Add graphs
        graphxml = self.graphs.getXML()
        xmlroot.append(graphxml)
        
        return xmlroot
    # - - - - - - - - -  - - - - - - - - -  - - - - - - - - -
    def validTag(self, tag):
        # Strip XML invalid characters from tag
        invalid = ['<', '>', '(', ')', '|']
        t = ''
        for c in tag:
            if c not in invalid:
                t +=c
        # Check for '|', if so prepend 'Mod'
        if '|' in tag:
            t = 'Mod'+t
        if '<' in tag:
            t = 'Mean'+t
        return t

    # - - - - - - - - -  - - - - - - - - -  - - - - - - - - - 
    def parseControlData(self, line):
        print("parseControlData", line)
        if 'RESOLUTION' in line:
            self.resolutionMax = line.split()[3]
    # - - - - - - - - -  - - - - - - - - -  - - - - - - - - - 
    def parseMTZfile(self, line):
        # Extract information about the input MTZ file
        # Sets:
        #  self.datasetsinfile  dictionary indexed on 1st column label
        self.datasetsinfile = {}   # indexed by 1st column label
        line = self.findLine('Number of Datasets', 0)
        self.ndatasetsinfile = int(line.split()[5])
        line = self.findLine('Dataset ID', 2)
        datasetinfo = []
        for i in range(self.ndatasetsinfile):
            line = self.lines[self.iline]
            datasetID = int(line.split()[0])
            self.iline += 2
            # Dataset name
            dname = self.lines[self.iline].split()[0]
            datasetinfo.append([datasetID, dname])
            self.iline += 3
            
        line = self.findLine(' Column Labels', 2)
        self.labels = line.split()
        line =  self.findLine('Associated datasets', 2)
        associateddatasets = [int(x) for x in line.split()]
        line = self.findLine('Space group', 0)

        # Associate labels with datasets
        # 1) get first column label for each dataset
        firstlabels = []
        currentid = 0
        for i, dtsid in enumerate(associateddatasets):
            if dtsid != currentid:
                # New ID
                firstlabels.append([dtsid, self.labels[i]])
                currentid = dtsid
        # 2) build self.datasetsinfile dictionary, indexed on first label
        for firstlabel in firstlabels:
            # For this dataset, we want datasetID, datasetname
            for dinfo in datasetinfo:
                if dinfo[0] == firstlabel[0]:
                    # found ID
                    self.datasetsinfile[firstlabel[1]] = dinfo

        # self.datasetsinfile if dictionary FirstLabel: [ID, dname]
                    
    # - - - - - - - - -  - - - - - - - - -  - - - - - - - - - 
    def parseDerivativeNames(self, line):
        # Get column names and datasetnames for datasets used
        #  Add to self.datasetsused dictionary
        fields = line.split()
        # "Native"
        idx1 = fields.index('FP=')
        lab1 = fields[idx1+1]
        if lab1 not in self.datasetsused:
            self.datasetsused[lab1] = self.datasetsinfile[lab1]
            self.nativedataset = lab1
        # "Derivative"
        idx2 = fields.index('FPH=')
        lab2 = fields[idx2+1]
        if lab2 in self.datasetsinfile:
            self.datasetsused[lab2] = self.datasetsinfile[lab2]
        else:
            # Derivative is not in its own dataset, so fake it
            self.datasetsused[lab2] = [lab2, 'DERIVATIVE']
        # this routine entered for each derivative,
        # so  self.derivativeIndex
        self.derivativeIndex.append(lab2)
        self.derivativeNames.append(self.datasetsused[lab2][1])
        # Number of dataset pairs analysed
        self.nderivatives = len(self.datasetsused) - 1
        hasanom = False
        if len(fields) > 6:
            # This derivative has anomalous
            hasanom = True
        self.anomused.append(hasanom)
        
    # - - - - - - - - -  - - - - - - - - -  - - - - - - - - - 
    def parseScales(self, line):
        # Get scale factors and [an]isotropic B-factors
        self.scales = []
        for i in range(self.nderivatives):
            line = self.findLine('Derivative', 0)
            scales = line.split()
            if len(scales) > 4:
                # anisotropic
                self.iline += 1
                line = self.findLine('Derivative', 0)
                scales = line.split()
            self.iline += 1
            self.scales.append(scales)
        
    # - - - - - - - - -  - - - - - - - - -  - - - - - - - - - 
    def parseDifferences(self, line):
        # Difference limits and list of Large differences
        self.iline -= 3  # back off a bit

        for ideriv in range(self.nderivatives):  # Loop derivatives
            # Isomorphous difference limits
            line = self.findLine('Isomorphous Differences')
            isolimits = self.getDifferenceLimits(line)
            anomlimits = None
            if self.anomused[ideriv]:
                # Anomalous limits
                line = self.findLine('Anomalous Differences')
                anomlimits = self.getDifferenceLimits(line)
            self.differenceLimits.append([isolimits, anomlimits])

        # Get list of large differences
        self.largeDifferences()
        
    # - - - - - - - - -  - - - - - - - - -  - - - - - - - - - 
    def getDifferenceLimits(self, line):
        #  1) RMSDIF
        line = self.findLine('Differences greater', 0)
        rmslimit = line.split()[3]
        #  2) acceptable difference
        line = self.nextLine()
        acceptablediff = line.split()[6]
        #  2) maximum difference
        line = self.nextLine()
        maxdiff = line.split()[2]
        limits = [rmslimit, acceptablediff, maxdiff]
        return limits
    # - - - - - - - - -  - - - - - - - - -  - - - - - - - - - 
    def largeDifferences(self):
        # Get list of large differences
        line = self.findLine('Large differences', 2)
        # List of columns in table
        labels = line.split()
        ncols = len(labels)
        self.largeDifferenceList = \
                LargeDifferenceList(self.nderivatives,
                                    self.derivativeNames,
                                    self.differenceLimits,
                                    labels)
        
        self.iline += 3
        self.iline += \
            self.largeDifferenceList.extractList(self.lines[self.iline:])

    # - - - - - - - - -  - - - - - - - - -  - - - - - - - - - 
    def parseNormalProbability(self, line):
        # Normal probability if present (only for one derivative)
        line = self.findLine("1/resol^",1)
        line = self.findLine("1/resol^")
        # Remove resolution and $$
        self.normalProbabilityHeaders = line.split()[2:-1]
        line = self.findLine("Total (cent)", 0)  # centric
        self.normalProbCentric = line.split()[2:]
        line = self.findLine("Total (acen)", 0)  # acentric
        self.normalProbAcentric = line.split()[2:]
        
    # - - - - - - - - -  - - - - - - - - -  - - - - - - - - - 
    def parseResolutionAnalysis(self, line):
        # Resolution analysis table, the totals for each derivative

        # Which derivative is this?
        l = line.split()
        lab = l[l.index('FPH=')+1]
        # Derivative number
        ideriv = self.derivativeIndex.index(lab)

        line = self.findLine("$$", 1)
        line = self.findLine("$$", 1)
        line = self.findLine("1/resol^")
        headers = line.split()[2:-1] # Remove resolution and $$
        self.resolutionHeaders.append(headers)
        line = self.findLine("THE TOTALS")
        values = line.split()[2:]
        self.resolutionTotalLines.append(values)
        # Store as dictionary
        indices = [0,3,4,5,6,7,8]
        if self.anomused[ideriv]:
            indices.extend([9,10,11,12,13])
        totalValues = {}
        for i in indices:
            lab = headers[i]
            value = values[i]
            totalValues[lab] = value
        self.resolutionTotals.append(totalValues) # for each derivative
    # - - - - - - - - -  - - - - - - - - -  - - - - - - - - - 
    def findLine(self, tag, nskip=0):
        while self.iline < self.nlines:
            line = self.lines[self.iline]
            if tag in line:
                #print("tag found", tag, line)
                self.iline += nskip
                return self.lines[self.iline]
            self.iline += 1
        return None
    # - - - - - - - - -  - - - - - - - - -  - - - - - - - - - 
    def nextLine(self, nskip=0):
        self.iline += nskip + 1
        if self.iline < self.nlines:
            line = self.lines[self.iline]
            return line
        return None
    # - - - - - - - - -  - - - - - - - - -  - - - - - - - - - 
    def addElement(self, containerXML, elementname, elementtext):
        #print 'addElement', elementname, type(elementtext), elementtext 
        e2 = ET.Element(elementname)
        e2.text = elementtext
        containerXML.append(e2)
        
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
class LargeDifferenceList():
    # List of large differences
    def __init__(self, nderivatives, derivativeNames,
                 differenceLimits, labels):
        self.nderivatives = nderivatives
        self.derivativeNames = derivativeNames
        self.differenceLimits = differenceLimits
        self.collabels = labels

    # - - - - - - - - -  - - - - - - - - -  - - - - - - - - - 
    def extractList(self, lines):
        self.data = []
        iline = 0
        while True:
            line = lines[iline]
            if "SUMMARY_END" in line:
                break
            if "=" not in line:
                self.data.append(line)
            iline += 1
        return iline

    # - - - - - - - - -  - - - - - - - - -  - - - - - - - - - 
    def makeXML(self):
        largeDiffXML = ET.Element('LARGE_DIFFERENCES')

        for i in range(self.nderivatives):
            e = ET.Element('Derivative')
            addElement(e, 'Name', self.derivativeNames[i])
            limits = self.differenceLimits[i]
            el = ET.Element('IsomorphousDifferences')
            self.formatLimits(limits[0], el)
            e.append(el)
            if limits[1] is not None:
                # Anomalous
                el = ET.Element('AnomalousDifferences')
                self.formatLimits(limits[1], el)
                e.append(el)

            largeDiffXML.append(e)

        # List of differences
        labels = ' '.join(self.collabels)
        addElement(largeDiffXML,'TableColumnLabels', labels)

        text = '\n'+'\n'.join(self.data)+'\n'
        addElement(largeDiffXML,'LargeDiffList', text)

        return largeDiffXML
    # - - - - - - - - -  - - - - - - - - -  - - - - - - - - - 
    def formatLimits(self, limits, xmlelement):
        addElement(xmlelement, 'AcceptableLimitRMS', limits[0])
        addElement(xmlelement, 'AcceptableLimit', limits[1])
        addElement(xmlelement, 'MaxDiff', limits[2])


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# - - - - - - - - -  - - - - - - - - -  - - - - - - - - - 
def addElement(containerXML, elementname, elementtext):
    #print 'addElement', elementname, type(elementtext), elementtext 
    e2 = ET.Element(elementname)
    e2.text = elementtext
    containerXML.append(e2)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

if __name__ == "__main__":
    logfile = sys.argv[1] 
    print(logfile)
    #mgr = makeGraphs(logfile)
    #grxml = mgr.getXML()
    
    logscrape = scaleitLogScraper()
    logscrape.process(logfile)

    print("\nData\n")
    xml = logscrape.makeXML()

    et = ET.ElementTree(xml)
    # and write out the XML
    ET.indent(et)
    et.write('scl.xml')

