from __future__ import print_function
from core.CCP4PluginScript import CPluginScript
from PySide2 import QtCore
import os,re,time,sys
from core import CCP4XtalData
#from lxml import etree
from xml.etree import ElementTree as ET
import math
from core import CCP4Modules,CCP4Utils
import base64

class ShelxCDEBase(CPluginScript):
    
    tagsAndKeys = {'N(data)':'NData','<I/sig>':'IOverSig','%Complete':'Completeness','<d"/sig>':'AnomalousSignal'}

    MAINTAINER = 'matin.noble@newcastle.co.uk'
    
    ERROR_CODES = {  200 : { 'description' : 'ShelxCD exited with error status' }, 201 : { 'description' : 'ShelxCD failed in mtz2various' },202 : { 'description' : ' ShelxCD failed in processOutputFiles' },203 : { 'description' : 'ShelxCD failed in f2mtz - Map file' },204 : { 'description' : 'ShelxCD failed in f2mtz - PHS file' },205 : { 'description' : 'ShelxCD failed in sftools - MAP file' },206 : { 'description' : 'ShelxCD failed in sftools - PHS file' }, 207 : { 'description' : 'Failed scraping Shelxd logfile' }, 208 : { 'description' : 'Failed scraping Shelxe logfile' }, 209 : { 'description' : 'Failed converting .hat to pdb' }, 210 : { 'description' : 'ShelxE ended with non zero status' }, 211 : { 'description' : 'ShelxE ended with non zero code' }, 212 : { 'description' : 'ShelxE BETA EXPIRED' }}
    
    PeriodicTable = [{'Name':'Hydrogen','Symbol':'H','Number':1},
                     {'Name':'Helium','Symbol':'He','Number':2},
                     {'Name':'Lithium','Symbol':'Li','Number':3},
                     {'Name':'Beryllium','Symbol':'Be','Number':4},
                     {'Name':'Boron','Symbol':'B','Number':5},
                     {'Name':'Carbon','Symbol':'C','Number':6},
                     {'Name':'Nitrogen','Symbol':'N','Number':7},
                     {'Name':'Oxygen','Symbol':'O','Number':8},
                     {'Name':'Fluorine','Symbol':'F','Number':9},
                     {'Name':'Neon','Symbol':'Ne','Number':10},
                     {'Name':'Sodium','Symbol':'Na','Number':11},
                     {'Name':'Magnesium','Symbol':'Mg','Number':12},
                     {'Name':'Aluminium','Symbol':'Al','Number':13},
                     {'Name':'Silicon','Symbol':'Si','Number':14},
                     {'Name':'Phosphorus','Symbol':'P','Number':15},
                     {'Name':'Sulfur','Symbol':'S','Number':16},
                     {'Name':'Chlorine','Symbol':'Cl','Number':17},
                     {'Name':'Argon','Symbol':'Ar','Number':18},
                     {'Name':'Potassium','Symbol':'K','Number':19},
                     {'Name':'Calcium','Symbol':'Ca','Number':20},
                     {'Name':'Scandium','Symbol':'Sc','Number':21},
                     {'Name':'Titanium','Symbol':'Ti','Number':22},
                     {'Name':'Vanadium','Symbol':'V','Number':23},
                     {'Name':'Chromium','Symbol':'Cr','Number':24},
                     {'Name':'Manganese','Symbol':'Mn','Number':25},
                     {'Name':'Iron','Symbol':'Fe','Number':26},
                     {'Name':'Cobalt','Symbol':'Co','Number':27},
                     {'Name':'Nickel','Symbol':'Ni','Number':28},
                     {'Name':'Copper','Symbol':'Cu','Number':29},
                     {'Name':'Zinc','Symbol':'Zn','Number':30},
                     {'Name':'Gallium','Symbol':'Ga','Number':31},
                     {'Name':'Germanium','Symbol':'Ge','Number':32},
                     {'Name':'Arsenic','Symbol':'As','Number':33},
                     {'Name':'Selenium','Symbol':'Se','Number':34},
                     {'Name':'Bromine','Symbol':'Br','Number':35},
                     {'Name':'Krypton','Symbol':'Kr','Number':36},
                     {'Name':'Rubidium','Symbol':'Rb','Number':37},
                     {'Name':'Strontium','Symbol':'Sr','Number':38},
                     {'Name':'Yttrium','Symbol':'Y','Number':39},
                     {'Name':'Zirconium','Symbol':'Zr','Number':40},
                     {'Name':'Niobium','Symbol':'Nb','Number':41},
                     {'Name':'Molybdenum','Symbol':'Mo','Number':42},
                     {'Name':'Technetium','Symbol':'Tc','Number':43},
                     {'Name':'Ruthenium','Symbol':'Ru','Number':44},
                     {'Name':'Rhodium','Symbol':'Rh','Number':45},
                     {'Name':'Palladium','Symbol':'Pd','Number':46},
                     {'Name':'Silver','Symbol':'Ag','Number':47},
                     {'Name':'Cadmium','Symbol':'Cd','Number':48},
                     {'Name':'Indium','Symbol':'In','Number':49},
                     {'Name':'Tin','Symbol':'Sn','Number':50},
                     {'Name':'Antimony','Symbol':'Sb','Number':51},
                     {'Name':'Tellurium','Symbol':'Te','Number':52},
                     {'Name':'Iodine','Symbol':'I','Number':53},
                     {'Name':'Xenon','Symbol':'Xe','Number':54},
                     {'Name':'Caesium','Symbol':'Cs','Number':55},
                     {'Name':'Barium','Symbol':'Ba','Number':56},
                     {'Name':'Lanthanum','Symbol':'La','Number':57},
                     {'Name':'Cerium','Symbol':'Ce','Number':58},
                     {'Name':'Praseodymium','Symbol':'Pr','Number':59},
                     {'Name':'Neodymium','Symbol':'Nd','Number':60},
                     {'Name':'Promethium','Symbol':'Pm','Number':61},
                     {'Name':'Samarium','Symbol':'Sm','Number':62},
                     {'Name':'Europium','Symbol':'Eu','Number':63},
                     {'Name':'Gadolinium','Symbol':'Gd','Number':64},
                     {'Name':'Terbium','Symbol':'Tb','Number':65},
                     {'Name':'Dysprosium','Symbol':'Dy','Number':66},
                     {'Name':'Holmium','Symbol':'Ho','Number':67},
                     {'Name':'Erbium','Symbol':'Er','Number':68},
                     {'Name':'Thulium','Symbol':'Tm','Number':69},
                     {'Name':'Ytterbium','Symbol':'Yb','Number':70},
                     {'Name':'Lutetium','Symbol':'Lu','Number':71},
                     {'Name':'Hafnium','Symbol':'Hf','Number':72},
                     {'Name':'Tantalum','Symbol':'Ta','Number':73},
                     {'Name':'Tungsten','Symbol':'W','Number':74},
                     {'Name':'Rhenium','Symbol':'Re','Number':75},
                     {'Name':'Osmium','Symbol':'Os','Number':76},
                     {'Name':'Iridium','Symbol':'Ir','Number':77},
                     {'Name':'Platinum','Symbol':'Pt','Number':78},
                     {'Name':'Gold','Symbol':'Au','Number':79},
                     {'Name':'Mercury','Symbol':'Hg','Number':80},
                     {'Name':'Thallium','Symbol':'Tl','Number':81},
                     {'Name':'Lead','Symbol':'Pb','Number':82},
                     {'Name':'Bismuth','Symbol':'Bi','Number':83},
                     {'Name':'Polonium','Symbol':'Po','Number':84},
                     {'Name':'Astatine','Symbol':'At','Number':85},
                     {'Name':'Radon','Symbol':'Rn','Number':86},
                     {'Name':'Francium','Symbol':'Fr','Number':87},
                     {'Name':'Radium','Symbol':'Ra','Number':88},
                     {'Name':'Actinium','Symbol':'Ac','Number':89},
                     {'Name':'Thorium','Symbol':'Th','Number':90},
                     {'Name':'Protactinium','Symbol':'Pa','Number':91},
                     {'Name':'Uranium','Symbol':'U','Number':92},
                     {'Name':'Neptunium','Symbol':'Np','Number':93},
                     {'Name':'Plutonium','Symbol':'Pu','Number':94},
                     {'Name':'Americium','Symbol':'Am','Number':95},
                     {'Name':'Curium','Symbol':'Cm','Number':96},
                     {'Name':'Berkelium','Symbol':'Bk','Number':97},
                     {'Name':'Californium','Symbol':'Cf','Number':98},
                     {'Name':'Einsteinium','Symbol':'Es','Number':99},
                     {'Name':'Fermium','Symbol':'Fm','Number':100},
                     {'Name':'Mendelevium','Symbol':'Md','Number':101},
                     {'Name':'Nobelium','Symbol':'No','Number':102},
                     {'Name':'Lawrencium','Symbol':'Lr','Number':103},
                     {'Name':'Rutherfordium','Symbol':'Rf','Number':104},
                     {'Name':'Dubnium','Symbol':'Db','Number':105},
                     {'Name':'Seaborgium','Symbol':'Sg','Number':106},
                     {'Name':'Bohrium','Symbol':'Bh','Number':107},
                     {'Name':'Hassium','Symbol':'Hs','Number':108},
                     {'Name':'Meitnerium','Symbol':'Mt','Number':109},
                     {'Name':'Darmstadtium','Symbol':'Ds','Number':110},
                     {'Name':'Roentgenium','Symbol':'Rg','Number':111},
                     {'Name':'Copernicium','Symbol':'Cp','Number':112},
                     {'Name':'Ununtrium','Symbol':'Uut','Number':113},
                     {'Name':'Ununquadium','Symbol':'Uuq','Number':114},
                     {'Name':'Ununpentium','Symbol':'Uup','Number':115},
                     {'Name':'Ununhexium','Symbol':'Uuh','Number':116},
                     {'Name':'Ununseptium','Symbol':'Uus','Number':117},
                     {'Name':'Ununoctium','Symbol':'Uuo','Number':118}]
    
    def prepareShelxc(self):
        
        #Test for location of shelxc executable:
        bin =  self.getCommand('shelxc')

        datasets = {}
        inp = self.container.inputData
        if self.container.controlParameters.MODE.__str__() == 'SAD':
            datasets = {'SAD':inp.SAD,'NAT':inp.NAT}
        elif self.container.controlParameters.MODE.__str__() == 'MAD':
            datasets = {'HREM':inp.HREM,'LREM':inp.LREM,'INFL':inp.INFL,'PEAK':inp.PEAK,'NAT':inp.NAT}
        elif self.container.controlParameters.MODE.__str__() == 'SIR':
            datasets = {'SIR':inp.SIR,'NAT':inp.NAT}
        elif self.container.controlParameters.MODE.__str__() == 'SIRAS':
            datasets = {'SIRA':inp.SIRA,'NAT':inp.NAT}
        elif self.container.controlParameters.MODE.__str__() == 'RIP':
            datasets = {'RIP':inp.RIP,'NAT':inp.NAT}
        elif self.container.controlParameters.MODE.__str__() == 'RIPAS':
            datasets = {'RIPA':inp.RIPA,'NAT':inp.NAT}
        
        for datasetName in datasets:
            from core.CCP4XtalData import CObsDataFile
            if datasets[datasetName].isSet():
                self.makeHklin([datasetName])
                mtzFilepath = os.path.join(self.getWorkDirectory(),'hklin.mtz')
                if datasets[datasetName].contentFlag == CObsDataFile.CONTENT_FLAG_IPAIR:
                    hklFilepath = os.path.join(self.getWorkDirectory(),'result_'+datasetName+'.sca')
                    result = self.makeScafile(infile=mtzFilepath, outfile=[hklFilepath,["Iplus","SIGIplus","Iminus","SIGIminus"]])
                elif datasets[datasetName].contentFlag == CObsDataFile.CONTENT_FLAG_FPAIR:
                    hklFilepath = os.path.join(self.getWorkDirectory(),'result_'+datasetName+'.sca')
                    result = self.makeScafile(infile=mtzFilepath, outfile=[hklFilepath,["Fplus","SIGFplus","Fminus","SIGFminus"]])
                elif datasets[datasetName].contentFlag == CObsDataFile.CONTENT_FLAG_IMEAN:
                    hklFilepath = os.path.join(self.getWorkDirectory(),'result_'+datasetName+'.hkl')
                    result = self.makeHklfile(infile=mtzFilepath, outfile=[hklFilepath,["I=I","SIGI=SIGI"]])
                elif datasets[datasetName].contentFlag == CObsDataFile.CONTENT_FLAG_FMEAN:
                    hklFilepath = os.path.join(self.getWorkDirectory(),'result_'+datasetName+'.hkl')
                    result = self.makeHklfile(infile=mtzFilepath, outfile=[hklFilepath,["FP=SIGF","SIGFP=SIGF"]])
                if result != CPluginScript.SUCCEEDED: return result
        return CPluginScript.SUCCEEDED

    def makeScafile( self, infile=None, outfile=[], parentPlugin=None ):
        bin =  os.path.normpath(os.path.join( CCP4Utils.getCCP4Dir(), 'bin', 'mtz2sca' ))
        name,colin = outfile
        logFile =  os.path.normpath(os.path.join(self.getWorkDirectory(), 'mtz2sca.log'))
        arglist = ['-p',colin[0],'-P',colin[1],'-m',colin[2],'-M',colin[3],infile,name]
        pid = CCP4Modules.PROCESSMANAGER().startProcess(bin,arglist,logFile=logFile)
        status = CCP4Modules.PROCESSMANAGER().getJobData(pid)
        exitCode = CCP4Modules.PROCESSMANAGER().getJobData(pid,'exitCode')
        if status == 0 and os.path.exists(outfile[0]):
            return CPluginScript.SUCCEEDED
        else:
            self.appendErrorReport(201)
            return CPluginScript.FAILED

    def makeHklfile( self, infile=None, outfile=[], parentPlugin=None ):
        bin =  os.path.normpath(os.path.join( CCP4Utils.getCCP4Dir(), 'bin', 'mtz2various' ))
        name,colin = outfile
        logFile =  os.path.normpath(os.path.join(self.getWorkDirectory(), 'mtz2various.log'))
        arglist = ['hklin',infile,'hklout',name]
        inputText = 'output SHELX\n'
        inputText += ('LABIN ' + ' '.join(colin) + '\n')
        pid = CCP4Modules.PROCESSMANAGER().startProcess(bin,arglist,inputText=inputText, logFile=logFile)
        status = CCP4Modules.PROCESSMANAGER().getJobData(pid)
        exitCode = CCP4Modules.PROCESSMANAGER().getJobData(pid,'exitCode')
        if status == 0 and os.path.exists(outfile[0]):
            return CPluginScript.SUCCEEDED
        else:
            self.appendErrorReport(201)
            return CPluginScript.FAILED

    def runShelxc( self ):
        # CpluginScript.getCommand checks preferences for shelxc directory
        bin =  self.getCommand('shelxc')
        logFile =  os.path.normpath(os.path.join(self.getWorkDirectory(), 'shelxc.log'))
        arglist = ['result']

        inputText = ''
        cell = self.container.inputData.SAD.fileContent.cell
        factor = 1.#180. / math.pi
        cellString, spaceGroupNumber, extantFile = self.cellString()
        inputText += cellString
        #From shelx docs

        # SPAG - the name of the space group. Only Sohncke space groups are permitted, but some common non-standard settings are allowed,
        #  e.g. 'P22121'. Embedded spaces are ignored. SPAG is used to generate the LATT and SYMM instructions that are written to name_fa.ins.
        #If the space group is specified as 'R3' or 'R32' the program checks the cell dimensions to see whether the hexagonal or
        #    primitive rhombohedral setting is required.

        inputText += 'SPAG %s\n'%extantFile.fileContent.spaceGroup.__str__().replace(' ','').replace(':H','').replace(':R','')
        try:
            inputText += 'SFAC %s\n'%self.container.controlParameters.SFAC
        except:
            print('#ShelxCDEBase runShelxc: SFAC unset or not relevant')
        
        #Simple float parameters
        for parameterName in ['RIPW','ASCA','DSCA']:
            if parameterName in self.container.controlParameters.contents():
                parameter = getattr(self.container.controlParameters, parameterName)
                if parameter.isSet():
                    inputText += (parameterName+' %8.2f\n'%(float(parameter)))
                              
        #Float range parameters
        for parameterName in ['MIND','SHEL']:
            if parameterName in self.container.controlParameters.contents():
                parameter = getattr(self.container.controlParameters, parameterName)
                if parameter.isSet():
                    inputText += (parameterName+' %8.2f %8.2f\n'%(float(parameter.start), float(parameter.end )))

        datasetNames = []
        if self.container.controlParameters.MODE.__str__() == 'SAD':
            datasetNames = ['SAD', 'NAT']
        elif self.container.controlParameters.MODE.__str__() == 'MAD':
            datasetNames = ['HREM','LREM','PEAK','INFL', 'NAT']
        elif self.container.controlParameters.MODE.__str__() == 'SIR':
            datasetNames = ['SIR','NAT']
        elif self.container.controlParameters.MODE.__str__() == 'SIRAS':
            datasetNames = ['SIRA','NAT']
        elif self.container.controlParameters.MODE.__str__() == 'RIP':
            datasetNames = ['RIP','NAT']
        elif self.container.controlParameters.MODE.__str__() == 'RIPAS':
            datasetNames = ['RIPA','NAT']

        from core.CCP4XtalData import CObsDataFile
        inp = self.container.inputData
        for datasetName in datasetNames:
            dataset = getattr(inp,datasetName)
            if dataset.isSet():
                if dataset.contentFlag == CObsDataFile.CONTENT_FLAG_IMEAN:
                    inputText += (datasetName+' result_'+datasetName+'.hkl\n')
                elif dataset.contentFlag == CObsDataFile.CONTENT_FLAG_FMEAN:
                    inputText += (datasetName+' result_'+datasetName+'.hkl\n')
                else:
                    inputText += (datasetName+' result_'+datasetName+'.sca\n')


        try:
            inputText += 'FIND %d\n'%int(self.container.controlParameters.FIND)
        except:
            print('#ShelxCDEBase runShelxc: FIND unset or not relevant')
        
        try:
            if self.container.controlParameters.NTRY.isSet():
                inputText += 'NTRY %d\n'%int(self.container.controlParameters.NTRY)
        except:
            print('#ShelxCDEBase runShelxc: NTRY unset or not relevant')

        pid = CCP4Modules.PROCESSMANAGER().startProcess(bin, arglist, inputText=inputText, logFile=logFile,cwd=self.getWorkDirectory())
        status = CCP4Modules.PROCESSMANAGER().getJobData(pid)
        exitCode = CCP4Modules.PROCESSMANAGER().getJobData(pid,'exitCode')

        if status == 0 and os.path.exists(os.path.join(self.getWorkDirectory(),'result_fa.hkl')):
            return CPluginScript.SUCCEEDED
        else:
            self.appendErrorReport(201,str(bin)+' '+str(arglist)+' '+str(inputText)+' '+str(logFile))
            return CPluginScript.FAILED

    def scrapeShelxcLog(self, xmlroot, logFile):
        shelxcNode = ET.SubElement(xmlroot, 'Shelxc')
        analysisNode = None
        with open(logFile,'r') as insFile:
            insLines = [line.strip() for line in insFile.readlines()]
            datasetNode = None
            e2TableLine = -1
            for line in insLines:
                tokens = line.split()
                if 'Reflections read from' in line:
                    datasetNode = ET.SubElement(shelxcNode,'Dataset')
                    self.subElementOfTypeWithText(datasetNode,'Name', tokens[4])
                    self.subElementOfTypeWithText(datasetNode,'ReflectionsRead', tokens[0])
                elif 'Unique reflections, highest resolution' in line:
                    self.subElementOfTypeWithText(datasetNode,'UniqueReflections', tokens[0])
                    self.subElementOfTypeWithText(datasetNode,'HighestResolution', tokens[5])
                elif 'for input to SHELXD/E' in line:
                    self.subElementOfTypeWithText(shelxcNode,'FAsWritten', tokens[0])
                elif 'for input to SHELXE' in line:
                    self.subElementOfTypeWithText(shelxcNode,'FsWritten', tokens[0])
                elif len(tokens) == 6 and 'Unique reflections, highest resolution' in line and datasetNode is not None:
                    self.subElementOfTypeWithText(datasetNode,'UniqueReflections', tokens[0])
                    self.subElementOfTypeWithText(datasetNode,'HighestResolutionRead', tokens[5])
                elif len(tokens) == 8 and 'Friedel pairs used on average for local' in line and datasetNode is not None:
                    self.subElementOfTypeWithText(datasetNode,'AverageFriedelsPerBin', tokens[0])
                elif len(tokens) > 0 and tokens[0] == 'Resl.' and datasetNode is not None:
                    analysisNode = ET.SubElement(datasetNode,'DataAnalysis')
                    for iBin in range (1, len(tokens)-1):
                        binNode = ET.SubElement(analysisNode,'Bin')
                        self.subElementOfTypeWithText(binNode,'LowRes', tokens[iBin])
                        self.subElementOfTypeWithText(binNode,'HighRes', tokens[iBin+1])
                elif len(tokens) > 0 and tokens[0] in ShelxCDEBase.tagsAndKeys and analysisNode is not None:
                    for iBin in range (1, len(tokens)):
                        binNode = analysisNode[iBin-1]
                        self.subElementOfTypeWithText(binNode,ShelxCDEBase.tagsAndKeys[tokens[0]], tokens[iBin])
                elif 'Correlation coefficients' in line:
                    datasetNode = None
                elif e2TableLine >= 0:
                    if len(tokens) == 0: e2TableLine = -1
                    else:
                        datasetNodes = shelxcNode.findall('.//Dataset')
                        if e2TableLine < len(datasetNodes):
                            self.subElementOfTypeWithText(datasetNodes[e2TableLine],'eSquaredMinus1',tokens[1])
                            e2TableLine += 1
                elif len(tokens) > 1 and tokens[0] == 'Dataset' and tokens[1] == '<|E^2-1|>':
                    e2TableLine = 0

    def handleShelxeLogChanged(self, logFilename):
        import datetime
        timeNow = datetime.datetime.now()
        if not hasattr(self,"lastScrapeTime"): self.lastScrapeTime = timeNow
        deltaTime = timeNow - self.lastScrapeTime
        if deltaTime.seconds > 0 or deltaTime.days > 0:
            self.scrapeShelxeLog(self.xmlroot)
            self.flushXML()
            self.lastScrapeTime = timeNow
    
    def handleShelxdLogChanged(self, logFilename):
        import datetime
        timeNow = datetime.datetime.now()
        if not hasattr(self,"lastScrapeTime"): self.lastScrapeTime = timeNow
        deltaTime = timeNow - self.lastScrapeTime
        if deltaTime.seconds > 0 or deltaTime.days > 0:
            self.scrapeShelxdLog(self.xmlroot)
            self.flushXML()
            self.lastScrapeTime = timeNow

    def subElementOfTypeWithText(self, parent, type, text):
        result = ET.SubElement(parent,type)
        result.text = str(text)
        return result

    def flushXML(self):
        tmpFileName = self.makeFileName('PROGRAMXML') + '_tmp'
        with open(tmpFileName,'w') as xmlFile:
            ET.indent(self.xmlroot)
            CCP4Utils.writeXML(xmlFile,ET.tostring(self.xmlroot))
        #This is attempt to fix error overwriting exisiting file on Windows
        PROGRAMXML = self.makeFileName('PROGRAMXML')
        #print 'flushXML',PROGRAMXML
        if os.path.exists(PROGRAMXML): os.remove(PROGRAMXML)
        self.renameFile(tmpFileName, PROGRAMXML)

    def pdbToRes( self, pdbDataFile, insFilePath):
        bin =  os.path.normpath(os.path.join( CCP4Utils.getCCP4Dir(), 'bin', 'coordconv' ))
        logFile =  os.path.normpath(os.path.join(self.getWorkDirectory(), 'coordconv.log'))
        fracFilePath = os.path.normpath(os.path.join(self.getWorkDirectory(), 'ha.frac'))
        arglist = ['XYZIN',pdbDataFile.fullPath.__str__(),'XYZOUT',fracFilePath]
        
        inputText = '''INPUT pdb
OUTPUT frac
'''
        pid = CCP4Modules.PROCESSMANAGER().startProcess(bin, arglist, inputText=inputText, logFile=logFile,cwd=self.getWorkDirectory())
        status = CCP4Modules.PROCESSMANAGER().getJobData(pid)
        exitCode = CCP4Modules.PROCESSMANAGER().getJobData(pid,'exitCode')

        if status != 0 or not os.path.exists(fracFilePath):
            self.appendErrorReport(201,str(bin)+' '+str(arglist)+' '+str(inputText)+' '+str(fracFilePath))
            return CPluginScript.FAILED

        linetypesToGrab = ['TITL','CELL','LATT','SYMM','SFAC','UNIT','HKLF','END']
        grabbedLines = {}

        with open(insFilePath,'r') as insFile:
            insLines = insFile.readlines()
            for linetype in linetypesToGrab:
                grabbedLines[linetype] = [insLine for insLine in insLines if insLine.startswith(linetype)]

        coords = []
        with open(fracFilePath,'r') as fracFile:
            fracLines = fracFile.readlines()
            for fracLine in fracLines:
                coord = {'NUMBER':fracLine[0:5].strip(), 'X':float(fracLine[5:15].strip()), 'Y':float(fracLine[15:25].strip()), 'Z':float(fracLine[25:35].strip()), 'BFAC':float(fracLine[35:45].strip()), 'OCC':float(fracLine[45:50].strip()), 'AtomicNr':int(fracLine[50:55].strip())}

                coords.append(coord)

        resFilePath = os.path.normpath(os.path.join(self.getWorkDirectory(),'result_fa.res'))
        
        # Get a list of atom cards to identify elements which should be on SFAC line
        atomLines = []
        oldElements = []
        if 'SFAC' in grabbedLines and len(grabbedLines['SFAC'][0].split()) > 1:
            oldElements = grabbedLines['SFAC'][0].split()[1:]
        elements = [element.upper() for element in oldElements]
        iAtom = 1
        for coord in coords:
            #MN CHECKME FIXME The following line forces unknown elements to be Seleniums
            if coord['AtomicNr'] == 0: coord['AtomicNr'] = 34
            element = [ element for element in ShelxCDEBase.PeriodicTable if element['Number'] == coord['AtomicNr']]
            elementSymbol = element[0]['Symbol'].upper()
            if elementSymbol not in elements: elements.append(elementSymbol)
            atomLine = '%s%d %d %10.6f%10.6f%10.6f%10.6f%10.6f\n'%(elementSymbol,iAtom,1,coord['X'],coord['Y'],coord['Z'],coord['OCC'],coord['BFAC']/100.)
            atomLines.append(atomLine)
            iAtom += 1
        
        with open(resFilePath, 'w') as resFile:
            for linetype in ['TITL','CELL','LATT','SYMM','UNIT']:
                if linetype in grabbedLines:
                    for grabbedLine in grabbedLines[linetype]:
                        resFile.write(grabbedLine)
            resFile.write('SFAC '+' '.join(elements)+'\n')
            iAtom = 1
            for atomLine in atomLines:
                resFile.write(atomLine)
            for linetype in ['HKLF','END']:
                if linetype in grabbedLines:
                    for grabbedLine in grabbedLines[linetype]:
                        resFile.write(grabbedLine)
                            
        return CPluginScript.SUCCEEDED

    def cellString(self):
        fileTypes = ['NAT','SAD','HREM','LREM','INFL','PEAK','SIR','SIRA','RIP','RIPA']
        extantFile = None
        
        for fileType in fileTypes:
            fileObj = getattr(self.container.inputData, fileType, None)
            if fileObj is not None and fileObj.isSet() and os.path.isfile(fileObj.fullPath.__str__()):
                extantFile  = fileObj
            if extantFile is not None: break

        if extantFile is None:
            return CPluginScript.FAILED
        
        cell = extantFile.fileContent.cell
        factor = 1.#180. / math.pi
        cellString = 'CELL %.3f %.3f %.3f %.3f %.3f %.3f\n' %(cell.a, cell.b, cell.c, factor*float(cell.alpha), factor*float(cell.beta), factor*float(cell.gamma),)
        spaceGroupNumber = extantFile.fileContent.spaceGroup.number()
        
        return cellString, spaceGroupNumber, extantFile

    def harvestPhsFile(self, phsFilePath, mapFilePath, phsOutFilePath=None):
        
        bin =  os.path.normpath(os.path.join( CCP4Utils.getCCP4Dir(), 'bin', 'f2mtz' ))
        logFile =  os.path.normpath(os.path.join(self.getWorkDirectory(), 'f2mtz_map.log'))
        tmpMtzFilePath = mapFilePath+'_tmp'
        arglist = ['HKLIN',phsFilePath,'HKLOUT',tmpMtzFilePath]
        
        inputText = "FORMAT '(F4.0,F4.0,F4.0,F9.2,F8.4,F8.1)'\n"
        inputText += "LABOUT H K L FOBS FOM PHI\n"
        inputText += "CTYPOUT H H H F W P\n"
        cellString, spaceGroupNumber, extantFile = self.cellString()
        inputText += cellString
        spgrName = extantFile.fileContent.spaceGroup.__str__().replace(' ','')
        cellTokens = cellString.split()
        if spgrName == 'H32' and cellTokens[1] == cellTokens[2] and cellTokens[1] == cellTokens[3] and cellTokens[4] == cellTokens[5] and cellTokens[4] == cellTokens[6]: spgrName = 'R32'
        elif spgrName == 'R32' and cellTokens[1] == cellTokens[2] and cellTokens[4] == '90.000' and cellTokens[6] == '120.000': spgrName = 'H32'

        inputText += 'SYMMETRY %s\n'%spgrName

        pid = CCP4Modules.PROCESSMANAGER().startProcess(bin, arglist, inputText=inputText, logFile=logFile,cwd=self.getWorkDirectory())
        status = CCP4Modules.PROCESSMANAGER().getJobData(pid)
        exitCode = CCP4Modules.PROCESSMANAGER().getJobData(pid,'exitCode')

        if status != 0 or not os.path.exists(tmpMtzFilePath):
            self.appendErrorReport(203,str(bin)+' '+str(arglist)+' '+str(inputText)+' '+str(logFile))
            return CPluginScript.FAILED
        
        #Separate out mFo,Phi map coefficents
        bin =  os.path.normpath(os.path.join( CCP4Utils.getCCP4Dir(), 'bin', 'sftools' ))
        logFile =  os.path.normpath(os.path.join(self.getWorkDirectory(), 'sftools_map.log'))
        argList = []
        
        inputText = "READ "+tmpMtzFilePath+"\n"
        inputText += "REDUCE\n"
        inputText += "SORT\n"
        inputText += "CALC COL F = COL 1 COL 2 *\n"
        inputText += "SET TYPE COL 4\n"
        inputText += "F\n"
        inputText += "WRITE "+mapFilePath+" COL 4 3\n"
        inputText += "STOP\n"

        pid = CCP4Modules.PROCESSMANAGER().startProcess(bin, arglist, inputText=inputText, logFile=logFile,cwd=self.getWorkDirectory())
        status = CCP4Modules.PROCESSMANAGER().getJobData(pid)
        exitCode = CCP4Modules.PROCESSMANAGER().getJobData(pid,'exitCode')

        if status != 0 or not os.path.exists(mapFilePath):
            self.appendErrorReport(205,str(bin)+' '+str(arglist)+' '+str(inputText)+' '+str(logFile))
            return CPluginScript.FAILED
                
        #Separate out PHI,FOM Phase probability distribution
        if phsOutFilePath is not None:
            bin =  os.path.normpath(os.path.join( CCP4Utils.getCCP4Dir(), 'bin', 'sftools' ))
            logFile =  os.path.normpath(os.path.join(self.getWorkDirectory(), 'sftools_phs.log'))
            argList = []
            
            inputText = "READ "+tmpMtzFilePath+"\n"
            inputText += "REDUCE\n"
            inputText += "SORT\n"
            inputText += "WRITE "+phsOutFilePath+" COL 3 2\n"
            inputText += "STOP\n"

            pid = CCP4Modules.PROCESSMANAGER().startProcess(bin, arglist, inputText=inputText, logFile=logFile,cwd=self.getWorkDirectory())
            status = CCP4Modules.PROCESSMANAGER().getJobData(pid)
            exitCode = CCP4Modules.PROCESSMANAGER().getJobData(pid,'exitCode')

            if status != 0 or not os.path.exists(phsOutFilePath):
                self.appendErrorReport(206,str(bin)+' '+str(arglist)+' '+str(inputText)+' '+str(logFile))
                return CPluginScript.FAILED

        return CPluginScript.SUCCEEDED

    def scrapeShelxdLog(self, xmlroot):
        if True:
            oldShelxdNodes = xmlroot.findall('Shelxd')
            for oldShelxNode in oldShelxdNodes:
                xmlroot.remove(oldShelxNode)
            shelxdNode = ET.SubElement(xmlroot,'Shelxd')
            with open(self.makeFileName('LOG')) as logFile:
                longlines = logFile.readlines()
                lines = [line.strip() for line in longlines]
                for line in lines:
                    #Good fun here: the number of the CPU was fixed format, so n> 10 bumped
                    #into the word CPU :-)
                    tokens = line.replace('CPU','CPU ').split()
                    if len(tokens) > 0 and tokens[0] == 'Try':
                        # Try     15, CPU 4, CC All/Weak 37.5 / 22.8, CFOM 60.3, best 60.7, PATFOM  45.67
                        tryNode = ET.SubElement(shelxdNode,'Try')
                        self.subElementOfTypeWithText(tryNode,'Number',tokens[1][:-1])
                        self.subElementOfTypeWithText(tryNode,'CCAll',tokens[6])
                        self.subElementOfTypeWithText(tryNode,'CCWeak',tokens[8][:-1])
                        self.subElementOfTypeWithText(tryNode,'CFOM',tokens[10][:-1])
                        try:
                            self.subElementOfTypeWithText(tryNode,'PATFOM',tokens[14])
                        except:
                            print(tokens)
        
            lstPath = os.path.normpath(os.path.join(self.getWorkDirectory(),'result_fa.lst'))
            if not os.path.isfile(lstPath): lstPath = os.path.normpath(os.path.join(self.getWorkDirectory(),'result_fa.lst.txt'))
            if os.path.isfile(lstPath):
                with open(lstPath,'rb') as lstFile:
                    lstText = lstFile.read()
                    lstNode = ET.SubElement(shelxdNode,'LstText')
                    #lstNode.text=etree.CDATA(lstText)
                    lstNode.text=base64.b64encode(lstText).decode()
            return CPluginScript.SUCCEEDED
        else:
            self.appendErrorReport(207,self.getWorkDirectory())
            return CPluginScript.FAILED
    

    def scrapeShelxeLog(self, xmlroot):
        try:
            oldShelxeNodes = xmlroot.findall('Shelxe')
            for oldShelxNode in oldShelxeNodes:
                xmlroot.remove(oldShelxNode)
            shelxeNode = ET.SubElement(xmlroot,'Shelxe')
            globalTraceNode = ET.SubElement(shelxeNode,'GlobalAutotracingCycle')
            numberNode = self.subElementOfTypeWithText(globalTraceNode,'Number','1')
            
            with open(self.makeFileName('LOG')) as logFile:
                longlines = logFile.readlines()
                lines = [line.strip() for line in longlines]
                for line in lines:
                    tokens = line.split()
                    if len(tokens) > 0 and 'for dens.mod.' in line:
                        # Try     15, CPU 4, CC All/Weak 37.5 / 22.8, CFOM 60.3, best 60.7, PATFOM  45.67
                        cycleNode = ET.SubElement(globalTraceNode,'Cycle')
                        self.subElementOfTypeWithText(cycleNode,'Wt',tokens[2][:-1])
                        self.subElementOfTypeWithText(cycleNode,'Contrast',tokens[5][:-1])
                        self.subElementOfTypeWithText(cycleNode,'Connect',tokens[8])
                        self.subElementOfTypeWithText(cycleNode,'Number',tokens[12])
                    elif 'Global autotracing cycle' in line:
                        globalTraceNode = ET.SubElement(shelxeNode,'GlobalAutotracingCycle')
                        numberNode = self.subElementOfTypeWithText(globalTraceNode,'Number',tokens[3])
                    elif 'This beta-test version has expired - please update it' in line:
                        betaExpiredNode = ET.SubElement(shelxeNode,'BETAEXPIRED')
                    elif 'Estimated mean FOM =' in line:
                        #Estimated mean FOM = 0.494   Pseudo-free CC = 52.89 %
                        overallFOMsNode = ET.SubElement(globalTraceNode,'OverallFOMs')
                        self.subElementOfTypeWithText(overallFOMsNode,'FOM',tokens[4])
                        self.subElementOfTypeWithText(overallFOMsNode,'PseudoFreeCC',tokens[8])
    
            lstPath = os.path.normpath(os.path.join(self.getWorkDirectory(),'result.lst'))
            for lstRoot in ['result.lst','result_i.lst','result.lst.txt','result_i.lst.txt']:
                lstPath = os.path.normpath(os.path.join(self.getWorkDirectory(),lstRoot))
                if os.path.isfile(lstPath):
                    with open(lstPath,'r') as lstFile:
                        lstText = lstFile.read()
                        lstNode = ET.SubElement(shelxeNode,'LstText')
                        #lstNode.text=etree.CDATA(lstText)
                        lstNode.text=base64.b64encode(lstText)
            return CPluginScript.SUCCEEDED
        except:
            self.appendErrorReport(208,self.getWorkDirectory())
            return CPluginScript.FAILED

    def invertCoordinates(self, inputPDBPath, outputPDBPath):
        inversion_operation = '-X,-Y,-Z'
        cellString, spaceGroupNumber, extantFile = self.cellString()
        if spaceGroupNumber==80:
            inversion_operation = '-X+1/2,-Y,-Z'
        elif spaceGroupNumber==98:
            inversion_operation = '-X+1/2,-Y,-Z+1/4'
        elif spaceGroupNumber==210:
            inversion_operation = '-X+1/4,-Y+1/4,-Z+1/4'
    
        bin =  os.path.normpath(os.path.join( 'pdbcur' ))
        logFile =  os.path.normpath(os.path.join(self.getWorkDirectory(), 'pdbcur.log'))
        arglist = ['XYZIN',inputPDBPath,'XYZOUT',outputPDBPath]
        
        inputText = 'SYMOP '+inversion_operation+'\n'
        inputText += 'symcommit\n'
        
        pid = CCP4Modules.PROCESSMANAGER().startProcess(bin, arglist, inputText=inputText, logFile=logFile,cwd=self.getWorkDirectory())
        status = CCP4Modules.PROCESSMANAGER().getJobData(pid)
        exitCode = CCP4Modules.PROCESSMANAGER().getJobData(pid,'exitCode')

        if status != 0 or not os.path.exists(outputPDBPath):
            self.appendErrorReport(201,str(bin)+' '+str(arglist)+' '+str(inputText)+' '+str(inputPDBPath))
            return CPluginScript.FAILED

        return CPluginScript.SUCCEEDED

    def resToPDB(self, inputResPath, outputPDBPath):
        elements = []
        coords = []
        with open(inputResPath,'r') as inputRes:
            lines = [line.strip() for line in inputRes.readlines()]
            for line in lines:
                tokens = line.split()
                if len(tokens) > 0 and tokens[0] == 'SFAC':
                    elements = [token.upper() for token in tokens[1:]]
                elif len(tokens) == 7:
                    coordElement = elements[int(tokens[1])-1]
                    coord = {'Element':coordElement,'X':float(tokens[2]),'Y':float(tokens[3]),'Z':float(tokens[4]),'OCC':float(tokens[5]),'B':float(tokens[6])}
                    #coord = {'element':coordElement,'xFrac':float(tokens[2]),'yFrac':float(tokens[3]),'zFrac':float(tokens[4]),'occupancy':float(tokens[5]),'tempFactor':float(tokens[6])}
                    coords.append(coord)

        """
        # This code (with the alternate definition of coord above) is better way to create PDB file
        from core import CCP4ModelData
        pdb = CCP4ModelData.CPdbData()
        cellString, spaceGroupNumber, extantFile = self.cellString()
        pdb.makeOneResPdb(atomDefList=coords,cell=extantFile.fileContent.cell,spaceGroup=extantFile.fileContent.spaceGroup,fileName=outputPDBPath)
        """
        
        outputFracPath = outputPDBPath+'.ha'
        with open(outputFracPath,'w') as outputHA:
            cellString, spaceGroupNumber, extantFile = self.cellString()
            outputHA.write(cellString)
            iCoord = 1
            for coord in coords:
                element = [ element for element in ShelxCDEBase.PeriodicTable if element['Symbol'].upper() == coord['Element']][0]
                #I confess I don't understand negative occupancies here
                outputHA.write ('ATOM %s %10.5f %10.5f %10.5f %10.5f %10.5f\n'%(element['Symbol'], coord['X'], coord['Y'], coord['Z'], coord['OCC'], coord['B']*100.))

        bin =  os.path.normpath(os.path.join( 'coordconv' ))
        logFile =  os.path.normpath(os.path.join(self.getWorkDirectory(), 'coordconv_2.log'))
        arglist = ['XYZIN',outputFracPath,'XYZOUT',outputPDBPath]
        
        inputText = '''INPUT HA
OUTPUT pdb
'''
        inputText += self.cellString()[0]
        pid = CCP4Modules.PROCESSMANAGER().startProcess(bin, arglist, inputText=inputText, logFile=logFile,cwd=self.getWorkDirectory())
        status = CCP4Modules.PROCESSMANAGER().getJobData(pid)
        exitCode = CCP4Modules.PROCESSMANAGER().getJobData(pid,'exitCode')

        if status != 0 or not os.path.exists(outputPDBPath):
            self.appendErrorReport(201,str(bin)+' '+str(arglist)+' '+str(inputText)+' '+str(outputFracPath))
            return CPluginScript.FAILED
        

        return CPluginScript.SUCCEEDED



    def txtOutputFiles(self):
        # Add '.txt' extension to files so will display correctly in i2 and desktop tools
        import glob,shutil
        # also : pha phs hat from shelxe
        for ext in [ '*.lst', '*.res', '*.ins', '*.frac' ]:
          fileList = glob.glob(os.path.join(self.getWorkDirectory(),ext))
          for fName in fileList:
            try:
              shutil.move(fName,fName+'.txt')
            except:
              print('Error moving',str(fName),'to',str(fName+'.txt'))
