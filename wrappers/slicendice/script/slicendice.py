from __future__ import print_function
import os
import shutil
import json
import multiprocessing
#from lxml import etree
from xml.etree import ElementTree as ET
from core import CCP4Utils
from core import CCP4XtalData
from core import CCP4File
from core.CCP4PluginScript import CPluginScript
from core.CCP4Modules import PROCESSMANAGER
from core import CCP4ErrorHandling
from core.CCP4ErrorHandling import *

class slicendice(CPluginScript):

    TASKTITLE='SliceNDice'
    TASKNAME = 'slicendice'
    TASKMODULE= ['alpha_fold']
    TASKCOMMAND = 'slicendice'
    TASKVERSION= 0.1
    MAINTAINER = 'ronan.keegan@stfc.ac.uk'
    PERFORMANCECLASS = 'CRefinementPerformance'
    
    ERROR_CODES = {19121 : {'description' : 'SliceNDice, Json Data file not found. ' \
                                                                           'Please check the SliceNDice log file for details.'},
                   19122 : {'description' : 'SliceNDice, No solution found in json file. ' \
                                                                           'Please check the SliceNDice log file for details.'},}

    def processInputFiles(self):
        error = None
        self.hklin = None
        self.xyzin = None
        dataObjects = []
        #Append Observation with representation dependent on whether we are detwining on Is or not
        dataObjects += [['F_SIGF',CCP4XtalData.CObsDataFile.CONTENT_FLAG_FMEAN]]
        #Include FreeRflag if called for
        if self.container.inputData.FREERFLAG.isSet():
            dataObjects += ['FREERFLAG']
        self.hklin,error = self.makeHklin(dataObjects)
        if error.maxSeverity()>CCP4ErrorHandling.SEVERITY_WARNING:
            return CPluginScript.FAILED
        else:
            return CPluginScript.SUCCEEDED
    
    def makeCommandAndScript(self):

        inp = self.container.inputData
        par = self.container.controlParameters
        mod = self.container.modelParameters
        gui = self.container.guiParameters
        out = self.container.outputData
    
        # Set the max number of processors
        MAXPROC=multiprocessing.cpu_count()  
    
        if self.container.inputData.XYZIN.isSet():
            self.xyzin = self.container.inputData.XYZIN.__str__()
    
        self.appendCommandLine( [ '--xyzin', self.xyzin ] )
        self.appendCommandLine( [ '--hklin', self.hklin ] )
        seqFile = os.path.join(self.workDirectory,'SEQIN.fasta')
        inp.ASUIN.writeFasta(fileName=seqFile)
        self.appendCommandLine( [ '--seqin', seqFile ] )
        if self.container.modelParameters.BFACTOR_TREATMENT:
            self.appendCommandLine( [ '--bfactor_column', str(self.container.modelParameters.BFACTOR_TREATMENT) ] )
            if str(self.container.modelParameters.BFACTOR_TREATMENT) == "plddt":
                self.appendCommandLine( [ '--plddt_threshold', str(self.container.modelParameters.PLDDT_THRESHOLD) ] )
            elif str(self.container.modelParameters.BFACTOR_TREATMENT) == "rms":
                self.appendCommandLine( [ '--rms_threshold', str(self.container.modelParameters.RMS_THRESHOLD) ] )
        if self.container.modelParameters.MIN_SPLITS:
            self.appendCommandLine( [ '--min_splits', str(self.container.modelParameters.MIN_SPLITS) ] )
        if self.container.modelParameters.MAX_SPLITS:
            self.appendCommandLine( [ '--max_splits', str(self.container.modelParameters.MAX_SPLITS) ] )
        if self.container.controlParameters.NPROC:
            self.appendCommandLine( [ '--nproc', str(self.container.controlParameters.NPROC) ] )
        if self.container.controlParameters.NCYC:
            self.appendCommandLine( [ '--ncyc_refmac', str(self.container.controlParameters.NCYC) ] )
        if self.container.controlParameters.NO_MOLS:
            self.appendCommandLine( [ '--no_mols', str(self.container.controlParameters.NO_MOLS) ] )
        self.appendCommandScript( "" )
        # Create output xml file
        self.xmlout = self.makeFileName('PROGRAMXML')
        return CPluginScript.SUCCEEDED

#    """
#    def postProcess( self, processId=-1, data={} ) :
#      out = self.container.outputData
#      xmlout = str( self.makeFileName( 'PROGRAMXML' ) )
#      self.reportStatus(0)
#    """
#

    # process one or more output files
    # also writes the XML file, previously done by postProcess()

    def processOutputFiles(self):
        # Load Json
        try:
            jsfloc = os.path.join(self.getWorkDirectory(), "slicendice_0", "slicendice_results.json")
            jfi = open(jsfloc)
            jdd = json.load(jfi)
            jfi.close()
        except:
            # Failed to find a solution in the json file.
            self.appendErrorReport(19121)
            print("SlicenDice: NO json output found.")
            return CPluginScript.FAILED
        # Get the best soln from the json file
        print("HERE:")
        try:
            #rfrl = list(jdd.get('final_r_free').values())
            lowest_rfree=1.0
            best_split=None
            for split in jdd['dice'].keys():
                if jdd['dice'][split]['final_r_free'] <= lowest_rfree:
                    lowest_rfree=jdd['dice'][split]['final_r_free']
                    best_split=split
            #rfrl = jdd['dice'][best_split]['final_r_free']
            #lindx = rfrl.index(min(rfrl))
            #bskey = list(jdd.get('final_r_free').keys())[lindx]
        except:
            # Failed to find a solution in the json file.
            self.appendErrorReport(19122)
            print("SlicenDice: NO solution found in the json outfile.")
            return CPluginScript.FAILED
        # Need to copy the files into the actual project directory - cannot be a sub-directory. Not entirely sure why but...
        #xyz = os.path.normpath(jdd.get('xyzout').get(bskey))
        #hkl = os.path.normpath(jdd.get('hklout').get(bskey))
        xyz = os.path.normpath(jdd['dice'][best_split]['xyzout'])
        hkl = os.path.normpath(jdd['dice'][best_split]['hklout'])
        xyzout = os.path.join(self.getWorkDirectory(),os.path.basename(xyz))
        hklout = os.path.join(self.getWorkDirectory(),os.path.basename(hkl))

        if os.path.isfile(xyz): 
            shutil.copy2(xyz, xyzout)
            self.container.outputData.XYZOUT=xyzout
        if os.path.isfile(hkl): 
            shutil.copy2(hkl, hklout)
            self.container.outputData.HKLOUT=hklout
        
        # Need to set the expected content flag  for phases data
        self.container.outputData.XYZOUT.annotation     = 'Model from SliceNDice refinement'
        self.container.outputData.FPHIOUT.annotation    = 'Weighted map from SliceNDice refinement'
        self.container.outputData.DIFFPHIOUT.annotation = 'Weighted difference map from SliceNDice refinement'

        # Split out data objects that have been generated. Do this after applying the annotation, and flagging
        # above, since splitHklout needs to know the ABCDOUT contentFlag
        
        outputFiles = ['FPHIOUT', 'DIFFPHIOUT']
        outputColumns = ['FWT,PHWT', 'DELFWT,PHDELWT']
        #if self.container.controlParameters.PHOUT:
        #    outputFiles+=['ABCDOUT']
        #    outputColumns+=['HLACOMB,HLBCOMB,HLCCOMB,HLDCOMB']
        error = self.splitHklout(outputFiles,outputColumns,infile=hklout)
        if error.maxSeverity() > CCP4ErrorHandling.SEVERITY_WARNING:
            return CPluginScript.FAILED

        #Set performance indicators
        try:
            #bid = jdd.get('split_id').get(bskey)
            #rrfr = [jdd.get('final_r_fact').get(bskey), jdd.get('final_r_free').get(bskey)]
            bid = best_split.split("_")[-1]
            rrfr = [jdd['dice'][best_split]['final_r_fact'], jdd['dice'][best_split]['final_r_free']]
            self.container.outputData.PERFORMANCEINDICATOR.RFactor.set(str(rrfr[0]))
            self.container.outputData.PERFORMANCEINDICATOR.RFree.set(str(rrfr[1]))
        except:
            print("Failed to load r-factors from json log")
        # xml info
        rootNode = ET.Element("SliceNDice")
        xmlRI = ET.SubElement(rootNode, "RunInfo")
        xmlbcyc = ET.SubElement(xmlRI, "Best")
        ET.SubElement(xmlbcyc, "bid").text = str(bid)
        ET.SubElement(xmlbcyc, "R").text = str(rrfr[0])
        ET.SubElement(xmlbcyc, "RFree").text = str(rrfr[1])
        # Get solns & save
        #for key in jdd.get('split_id').keys():
        for key in jdd['dice'].keys():
            xmlcyc = ET.SubElement(xmlRI, "Sol")
            ET.SubElement(xmlcyc, "SolID").text = str(key.split("_")[-1])
            ET.SubElement(xmlcyc, "llg").text = str(jdd['dice'][key]['phaser_llg'])
            ET.SubElement(xmlcyc, "tfz").text = str(jdd['dice'][key]['phaser_tfz'])
            ET.SubElement(xmlcyc, "srf").text = str(jdd['dice'][key]['final_r_fact'])
            ET.SubElement(xmlcyc, "sre").text = str(jdd['dice'][key]['final_r_free'])
            #ET.SubElement(xmlcyc, "SolID").text = str(jdd.get('split_id').get(key))
            #ET.SubElement(xmlcyc, "llg").text = str(jdd.get('phaser_llg').get(key))
            #ET.SubElement(xmlcyc, "tfz").text = str(jdd.get('phaser_tfz').get(key))
            #ET.SubElement(xmlcyc, "srf").text = str(jdd.get('final_r_fact').get(key))
            #ET.SubElement(xmlcyc, "sre").text = str(jdd.get('final_r_free').get(key))
        # Save xml
        xmlfile = open(self.xmlout, 'wb')
        ET.indent(rootNode)
        xmlString= ET.tostring(rootNode)
        xmlfile.write(xmlString)
        xmlfile.close()
        return CPluginScript.SUCCEEDED

#------------------------------------------------------------------------------------
import unittest

class testslicendice( unittest.TestCase ) :

#- def setUp( self ) :
#- def tearDown( self ) :
#- def test2( self ) :

   def test1( self ) :
      from core.CCP4Utils import getCCP4I2Dir
      xmlInput = os.path.join( getCCP4I2Dir(), 'wrappers', 'slicendice', 'test_data', 'test1'+'.params.xml' )
      self.wrapper = slicendice( name='job' )
      self.wrapper.container.loadDataFromXml( xmlInput )
      self.wrapper.setWaitForFinished( 1000000 )

      pid = self.wrapper.process()
      self.wrapper.setWaitForFinished( -1 )
      if len(self.wrapper.errorReport)>0:
         print(self.wrapper.errorReport.report())

def TESTSUITE() :

   suite = unittest.TestLoader().loadTestsFromTestCase( testslicendice )
   return suite

def testModule() :

   suite = TESTSUITE()
   unittest.TextTestRunner( verbosity=2 ).run( suite )

