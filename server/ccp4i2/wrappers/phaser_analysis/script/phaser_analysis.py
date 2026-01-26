"""
Run Phaser to analyse merged data, mainly to get the plot of
 average information content per reflection as a function of resolution,
 in order to estimate the useful resolution

See:
  Measuring and using information gained by observing diffraction data
  R.J. Read, R.D. Oeffner and A.J. McCoy  Acta Cryst. (2020). D76, 238-247
  https://doi.org/10.1107/S2059798320001588
"""

import math
import os

import phaser
from lxml import etree

from ccp4i2.core.CCP4PluginScript import CPluginScript


class CallbackObject(object):
    # This will become the output object for the Phaser run
    #  so that the loggraphs are separated from the rest of the logfile
    # Other entries are dummies (I think)
    def __init__(self):
        self.loggraphs = []

    def loggraph (self, title, data):
        # Accumulate loggraphs
        #print(":loggraph ", title)
        self.loggraphs.append([title,data])
    
    def getloggraphs(self):
        return self.loggraphs

    # Dummies
    def startProgressBar (self, label, size):
        pass
    def incrementProgressBar (self):
        pass
    def endProgressBar (self):
        pass
    def warn (self, message):
        pass
    def call_back (self, message, data):
        pass

###### end class CallbackObject


class phaser_analysis(CPluginScript):
    TASKMODULE = 'test'      # Where this plugin will appear on the gui
    TASKTITLE = 'phaser_analysis' # A short title for gui menu
    TASKNAME = 'phaser_analysis'   # Task name - should be same as class name
    TASKVERSION= 0.0               # Version of this plugin
    MAINTAINER = 'pre@mrc-lmb.cam.ac.uk'
 
    def startProcess(self):
        # Where are we?
        self.xmlout = self.makeFileName( 'PROGRAMXML' )

        # Data from input
        filename = self.container.inputData.HKLIN.fullPath
        print("Filename: ",filename)

        self.pxdname = str(self.container.inputData.PXDNAME)

        threshold = self.container.controlParameters.INFOCONTENTTHRESHOLD
        self.threshold = float(threshold)

        inputObject = phaser.InputMR_DAT()
        inputObject.setHKLI(str(filename))
        inputObject.setMUTE(True)
        resultObject = phaser.runMR_DAT(inputObject)
        if resultObject == CPluginScript.FAILED:
            self.reportStatus(CPluginScript.FAILED)
            return CPluginScript.FAILED

        inputObject = phaser.InputNCS()
        inputObject.setSPAC_HALL(resultObject.getSpaceGroupHall())
        inputObject.setCELL6(resultObject.getUnitCell())
        inputObject.setREFL_DATA(resultObject.getREFL_DATA())
        inputObject.setINFO(True)
        inputObject.setMUTE(True)
        inputObject.setHKLO(False) # No output mtz

        # The call-back here is just (?) to separate out the loggraphs
        outputObject = phaser.Output()
        callbackobject = CallbackObject()
        outputObject.setPhenixCallback(callbackobject)

        try:
            self.resultObject = phaser.runNCS(inputObject, outputObject)
        except RuntimeError as e:
            print("Phaser failed")
            self.appendErrorReport(105, str(e))
            self.reportStatus(CPluginScript.FAILED)
            return CPluginScript.FAILED

        if not self.resultObject.Success():
            self.appendErrorReport(105, self.resultObject.ErrorName() + '-' + self.resultObject.ErrorMessage())
            self.reportStatus(CPluginScript.FAILED)
            return CPluginScript.FAILED

        # Get logfile and loggraphs
        self.logfile = self.resultObject.logfile()
        logfilename = self.makeFileName('LOG')
        with open(str(logfilename), 'w') as f:
            f.write(self.logfile)
        
        # Retrieve all loggraphs, a list of [title, data] pairs
        self.loggraphs = callbackobject.getloggraphs()

        self.makeXML(self.xmlout)  # makes self.xmlroot and writes out
        self.reportStatus(CPluginScript.SUCCEEDED)

        return CPluginScript.SUCCEEDED

    # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
    def makeXML(self, xmlout):
        from .phaser_analysis_utils import AnalysisLog

        # Uses self.logfile and self.loggraphs
        self.xmlroot = etree.Element('PHASER_ANALYSIS')
        self.xmlroot.set('name', self.pxdname)  # dataset name

        # Save all loggraphs to file (for now anyway)
        loggraphfilename = os.path.join(self.workDirectory, 'loggraphs.log')
        with open(loggraphfilename, 'w') as f:
            for loggraph in self.loggraphs:
                f.write(loggraph[1])

        # phaser.runNCS produces a number of loggraphs in the callback
        # These are the one we want keep, identified by title and a short label

        # Wanted graphs: title, name
        found = self.makeXMLloggraph('Weighted second moments after tNCS (best) refinement',
                                     '2ndmoment')
        # if "after tNCS" graph found, don't get the alternative
        if not found:
            self.makeXMLloggraph('Weighted second moments after Anisotropy Correction',
                                 '2ndmoment')
        # always get this one
        self.makeXMLloggraph('Information content of data', 'informationcontent')

        # XML from log file
        print("Analyse logfile")
        analyselogfile = AnalysisLog(self.logfile)
        self.xmlroot.append(analyselogfile.getXML())
        # Get maximum resolution of data from XML (from logfile)
        self.maxresolution = analyselogfile.getmaxresolution()

        self.resolutionlimit(self.threshold)

        with open(xmlout, 'wb') as f:
            f.write(etree.tostring(self.xmlroot, pretty_print=True))

    # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
    def makeXMLloggraph(self, graphtitle, gname):
        from .phaser_analysis_utils import Makexmlgraph

        # Make XML version of loggraph with title graphtitle and name gname
        # Returns True if this graph was found

        for loggraph in self.loggraphs:
            title = loggraph[0]
            # Do we want to keep this one? search for title
            if graphtitle  in title:
                #   converted to Makexmlgraph objects,
                #      named in a dictionary
                lgxml = Makexmlgraph(gname,loggraph[1])
                self.xmlroot.append(lgxml.getXML())
                return True

        #print("***graph not found", graphtitle)
        return False

    # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
    def resolutionlimit(self, threshold):
        from .phaser_analysis_utils import AnalyseGraph, addElement

        # estimate resolution limit from XML graph of information content
        # Add result to self.xmlroot
        
        analysegraphs =  AnalyseGraph(self.xmlroot)

        # Information content table
        name = 'informationcontent'
        xcol = '1/d^2'
        ycol = 'Selected'
        
        #  allOk = -1 no valid points, +1 all valid, else = 0 
        allOK, dslimit = analysegraphs.getresolutionlimit(name, xcol, ycol, threshold)
        reslimit = 1.0/math.sqrt(float(dslimit))

        # Store for access from pipeline
        self.allOK = allOK
        self.reslimit = reslimit

        message = None
        if allOK > 0:
            message = '== maximum resolution'
            # If data go to edge, then limit estimate is last bin value
            #  replace with maximum value in file
            reslimit = float(self.maxresolution)
        elif allOK < 0:
            message = 'No acceptable data'

        resolutionxml = etree.Element('ResolutionEstimate', type=name)
        addElement(resolutionxml, 'Threshold', "{:6.2f}".format(threshold))
        addElement(resolutionxml, 'Columns', xcol+' '+ycol)
        addElement(resolutionxml, 'ResolutionLimitEstimate', "{:6.2f}".format(reslimit))
        if message is not None:
            addElement(resolutionxml, 'Message', message)

        self.xmlroot.append(resolutionxml)

    # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
    def getresolutionlimit(self):
        return self.allOK, self.reslimit

    # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
    def getXML(self):
        return self.xmlroot
