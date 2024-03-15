from __future__ import print_function

"""
    import_merged.py: CCP4 GUI Project
     Copyright (C) 2015 STFC

     This library is free software: you can redistribute it and/or
     modify it under the terms of the GNU Lesser General Public License
     version 3, modified in accordance with the provisions of the 
     license to address the requirements of UK law.
 
     You should have received a copy of the modified GNU Lesser General 
     Public License along with this library.  If not, copies may be 
     downloaded from http://www.ccp4.ac.uk/ccp4license.php
 
     This program is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU Lesser General Public License for more details.
"""

import sys
import os,shutil
from PySide2 import QtCore
from core.CCP4PluginScript import CPluginScript
from core import CCP4Utils
from core.CCP4ErrorHandling import *
from lxml import etree
from pipelines.aimless_pipe.script.aimless_pipe_utils import *

from  pipelines.import_merged.script.mmcifutils import *
from  pipelines.import_merged.script.mmcifconvert import *
from  pipelines.import_merged.script.importutils import *
from  pipelines.import_merged.script.mtzimport import *

class import_merged(CPluginScript):

    TASKNAME = 'import_merged'
    MAINTAINER = 'liz.potterton@york.ac.uk'
    WHATNEXT = []
    ERROR_CODES = { 301 : { 'description' : 'No output file found after conversion program' },
                    302 : { 'description' : 'Output from conversion does not contain recognised reflection or FreeR set data' }
                    }
    # Note - preserving the HKLOUT by changing severity from the system default of 1 to 5 and
    # beware issues with caseinsensitivity
    PURGESEARCHLIST = [ [ 'HKLIN*.mtz' , 1 ],
                        ['aimless_pipe%*/HKLOUT*.mtz', 1],
                        [ 'hklout.mtz' , 5 ],    
                        [ 'HKLOUT.mtz' , 5 ]
                      ]
    #------------------------------------------------------------------------
    def process(self):
      self.container.inputData.HKLIN.loadFile()
      self.fformat = self.container.inputData.HKLIN.getFormat()
      # Type(self.fformat) can be either [mtz] <class 'CCP4Data.CString'> or [mmcif] str  WHY? 
      print("process self.fformat", type(self.fformat))
      merged = self.container.inputData.HKLIN.getMerged()
      self.isintensity = 0  # unknown I or F
      
      obsout = None

      self.x2mtz = None
      self.mmcifXML = None
      self.resolutioncutoff = False
      #print("IDRR", self.container.inputData.RESOLUTION_RANGE_SET)
      if self.container.inputData.RESOLUTION_RANGE_SET:
          self.resolutioncutoff = True

      if str(self.fformat) == 'mtz':
          # MTZ file: if there is a resolution cutoff specified,
          # cannot use x2mtz to run cmtzsplit
          if not self.resolutioncutoff:
              self.x2mtz = self.makePluginObject('x2mtz')
      #elif str(self.fformat) == 'mmcif':
      #    pass
      # # self.x2mtz = self.makePluginObject('cif2mtz')  # not needed now
      elif str(self.fformat) in [ 'shelx' ]:
        self.x2mtz = self.makePluginObject('convert2mtz')
      elif str(self.fformat) in [ 'sca' ]:
        self.x2mtz = self.makePluginObject('scalepack2mtz')

      if self.x2mtz is not None:
          #  Copy parameters to x2mtz sub-object
          self.x2mtz.container.inputData.copyData(otherContainer=self.container.inputData)
          self.x2mtz.container.outputData.copyData(otherContainer=self.container.outputData,dataList=['HKLOUT','OBSOUT'])

      self.freeRcompleteTried = True
      self.importXML = None
      self.freeout = None
      if self.fformat == 'mtz':
        if self.resolutioncutoff:
            self.importXML = ET.Element('IMPORT_LOG')  # information about the import step
            status = self.importmtz()
            self.makeReportXML(self.importXML)  # add initial stuff for XML into self.importXML
            self.process1(status)
            return  # Probably doesnt get here
        # No resolution cutoff
        # Just call the processOutputFiles() to convert to mini mtzs
        self.x2mtz.container.outputData.HKLOUT.set(self.container.inputData.HKLIN)
        # Pick up dataset name (and crystal name if available)
        fcontent = self.container.inputData.HKLIN.getFileContent()
        #print 'input file content', type(fcontent), fcontent
        #print 'columns', fcontent.listOfColumns
        # check if selection columns are intensity or amplitudes
        # return +1 if intensity, -1 if amplitude, 0 if unknown
        self.isintensity = self.isIntensity(self.container.inputData.HKLIN_OBS_COLUMNS,
                                            fcontent.listOfColumns)
        self.importXML = ET.Element('IMPORT_LOG')  # information about the import step
        self.makeReportXML(self.importXML)  # add initial stuff for XML into self.importXML

        if len(fcontent.datasets)>=2:
          #print fcontent.datasets[1]
          # self.crystalName = ###  set this, but how?
          self.container.inputData.DATASETNAME = fcontent.datasets[1]


        if self.container.controlParameters.SKIP_FREER:
            self.x2mtz.container.outputData.FREEOUT.set(self.container.outputData.FREEOUT)
        self.x2mtz.checkOutputData()
        ret = self.x2mtz.processOutputFiles()  # runs cmtzsplit
        if sys.platform != 'win32': # CCP4Utils.samefile() doesn't work (r1728)
          self.x2mtz.reportStatus(ret)
        self.container.outputData.OBSOUT.set(self.x2mtz.container.outputData.OBSOUT)
        if self.importXML is not None:
            freeRcolumnLabel = str(self.x2mtz.freeRcolumnLabel)
            if freeRcolumnLabel is not None:
                self.addElement(self.importXML, 'freeRcolumnLabel',
                                freeRcolumnLabel)
            self.outputLogXML(self.importXML)  # send self.importXML to program.xml
        self.process1({'finishStatus': ret })
      else:
          # not MTZ
          self.importXML = ET.Element('IMPORT_LOG')  # information about the import step
          #  +1 if intensity, -1 if amplitude, 0 if unknown
          self.isintensity = 0
          if str(self.fformat) in [ 'sca' ]:
              self.isintensity = +1  # scalepack files are intensity
          if self.container.inputData.MMCIF_SELECTED_ISINTENSITY:
              self.isintensity = self.container.inputData.MMCIF_SELECTED_ISINTENSITY

          self.makeReportXML(self.importXML)  # add initial stuff for XML into self.importXML
          self.outputLogXML(self.importXML)  # send self.importXML to program.xml

          # mmCIF, direct import
          if str(self.fformat) == 'mmcif':
              status = self.convertmmcif()
              self.process1(status)
              return  # Probably doesnt get here

          self.connectSignal(self.x2mtz,'finished',self.process1slot)
          ret = self.x2mtz.process()
        
    #------------------------------------------------------------------------
    @QtCore.Slot(dict)
    def process1slot(self,status):
        self.process1(status)

    def process1(self,status, completeFreeR=True):
        'if completeFreeR False, always generate new FreeR (for 2nd attempt)'
        #print('process1',type(status),status)
        if status is not None and status.get('finishStatus') == CPluginScript.FAILED:
            self.reportStatus(status)
            return
      
        # Is FreeR generation switched off?
        if self.container.controlParameters.SKIP_FREER:
            # No freeR generation, leave as is, eg from STARANISO
            self.process2(CPluginScript.SUCCEEDED)

        if not self.container.inputData.HASFREER:
            completeFreeR = False   # no valid FreeR data

        # Create or complete a freer set
        self.freerflag = self.makePluginObject('freerflag')
        if self.x2mtz is not None:
            self.freerflag.container.inputData.F_SIGF = \
                            self.x2mtz.container.outputData.OBSOUT
        else:
            self.freerflag.container.inputData.F_SIGF = \
                        self.container.outputData.OBSOUT
            
        #print 'import_merged.process1',self.x2mtz.container.outputData.FREEOUT,self.x2mtz.container.outputData.FREEOUT.exists()
        newfreer = 'True'
        freeRsource = None
        if completeFreeR:
            # Cases if completeFreeR == True:
            #  1) inputData.FREERFLAG is set (FreeR from separate object|file), complete this one, or
            #     freeRsource = 'Explicit'
            #  2) FREEOUT.exists from main import, complete this
            #     freeRsource = 'Input'
            #  3) else generate new one
            if self.container.inputData.FREERFLAG.isSet():  # case (1)
                self.freerflag.container.inputData.FREERFLAG = self.container.inputData.FREERFLAG
                freeRsource = 'Explicit'
                #  Check compatible cells etc
                #print "****"
                #print "FR", self.freerflag.container.inputData.FREERFLAG.fileContent
                #print "OBSOUT", self.container.outputData.OBSOUT.fileContent
                #print "HKLOUT", self.container.outputData.HKLOUT.fileContent

                tolerance = None # use default
                cellcheck = \
                          CellCheck(self.freerflag.container.inputData.F_SIGF.fileContent,
                                    self.freerflag.container.inputData.FREERFLAG.fileContent,
                                    tolerance)
                cellsAreTheSame, freerReportXML = cellcheck.checks()
                if not cellsAreTheSame['validity']:
                    completeFreeR = False

                if freerReportXML is not None:
                    self.importXML.append(freerReportXML)
            elif self.freeout is not None:
                # FREEOUT from mmcif
                # A freeR set has been imported, so extend/complete it
                self.freerflag.container.inputData.FREERFLAG.set(self.freeout)
                freeRsource = 'Input'

            elif self.x2mtz.container.outputData.FREEOUT.exists():
                # A freeR set has been imported, so extend/complete it
                self.freerflag.container.inputData.FREERFLAG = self.x2mtz.container.outputData.FREEOUT                
                freeRsource = 'Input'
            else:
                completeFreeR = False

        if completeFreeR:
            self.freerflag.container.controlParameters.GEN_MODE = 'COMPLETE'
            self.freerflag.container.controlParameters.COMPLETE = True
            self.freerflag.container.controlParameters.CUTRESOLUTION = \
                          self.container.controlParameters.CUTRESOLUTION
            newfreer = 'False'
          
        self.freerflag.container.controlParameters.FRAC = \
                         self.container.controlParameters.FREER_FRACTION
        self.addElement(self.importXML, 'newFreeR', newfreer) 
        if freeRsource is not None:
            self.addElement(self.importXML, 'freeRsource', freeRsource) 
        self.outputLogXML(self.importXML)  # send self.importXML to program.xml
        self.freerflag.container.outputData.FREEROUT.setFullPath(str(self.container.outputData.FREEOUT))

        self.connectSignal(self.freerflag,'finished',self.process2)
        status = self.freerflag.process()
        
    #------------------------------------------------------------------------
    @QtCore.Slot(dict)
    def process2(self,status):
      freerOK = True
      doFreeR = True
      if self.container.controlParameters.SKIP_FREER:
        doFreeR = False
      else:
        if status is not None and status.get('finishStatus') == CPluginScript.FAILED:
          # FreeR run has failed, create error message and continue
          self.addElement(self.importXML, 'FreeRfailed', 'True')
          self.outputLogXML(self.importXML)  # send self.importXML to program.xml
          # try again
          if self.freeRcompleteTried:
              print("trying again")
              self.freeRcompleteTried = False
              self.process1(None, False)  # try to make a new FreeR set
          else:
              # failed a 2nd time, make error and continue
              print("failed again")
              self.addElement(self.importXML, 'FreeRfailed', 'Again')
              self.outputLogXML(self.importXML)  # send self.importXML to program.xml
              freerOK = False
              self.reportStatus(status)
              return

        # If FreeR is OK, then create annotation
        if freerOK:
          #print "FREEROUT content",self.freerflag.container.outputData.FREEROUT.fileContent
          if self.freerflag.container.outputData.FREEROUT.fileContent.spaceGroup.isSet():
              sgname = self.freerflag.container.outputData.FREEROUT.fileContent.spaceGroup.__str__()
          else:
              sgname = 'Unk'

          highresFRformatted = "%7.2f" % float(self.freerflag.container.outputData.FREEROUT.fileContent.resolutionRange.high)
          title ='FreeR - Spg:'+str(sgname).strip()+';Resln:'+highresFRformatted.strip() + "A;"
          try:
              title = title + "Cell:"+self.freerflag.container.outputData.FREEROUT.fileContent.cell.guiLabel()
          except Exception as e:
              print('Error writing cell parameters',e)

          #print "FreeR title:",title
          self.container.outputData.FREEOUT.annotation = title

      # Add x2mtz XML if present
      if self.x2mtz is not None:
          x2mtz_XMLpath = self.x2mtz.makeFileName('PROGRAMXML')
          x2mtz_Etree = ET.parse(x2mtz_XMLpath)
          x2mtz = x2mtz_Etree.getroot()
          self.importXML.append(x2mtz)
      if self.mmcifXML is not None:
          self.importXML.append(self.mmcifXML)

      # add in FreeR XML
      if doFreeR:
          self.importXML.append(self.freerflag.getXML())
      else:
          self.addElement(self.importXML, 'freeRsource', 'None') 

      # Run aimless for a report on data quality
      self.aimlesspipe = self.makePluginObject('aimless_pipe',pluginTitle='DR run for data analysis')
      unmergedList = self.aimlesspipe.container.inputData.UNMERGEDFILES
      #print '\nunmergedList 0',    unmergedList
      if len(unmergedList)==0: unmergedList.addItem()
      # Always do analysis on the file which is saved as pipeline output
      unmergedList[0].file.set(self.container.outputData.OBSOUT.__str__())
      ##  earlier versions in some cases analysed the input file
      #if self.fformat in ['mmcif']:
      #    unmergedList[0].file.set(self.container.outputData.HKLOUT.__str__())
      #      elif  self.fformat in ['mtz']:
      #          unmergedList[0].file.set(self.container.outputData.OBSOUT.__str__())
      #      else:
      #          unmergedList[0].file.set(self.container.inputData.HKLIN.__str__())
      #print 'unmergedList 1',    unmergedList
      xname = self.filteredName(str(self.container.inputData.CRYSTALNAME), 'X')
      #xname = CCP4Utils.safeOneWord(str(self.container.inputData.CRYSTALNAME))
      dname = self.filteredName(str(self.container.inputData.DATASETNAME), 'D')
      #dname = CCP4Utils.safeOneWord(str(self.container.inputData.DATASETNAME))
      unmergedList[0].crystalName.set(xname)
      unmergedList[0].dataset.set(dname)
      #print 'unmergedList 2',    unmergedList
      unmergedList[0].cell.set(self.container.inputData.SPACEGROUPCELL.cell)
      #print 'unmergedList 3',    unmergedList
      #print 'self.container.inputData',self.container.inputData
      #print 'self.container.inputData.SPACEGROUPCELL',\
      #      self.container.inputData.SPACEGROUPCELL
      unmergedList[0].wavelength.set(self.container.inputData.WAVELENGTH)

      # parameters for Pointless
      self.aimlesspipe.container.controlParameters.MODE = 'CHOOSE'
      self.aimlesspipe.container.controlParameters.CHOOSE_MODE = 'SPACEGROUP'
      self.aimlesspipe.container.controlParameters.CHOOSE_SPACEGROUP = \
            self.container.inputData.SPACEGROUPCELL.spaceGroup
      # parameters for Aimless
      self.aimlesspipe.container.controlParameters.SCALING_PROTOCOL = 'CONSTANT'
      self.aimlesspipe.container.controlParameters.ONLYMERGE = True
      self.aimlesspipe.container.controlParameters.ANALYSIS_MODE = True
      self.aimlesspipe.container.controlParameters.OUTPUT_UNMERGED = False
      self.aimlesspipe.container.controlParameters.SDCORRECTION_OVERRIDE = True
      self.aimlesspipe.container.controlParameters.SDCORRECTION_REFINE = False
      self.aimlesspipe.container.controlParameters.SDCORRECTION_SET = True
      self.aimlesspipe.container.controlParameters.SDCORRECTION_SDFAC = 1.0
      self.aimlesspipe.container.controlParameters.SDCORRECTION_SDB = 0.0
      self.aimlesspipe.container.controlParameters.SDCORRECTION_SDADD = 0.0
      

#  Probably shouldn't run Phaser but try anyway
#      if self.isintensity < 0:
#         print("* Fs input, don't run Phaser")
#          self.aimlesspipe.container.controlParameters.DOPHASERANALYSIS = False

      tempXML = self.importXML
      self.addElement(tempXML, "DRPIPE_RUNNING", "True") 
      self.outputLogXML(tempXML)  # SEND self.importXML to program.xml
      #  Start data reduction
      self.connectSignal(self.aimlesspipe,'finished',self.nearlyDone)
      print("starting aimless_pipe")
      self.aimlesspipe.process()

      return CPluginScript.SUCCEEDED

    #------------------------------------------------------------------------
    @QtCore.Slot(dict)
    def nearlyDone(self,status):
      print('import_merged.nearlyDone')
      self.container.outputData.OBSOUT.setContentFlag(reset=True)
      import shutil
      try:
          # XML data: We have
          #   a) self.importXML etree element report on the import step
          #   b) aimless pipe report
          # so put these together into an IMPORT_MERGED block
          aimless_pipe_XMLpath = self.aimlesspipe.makeFileName('PROGRAMXML')
          aimless_pipe_Etree = ET.parse(aimless_pipe_XMLpath)
          aimless_pipe = aimless_pipe_Etree.getroot()
          #print 'aimless_pipe', type(aimless_pipe), aimless_pipe
          self.outputLogXML(self.importXML, aimless_pipe)  # and output it
      except:
        pass

      print('import_merged.nearlyDone, finished')
      self.reportStatus(status)

    #------------------------------------------------------------------------
    def isIntensity(self, selectedcolumns, listOfColumns):
        """
        check if selection columns are intensity or amplitudes
        return +1 if intensity, -1 if amplitude, 0 if unknown
        selectedcolumns  string of columns labels
        listOfColumns  list eg
        [{'groupIndex': '1', 'columnLabel': 'F_xe1a', 'columnType': 'F', 'dataset': 'xe1a'},
        """
        #print 'isIntensity',selectedcolumns,type(selectedcolumns)
        if not selectedcolumns.isSet():
          selcol1 = str(listOfColumns[0].get()['columnLabel'])
        else:
          selected = selectedcolumns.split(',')
          selcol1 = selected[0] # first label
        #print 'selected', selectedcolumns, selcol1
        columntype = None
        for col in listOfColumns:
            column = col.get()
            if column['columnLabel'] == selcol1:
                #print 'found', column['columnLabel'], column['columnType']
                columntype = column['columnType']
        isintensity = 0
        if columntype is None:
            print("Unrecognised column " + selcol1)
        elif (columntype == 'J') or (columntype == 'K'):
            isintensity = +1
        elif (columntype == 'F') or (columntype == 'G'):
            isintensity = -1

        return isintensity

    #------------------------------------------------------------------------
    def outputLogXML(self, x1XML, x2XML=None):
      'output x1XML and optionally x2XML to program.xml'
      #print "outputLogXML", x1XML
      rootXML = ET.Element('IMPORT_MERGED') # Global XML for everything
      rootXML.append(x1XML)
      if x2XML is not None:
          rootXML.append(x2XML)
      with open (self.makeFileName('PROGRAMXML'),"w") as outputXML:
          ET.indent(rootXML)
          CCP4Utils.writeXML(outputXML,ET.tostring(rootXML))

    #------------------------------------------------------------------------
    def makeReportXML(self, containerXML):
        'Make initial report XML'
        #print("mrx 1",  type(self.fformat),  self.fformat)
        filename = str(self.container.inputData.HKLIN)  # Input file
        self.addElement(containerXML, 'filename', filename)
        ffmt = str(self.fformat)
        self.addElement(containerXML, 'fileformat', ffmt)

        #print 'type merged', type(self.container.inputData.HKLIN.getMerged())
        if self.container.inputData.HKLIN.getMerged():
            self.addElement(containerXML, 'merged', 'True')
        else:
            self.addElement(containerXML, 'merged', 'False')
        if self.fformat == 'mtz':
            self.addElement(containerXML, 'columnlabels',
                            self.container.inputData.HKLIN_OBS_COLUMNS.get())
        # = +1 if intensity, -1 if amplitude, 0 if unknown
        if self.isintensity == 0:
            IorFtype = 'Unknown'
        elif self.isintensity > 0:
            IorFtype = 'Intensity'
        elif self.isintensity < 0:
            IorFtype = 'Amplitude'
        self.addElement(containerXML, 'IorFtype', IorFtype)

        if self.container.controlParameters.STARANISO_DATA:
            self.addElement(containerXML, 'StarAniso', 'True')
        
        if self.fformat == 'sca':
            # Scalepack
            resorange = self.makeResoRange()
            if resorange is not None:
                # XML version of resolution range, a tuple of (dmax, dmin)
                resoxml = ET.Element('ResolutionRange')
                resoxml.set('id', 'cutresolution')
                if resorange[0] > 0.0:
                    self.addElement(resoxml, 'min', "{:.3f}".format(resorange[0]))
                if resorange[1] > 0.0:
                    self.addElement(resoxml, 'max', "{:.3f}".format(resorange[1]))
                containerXML.append(resoxml)
        # mmCIF things
        if ffmt == 'mmcif':
            if self.container.inputData.MMCIF_SELECTED_BLOCK.isSet():
                self.addElement(containerXML, 'mmcifblock',
                                str(self.container.inputData.MMCIF_SELECTED_BLOCK))
            if self.container.inputData.MMCIF_SELECTED_DETAILS.isSet():
                self.addElement(containerXML, 'mmcifblockdetails',
                                str(self.container.inputData.MMCIF_SELECTED_DETAILS))
            if self.container.inputData.MMCIF_SELECTED_INFO.isSet():
                self.addElement(containerXML, 'mmcifblockinfo',
                                str(self.container.inputData.MMCIF_SELECTED_INFO))
            if self.container.inputData.MMCIF_SELECTED_COLUMNS.isSet():
                self.addElement(containerXML, 'mmcifblockcolumns',
                                str(self.container.inputData.MMCIF_SELECTED_COLUMNS))
            
    #------------------------------------------------------------------------
    def addElement(self, containerXML, elementname, elementtext):
        #print 'addElement', elementname, type(elementtext), elementtext 
        e2 = ET.Element(elementname)
        e2.text = elementtext
        containerXML.append(e2)

    #------------------------------------------------------------------------
    def filteredName(self,name, default=None):
        """ filtered to remove spaces and other rubbish """
        if name is None: return ""
        if len(name) == 0:
            if default is not None: return default
            return ""
        #  funny things with backslash
        if (name.find('\\')):
            name = name.split('\\')[0]
            # only letters and numbers, no spaces or other stuff
        return CCP4Utils.safeOneWord(name)

    #------------------------------------------------------------------------
    def convertmmcif(self):
        # Convert an mmcif file to a data file (OBSOUT) and
        #  (if present) a FreeR file
        # Uses Gemmi

        filename = str(self.container.inputData.HKLIN)  # Input file
        blockname = None    # selected block
        if self.container.inputData.MMCIF_SELECTED_BLOCK.isSet():
            blockname = str(self.container.inputData.MMCIF_SELECTED_BLOCK)

        outfile = str(self.container.outputData.OBSOUT)

        if self.container.controlParameters.SKIP_FREER:
            freerfile = str(self.container.outputData.FREEOUT)
        else:
            # Temporary place for FreeR in job_1 subdirectory
            wd = self.getWorkDirectory()
            wdir = os.path.join(wd, 'job_1')
            if not os.path.exists(wdir):
                os.mkdir(wdir, 0o777)
                freerfile = os.path.join(wdir, 'FREEOUT.mtz')

        self.freeout = freerfile
        reducehkl = True  # for now

        #print("convertmmcif files", outfile, freerfile)
        resorange = self.makeResoRange()
            
        convertcif = ConvertCIF(filename, blockname,
                                outfile, freerfile, reducehkl, resorange)

        self.mmcifXML = convertcif.getXML()
        status = {'finishStatus':CPluginScript.FAILED}
        if convertcif.getstatus():
            contentFlag = convertcif.contentFlag()
            self.container.outputData.OBSOUT.contentFlag.set(contentFlag)
            status = {'finishStatus':CPluginScript.SUCCEEDED}

        return status

    # ---------------------------------------------------------------------------
    def makeResoRange(self):
        resorange = None
        if self.container.inputData.RESOLUTION_RANGE:
            dmax = 0.0
            dmin = 0.0
            r1 = self.container.inputData.RESOLUTION_RANGE.start
            if r1.isSet():
                dmax = float(r1)
            r2 = self.container.inputData.RESOLUTION_RANGE.end
            if r2.isSet():
                dmin = float(r2)
            if r1.isSet() or r2.isSet():
                resorange = (dmax, dmin)
        return resorange

    #------------------------------------------------------------------------
    def importmtz(self):
        # Import an mtz file to a data file (OBSOUT) and
        #  (if present) a FreeR file, with optional resolution cutoffs
        # Uses Gemmi

        outputData = self.container.outputData
        filename = str(self.container.inputData.HKLIN)  # Input file
        outfile = str(outputData.OBSOUT)

        # Input file column labels for observed data and
        #   FreeR (blank '' if no FreeR)
        obsColLabels, freeRcolumnLabel = self.columnthings(filename)
        obsColLabels = list(obsColLabels.split(','))

        if freeRcolumnLabel == '':
            freerfile = None
            self.container.inputData.HASFREER.set(False)
        else:
            if self.container.controlParameters.SKIP_FREER:
                freerfile = str(self.container.outputData.FREEOUT)
            else:
                # Temporary place for FreeR in job_1 subdirectory
                wd = self.getWorkDirectory()
                wdir = os.path.join(wd, 'job_1')
                if not os.path.exists(wdir):
                    os.mkdir(wdir, 0o777)
                    freerfile = os.path.join(wdir, 'FREEOUT.mtz')

        self.freeout = freerfile
        reducehkl = True  # for now
        resorange = self.makeResoRange()

        mtzimport = ImportMTZ(filename, outfile, freerfile,
                              obsColLabels, int(self.contentFlag),
                              freeRcolumnLabel,
                              resorange)

        self.mtzXML = mtzimport.getXML()
        self.importXML.append(self.mtzXML)
        status = {'finishStatus':CPluginScript.FAILED}
        if mtzimport.getstatus():
            status = {'finishStatus':CPluginScript.SUCCEEDED}
        status = {'finishStatus':CPluginScript.SUCCEEDED}
        return status

    # -------------------------------------------------------------------------
    def columnthings(self, filename):
        #  Sort out which columns are wanted, cf x2mtz.py
        print('HKLIN_OBS_COLUMNS', self.container.inputData.HKLIN_OBS_COLUMNS)
        print('HKLIN_OBS_CONTENT_FLAG', self.container.inputData.HKLIN_OBS_CONTENT_FLAG)
        #  HKLIN_OBS_COLUMNS   list of wanted input columns, this should be set
        inputData = self.container.inputData
        outputData = self.container.outputData
        inputData.HKLIN_OBS.fileContent.loadFile(filename)
        columnGroups = \
              inputData.HKLIN_OBS.fileContent.getColumnGroups()
        iBestObs, ifree = self.bestcolumns(columnGroups)
        if inputData.HKLIN_OBS_COLUMNS.isSet():
            obsColLabels = str(self.container.inputData.HKLIN_OBS_COLUMNS)
            self.contentFlag = self.container.inputData.HKLIN_OBS_CONTENT_FLAG
        else:
            # Obs columns not already set (should not happen)
            self.contentFlag = columnGroups[iBestObs].contentFlag
            outputData.OBSOUT.contentFlag.set(self.contentFlag)
            outputData.OBSOUT.annotation.set(\
                outputData.OBSOUT.qualifiers('guiLabel')+' from '+\
                inputData.HKLIN.stripedName())
            obsColLabels = str(columnGroups[iBestObs].columnList[0].columnLabel)
            for col in columnGroups[iBestObs].columnList[1:]:
                obsColLabels += ','+str(col.columnLabel)

        #  +1 if intensity, -1 if amplitude, 0 if unknown
        self.isintensity = +1
        if self.contentFlag%2 == 0:
            # F amplitudes
            self.isintensity = -1
        freeRcolumnLabel = ''
        if ifree >= 0:
            #  Freer column, if present
            outputData.FREEOUT.annotation.set(\
                outputData.FREEOUT.qualifiers(\
                    'guiLabel')+' from '+inputData.HKLIN.stripedName())
            freeRcolumnLabel = str(columnGroups[ifree].columnList[0].columnLabel)
        return obsColLabels, freeRcolumnLabel

    # -----------------------------------------------------------------------
    def bestcolumns(self, columnGroups):
        # Try to figure the 'best' obs data and freer data in the input file
        #for cg in columnGroups:
        #    print('bestcolumns cg',cg.get())
        #   print("* type",cg.columnGroupType)

        ifree = -1
        iBestObs = -1
        for ii in range(len(columnGroups)):
            if columnGroups[ii].columnGroupType == 'FreeR':
                ifree = ii  #  Index for FreeR
            elif columnGroups[ii].columnGroupType == 'Obs':
              # print( "@cg ", columnGroups[ii].contentFlag, columnGroups[iBestObs].contentFlag)
              if iBestObs<0 or \
               columnGroups[ii].contentFlag<columnGroups[iBestObs].contentFlag:
                  # Now the "best" obs column based on content type
                  iBestObs = ii
        return iBestObs, ifree   # indices to columngroups for Obs and Free

"""
# Function to return list of names of exportable MTZ(s)
def exportJobFile(jobId=None,mode=None):
    import os
    from core import CCP4Modules
    from core import CCP4XtalData

    jobDir = CCP4Modules.PROJECTSMANAGER().jobDirectory(jobId=jobId,create=False)
    exportFile = os.path.join(jobDir,'exportMtz.mtz')
    if os.path.exists(exportFile): return exportFile

    childJobs = CCP4Modules.PROJECTSMANAGER().db().getChildJobs(jobId=jobId,details=True)
    print 'import_merged.exportMtz',childJobs
    truncateOut = None
    freerflagOut = None
    for jobNo,subJobId,taskName  in childJobs:
      if taskName == 'aimless_pipe':
         aimlessChildJobs = CCP4Modules.PROJECTSMANAGER().db().getChildJobs(jobId=subJobId,details=True)
         print 'import_merged.exportMtz aimlessChildJobs',aimlessChildJobs
         for jobNo0,subJobId0,taskName0  in aimlessChildJobs:
           if taskName0 == 'ctruncate':
             truncateOut = os.path.join( CCP4Modules.PROJECTSMANAGER().jobDirectory(jobId=subJobId0,create=False),'HKLOUT.mtz')
             if not os.path.exists(truncateOut): truncateOut = None
     
    freerflagOut = os.path.join( CCP4Modules.PROJECTSMANAGER().jobDirectory(jobId=jobId,create=False),'FREEOUT.mtz')
    if not os.path.exists(freerflagOut): freerflagOut = None
    if truncateOut is None: return None
    if freerflagOut is None: return truncateOut

    print 'aimless_pipe.exportJobFile  runCad:',exportFile,[ freerflagOut ]
    

    m = CCP4XtalData.CMtzDataFile(truncateOut)
    #print m.runCad.__doc__   #Print out docs for the function
    outfile,err = m.runCad(exportFile,[ freerflagOut ] )
    print 'aimless_pipe.exportJobFile',outfile,err.report()
    return   outfile                                                   
 
def exportJobFileMenu(jobId=None):
    # Return a list of items to appear on the 'Export' menu - each has three subitems:
    # [ unique identifier - will be mode argument to exportJobFile() , menu item , mime type (see CCP4CustomMimeTypes module) ]
    return [ [ 'complete_mtz' ,'MTZ file' , 'application/CCP4-mtz' ] ]
                                                
"""
