"""
    nautilus_build_refine.py: CCP4 GUI Project
    
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

import os, shutil
#from lxml import etree
from xml.etree import ElementTree as ET
from core import CCP4Utils
from core.CCP4PluginScript import CPluginScript

class nautilus_build_refine(CPluginScript):

    TASKMODULE="model_building"
    TASKTITLE = 'Autobuild RNA/DNA (Nautilus pipeline)'
    SHORTTASKTITLE='NAUTILUS'
    TASKNAME = 'nautilus_build_refine'
    TASKVERSION= 0.1
    MAINTAINER = 'kevin.cowtan@york.ac.uk'
    ERROR_CODES = { 201 : {'description' : 'Failed to analyse output files' } }
    PURGESEARCHLIST = [ ]
    RUNEXTERNALPROCESS = False
    WHATNEXT = ['coot_rebuild','prosmart_refmac','nautilus_build_refine']


    def __init__(self, *args, **kws):
        super(nautilus_build_refine, self).__init__(*args, **kws)

    #The startProcess method is where you build in the pipeline logic
    def startProcess(self, command, **kws):

        # XML
        self.pipelinexmlfile = self.makeFileName(format='PROGRAMXML')
        self.xmlroot = ET.Element("NautilusBuildRefineResult")
        self.xmlcyc = ET.SubElement(self.xmlroot,"BuildRefineCycle")
        ET.SubElement(self.xmlcyc,"Number").text = str(1)

        self.ncycles = int( self.container.controlParameters.ITERATIONS )

        # FIRST CYCLE

        #  NAUTILUS
        nautilusPlugin = self.makePluginObject("nautilus")
        for attr in ("F_SIGF","ABCD","FREERFLAG","ASUIN"):
          setattr(nautilusPlugin.container.inputData,attr,getattr(self.container.inputData,attr))
        if self.container.controlParameters.XYZIN_MODE:
          nautilusPlugin.container.inputData.XYZIN = self.container.inputData.XYZIN
        nautilusPlugin.container.controlParameters.CYCLES = self.container.controlParameters.NAUTILUS_CYCLES
        nautilusPlugin.container.controlParameters.RESOLUTION = self.container.controlParameters.NAUTILUS_RESOLUTION
        nautilusPlugin.container.controlParameters.ANISOTROPY_CORRECTION = self.container.controlParameters.NAUTILUS_ANISOTROPY_CORRECTION
        status = nautilusPlugin.process()
        self.currentCoordinates = nautilusPlugin.container.outputData.XYZOUT
        self.copyNautilusXML(nautilusPlugin)
        #  REFMAC
        refmacPlugin = self.makePluginObject('refmac')
        refmacPlugin.container.inputData.XYZIN.set( self.currentCoordinates )
        for attr in ("F_SIGF","ABCD","FREERFLAG","DICT"):
          setattr(refmacPlugin.container.inputData,attr,getattr(self.container.inputData,attr))
        refmacPlugin.container.controlParameters.NCYCLES = self.container.controlParameters.REFMAC_CYCLES
        refmacPlugin.container.controlParameters.USE_LOCAL_SYMMETRY = self.container.controlParameters.REFMAC_LOCAL_NCS
        refmacPlugin.container.controlParameters.HYDROGENS = 'NO'
        refmacPlugin.container.controlParameters.PHOUT = True
        status = refmacPlugin.process()
        self.currentCoordinates = refmacPlugin.container.outputData.XYZOUT
        self.currentMap = refmacPlugin.container.outputData.FPHIOUT
        self.currentABCD = refmacPlugin.container.outputData.ABCDOUT
        self.copyRefmacXML(refmacPlugin)

        # SUBSEQUENT CYCLES
        for cyc in range(1,self.ncycles):
          self.xmlcyc = ET.SubElement(self.xmlroot,"BuildRefineCycle")
          ET.SubElement(self.xmlcyc,"Number").text = str(cyc+1)

          #  NAUTILUS
          nautilusPlugin = self.makePluginObject("nautilus")
          for attr in ("F_SIGF","FREERFLAG","ASUIN"):
            setattr(nautilusPlugin.container.inputData,attr,getattr(self.container.inputData,attr))
          nautilusPlugin.container.controlParameters.CYCLES = self.container.controlParameters.NAUTILUS_CYCLES
          nautilusPlugin.container.controlParameters.RESOLUTION = self.container.controlParameters.NAUTILUS_RESOLUTION
          nautilusPlugin.container.controlParameters.ANISOTROPY_CORRECTION = self.container.controlParameters.NAUTILUS_ANISOTROPY_CORRECTION
          nautilusPlugin.container.inputData.XYZIN.set( self.currentCoordinates )
          nautilusPlugin.container.inputData.FWT_PHWT_IN.set( self.currentMap )
          nautilusPlugin.container.inputData.ABCD.set( self.currentABCD )
          status = nautilusPlugin.process()
          self.currentCoordinates = nautilusPlugin.container.outputData.XYZOUT
          self.copyNautilusXML(nautilusPlugin)
          #  REFMAC
          refmacPlugin = self.makePluginObject('refmac')
          refmacPlugin.container.inputData.XYZIN.set( self.currentCoordinates )
          for attr in ("F_SIGF","ABCD","FREERFLAG","DICT"):
            setattr(refmacPlugin.container.inputData,attr,getattr(self.container.inputData,attr))
          refmacPlugin.container.controlParameters.NCYCLES = self.container.controlParameters.REFMAC_CYCLES
          refmacPlugin.container.controlParameters.USE_LOCAL_SYMMETRY = self.container.controlParameters.REFMAC_LOCAL_NCS
          refmacPlugin.container.controlParameters.HYDROGENS = 'NO'
          refmacPlugin.container.controlParameters.PHOUT = True
          status = refmacPlugin.process()
          self.currentCoordinates = refmacPlugin.container.outputData.XYZOUT
          self.currentMap = refmacPlugin.container.outputData.FPHIOUT
          self.currentABCD = refmacPlugin.container.outputData.ABCDOUT
          self.currentDIFFPHI = refmacPlugin.container.outputData.DIFFPHIOUT
          self.copyRefmacXML(refmacPlugin)

        # DONE
        return CPluginScript.SUCCEEDED


    def processOutputFiles(self):
        shutil.copyfile( str(self.currentCoordinates), str(self.container.outputData.XYZOUT.fullPath) )
        self.container.outputData.XYZOUT.annotation.set('Model built by '+self.TASKTITLE)
        shutil.copyfile( str(self.currentMap), str(self.container.outputData.FPHIOUT.fullPath) )
        self.container.outputData.FPHIOUT.annotation.set('Map generated by '+self.TASKTITLE)
        shutil.copyfile( str(self.currentDIFFPHI), str(self.container.outputData.DIFFPHIOUT.fullPath) )
        self.container.outputData.DIFFPHIOUT.annotation.set('Difference map generated by '+self.TASKTITLE)
        shutil.copyfile( str(self.currentABCD), str(self.container.outputData.ABCDOUT.fullPath) )
        self.container.outputData.ABCDOUT.annotation.set('Phases generated by '+self.TASKTITLE)
        return CPluginScript.SUCCEEDED


    # Utility functions
    def copyNautilusXML(self,plugin):
        self.xmlcyc.append( CCP4Utils.openFileToEtree(plugin.makeFileName('PROGRAMXML')) )
        f = open( self.pipelinexmlfile,'w')
        CCP4Utils.writeXML(f,ET.tostring(self.xmlroot))

    def copyRefmacXML(self,plugin):
        rxml = CCP4Utils.openFileToEtree(plugin.makeFileName('PROGRAMXML'))
        rstats = rxml.findall(".//REFMAC/Overall_stats/stats_vs_cycle")
        if len(rstats)>0:
          refele = ET.SubElement(self.xmlcyc,'RefmacResult')
          for node in rstats[0].findall("new_cycle[last()]/r_factor | new_cycle[last()]/r_free | new_cycle[last()]/rmsBOND |  new_cycle[last()]/rmsANGLE"):
            node.text = str(node.text).strip()
            if node.tag == 'rmsBOND':
              node.text = str(100*float(node.text))
              node.tag='rmsBONDx100'
            refele.append( node )
        f = open( self.pipelinexmlfile,'w')
        CCP4Utils.writeXML(f,ET.tostring(self.xmlroot))

