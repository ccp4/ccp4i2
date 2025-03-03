"""
    sheetbend.py: CCP4 GUI Project
    
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

import os
from core.CCP4PluginScript import CPluginScript
from lxml import etree
import pathlib

class sheetbend(CPluginScript):
    TASKNAME = 'sheetbend'   # Task name - should be same as class name and match pluginTitle in the .def.xml file
    TASKVERSION= 0.1               # Version of this plugin
    MAINTAINER = 'kevin.cowtan@york.ac.uk'
    ERROR_CODES = { 201 : {'description' : 'Failed to analyse output files' },
                    202 : {'description' : 'Failed applying selection ot PDB file' }
                    }
    PURGESEARCHLIST = [ [ 'hklin.mtz' , 0 ], ['log_mtzjoin.txt', 0] ]
    TASKCOMMAND="csheetbend"
    
    def __init__(self, *args, **kws):
        super(sheetbend, self).__init__(*args, **kws)

    def processInputFiles(self):
        from core import CCP4XtalData
        colgrps = [ ['F_SIGF', CCP4XtalData.CObsDataFile.CONTENT_FLAG_FMEAN] ]
        if self.container.inputData.FREERFLAG.isSet(): colgrps.append( 'FREERFLAG' )
        self.hklin, columns, error = self.makeHklin0( colgrps )
        from core import CCP4ErrorHandling
        if error.maxSeverity()>CCP4ErrorHandling.SEVERITY_WARNING:
            return CPluginScript.FAILED
        #Preprocess coordinates to extract a subset
        self.selectedCoordinatesPath = os.path.join(self.getWorkDirectory(), "selected_xyzin.pdb")
        self.container.inputData.XYZIN.loadFile()
        if self.container.inputData.XYZIN.isMMCIF():
            self.selectedCoordinatesPath = str(pathlib.Path(self.selectedCoordinatesPath).with_suffix('.cif'))
        oldFullPath = pathlib.Path(self.container.outputData.XYZOUT.fullPath.__str__())
        if self.container.inputData.XYZIN.isMMCIF():
            self.container.outputData.XYZOUT.setFullPath(str(oldFullPath.with_suffix('.cif')))

        self.container.inputData.XYZIN.getSelectedAtomsPdbFile(self.selectedCoordinatesPath)
        return CPluginScript.SUCCEEDED

    def makeCommandAndScript(self,**kw):
        self.appendCommandLine([ '-stdin' ])
        self.appendCommandScript( 'mtzin '  + self.hklin )
        self.appendCommandScript( 'pdbin '  + self.selectedCoordinatesPath )
        self.appendCommandScript( 'pdbout %s'%(str(self.container.outputData.XYZOUT.fullPath)) )
        self.appendCommandScript( "colin-fo F_SIGF_F,F_SIGF_SIGF" )
        if self.container.inputData.FREERFLAG.isSet():          self.appendCommandScript( "colin-free FREERFLAG_FREER" )
        if ( self.container.controlParameters.REFINE_COORD   ): self.appendCommandScript( "coord" )
        if ( self.container.controlParameters.REFINE_U_ISO   ): self.appendCommandScript( "u-iso" )
        if ( self.container.controlParameters.REFINE_U_ANISO ): self.appendCommandScript( "u-aniso" )
        if ( self.container.controlParameters.POSTREFINE_COORD   ): self.appendCommandScript( "postrefine-coord" )
        if ( self.container.controlParameters.POSTREFINE_U_ISO   ): self.appendCommandScript( "postrefine-u-iso" )
        if ( self.container.controlParameters.POSTREFINE_U_ANISO ): self.appendCommandScript( "postrefine-u-aniso" )
        if ( self.container.controlParameters.PSEUDO_REGULARIZE  ): self.appendCommandScript( "pseudo-regularize" )
        self.appendCommandScript( "cycles %s"%(str(self.container.controlParameters.CYCLES)) )
        self.appendCommandScript( "refine-regularize-cycles %s"%(str(self.container.controlParameters.REFINE_REGULARIZE_CYCLES)) )
        self.appendCommandScript( "resolution-by-cycle %s"%(str(self.container.controlParameters.RESOLUTION)) )
        self.appendCommandScript( "radius-scale %s"%(str(self.container.controlParameters.RADIUS_SCALE)) )
        self.appendCommandScript( "xmlout %s"%(self.makeFileName('PROGRAMXML')) )
        return CPluginScript.SUCCEEDED

    def processOutputFiles(self):        
        self.xmlroot = etree.parse(self.makeFileName('PROGRAMXML')).getroot()
        return CPluginScript.SUCCEEDED

