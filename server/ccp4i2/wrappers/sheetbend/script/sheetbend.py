import os
import pathlib

from lxml import etree

from ccp4i2.core.CCP4PluginScript import CPluginScript


class sheetbend(CPluginScript):
    TASKNAME = 'sheetbend'   # Task name - should be same as class name and match pluginTitle in the .def.xml file
    TASKVERSION= 0.1               # Version of this plugin
    MAINTAINER = 'kathryn.cowtan@york.ac.uk'
    TASKCOMMAND="csheetbend"
    
    def __init__(self, *args, **kws):
        super(sheetbend, self).__init__(*args, **kws)

    def processInputFiles(self):
        from ccp4i2.core import CCP4XtalData
        colgrps = [ ['F_SIGF', CCP4XtalData.CObsDataFile.CONTENT_FLAG_FMEAN] ]
        if self.container.inputData.FREERFLAG.isSet(): colgrps.append( 'FREERFLAG' )
        self.hklin, columns, error = self.makeHklin0( colgrps )
        from ccp4i2.core import CCP4ErrorHandling
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

    def makeCommandAndScript(self):
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

