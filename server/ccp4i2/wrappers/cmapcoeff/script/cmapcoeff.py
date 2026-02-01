import os
import xml.etree.ElementTree as ET

from ccp4i2.core import CCP4ErrorHandling, CCP4Utils, CCP4XtalData
from ccp4i2.core.CCP4PluginScript import CPluginScript


class cmapcoeff(CPluginScript):

    TASKMODULE = 'wrappers'
    TASKTITLE = 'Prepare map coefficients'
    TASKNAME = 'cmapcoeff'
    TASKCOMMAND = 'cmapcoeff'
    WHATNEXT = [ 'coot_rebuild' ]
    MAINTAINER = 'jon.agirre@york.ac.uk'

    def processInputFiles ( self ):
        list_of_stuff = [ ]

        inp = self.container.inputData
        con = self.container.controlParameters

        if con.MAPTYPE == 'anom' :
            list_of_stuff.append ( [ 'F_SIGF1', CCP4XtalData.CObsDataFile.CONTENT_FLAG_FPAIR ] )
            if inp.ABCD1.isSet ( ) :
                list_of_stuff.append ( [ 'ABCD1', CCP4XtalData.CPhsDataFile.CONTENT_FLAG_HL ] )
        elif con.MAPTYPE == 'fobs' :
            list_of_stuff.append ( [ 'F_SIGF1', CCP4XtalData.CObsDataFile.CONTENT_FLAG_FMEAN ] )
            list_of_stuff.append ( [ 'ABCD1', CCP4XtalData.CPhsDataFile.CONTENT_FLAG_HL ] )
        else :
            list_of_stuff.append ( [ 'F_SIGF1', CCP4XtalData.CObsDataFile.CONTENT_FLAG_FMEAN ] )
            if inp.ABCD1.isSet ( ) :
                list_of_stuff.append ( [ 'ABCD1', CCP4XtalData.CPhsDataFile.CONTENT_FLAG_HL ] )
            if inp.F_SIGF2.isSet ( ) :
                list_of_stuff.append ( [ 'F_SIGF2', CCP4XtalData.CObsDataFile.CONTENT_FLAG_FMEAN ] )
                if inp.ABCD2.isSet ( ) :
                    list_of_stuff.append ( [ 'ABCD2', CCP4XtalData.CPhsDataFile.CONTENT_FLAG_HL ] )

        self.hklin, columns, error = self.makeHklin0 ( list_of_stuff )

    def makeCommandAndScript(self):

        inp = self.container.inputData
        out = self.container.outputData
        con = self.container.controlParameters

        self.hklout = os.path.join(self.workDirectory,"hklout.mtz")

        self.appendCommandLine([ '-stdin' ])

        self.appendCommandScript( '-mtzin ' + self.hklin )
        self.appendCommandScript( '-mtzout ' + self.hklout )

        self.appendCommandScript( '-colout i2' )

        if con.MAPTYPE == 'anom' :
            self.appendCommandScript ([ '-colin-fano F_SIGF1_Fplus,F_SIGF1_SIGFplus,F_SIGF1_Fminus,F_SIGF1_SIGFminus' ])
            self.appendCommandScript ([ '-colin-hl-1 ABCD1_HLA,ABCD1_HLB,ABCD1_HLC,ABCD1_HLD' ])

        elif con.MAPTYPE == 'fobs' :
            self.appendCommandScript ([ '-colin-fo-1 F_SIGF1_F,F_SIGF1_SIGF' ])
            self.appendCommandScript ([ '-colin-hl-1 ABCD1_HLA,ABCD1_HLB,ABCD1_HLC,ABCD1_HLD' ])

        else:
            self.appendCommandScript ([ '-colin-fo-1 F_SIGF1_F,F_SIGF1_SIGF' ])
            self.appendCommandScript ([ '-colin-hl-1 ABCD1_HLA,ABCD1_HLB,ABCD1_HLC,ABCD1_HLD' ])

            self.appendCommandScript ([ '-colin-fo-2 F_SIGF2_F,F_SIGF2_SIGF' ])
            if inp.ABCD2.isSet ( ) :
                self.appendCommandScript ([ '-colin-hl-2 ABCD2_HLA,ABCD2_HLB,ABCD2_HLC,ABCD2_HLD' ])

        if inp.RESOLUTION.isSet():
            self.appendCommandScript ([ '-resolution %s'%(str(inp.RESOLUTION)) ])

        if con.B_VALUE.isSet():
            self.appendCommandScript ([ '-b-value %s'%(str(con.B_VALUE)) ])

        if con.U_VALUE.isSet():
            self.appendCommandScript ([ '-u-value %s'%(str(con.U_VALUE)) ])

        if con.SCALE :
            if con.F1_TO_F2 :
                self.appendCommandScript ([ '-scale-fo-1 ' ])
            else :
                self.appendCommandScript ([ '-scale-fo-2 ' ])

        return CPluginScript.SUCCEEDED

    def processOutputFiles(self):

        error = self.splitHklout ( [ 'FPHIOUT' ],[ 'i2.F_phi.F,i2.F_phi.phi' ] )

        if error.maxSeverity ( ) > CCP4ErrorHandling.SEVERITY_WARNING:
            return CPluginScript.FAILED

        if self.container.controlParameters.MAPTYPE == "anom" :
            self.container.outputData.FPHIOUT.annotation.set('Anomalous difference map coefficients')
            self.container.outputData.FPHIOUT.subType.set(CCP4XtalData.CMapCoeffsDataFile.SUBTYPE_ANOM_DIFFERENCE)
        elif self.container.controlParameters.MAPTYPE == "fobsfobs" :
            self.container.outputData.FPHIOUT.annotation.set('Fobs - Fobs difference map coefficients')
            self.container.outputData.FPHIOUT.subType.set(CCP4XtalData.CMapCoeffsDataFile.SUBTYPE_DIFFERENCE)
        else :
            self.container.outputData.FPHIOUT.annotation.set('Fobs weighted map coefficients')
            self.container.outputData.FPHIOUT.subType.set(CCP4XtalData.CMapCoeffsDataFile.SUBTYPE_NORMAL)

        if self.container.controlParameters.MAP_OUTPUT :    # get ready for producing a map file

            self.cfftPlugin = self.makeCfftPlugin ( )
            error = self.cfftPlugin.process ( )

            if error == CPluginScript.FAILED:
                self.reportStatus ( error )
            else :
                self.container.outputData.MAPOUT = self.cfftPlugin.container.outputData.MAPOUT
                self.container.outputData.MAPOUT.annotation.set ( 'Calculated using ' + str ( self.container.outputData.FPHIOUT.annotation ) )

        with open ( self.makeFileName('PROGRAMXML'),'w' ) as xmlFile:
            xmlRoot = ET.Element('cmapcoeff')
            xmlString = ET.tostring ( xmlRoot, pretty_print=True )
            CCP4Utils.writeXML(xmlFile,xmlString)

        return CPluginScript.SUCCEEDED


    def makeCfftPlugin ( self ):

        cfftPlugin = self.makePluginObject( 'fft' )
        cfftPlugin.container.inputData.FPHIIN = self.container.outputData.FPHIOUT
        cfftPlugin.container.inputData.RESOLUTION = self.container.inputData.RESOLUTION
        cfftPlugin.container.inputData.INDEX_U = self.container.inputData.INDEX_U
        cfftPlugin.container.inputData.INDEX_V = self.container.inputData.INDEX_V
        cfftPlugin.container.inputData.INDEX_W = self.container.inputData.INDEX_W

        return cfftPlugin
