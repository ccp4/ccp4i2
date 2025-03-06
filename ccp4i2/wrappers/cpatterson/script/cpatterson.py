"""
     cpatterson.py: CCP4 GUI 2 Project
     Copyright (C) 2014 The University of York, 2024 STFC

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

from lxml import etree
from core.CCP4PluginScript import CPluginScript
from core import CCP4Utils

class cpatterson(CPluginScript):

    TASKMODULE = 'wrappers'
    TASKTITLE = 'Prepare map coefficients'
    TASKNAME = 'cpatterson'
    TASKCOMMAND = 'cpatterson'
    TASKVERSION= 0.0
    MAINTAINER = 'stuart.mcnicholas@york.ac.uk'

    def processInputFiles ( self ):
        from core import CCP4XtalData
        self.hklin,error = self.makeHklin([['F_SIGF',CCP4XtalData.CObsDataFile.CONTENT_FLAG_FMEAN]])

    def makeCommandAndScript(self):

        out = self.container.outputData

        self.appendCommandLine([ '-stdin' ])

        self.appendCommandScript([ '-mtzin ' + self.hklin ])
        self.appendCommandScript([ '-mapout %s'%(str(out.MAPOUT.fullPath)) ])
        self.appendCommandScript ([ '-colin-fo F,SIGF' ])

        return CPluginScript.SUCCEEDED

    def processOutputFiles(self):

        self.container.outputData.MAPOUT.annotation = 'Computed using ' + str(self.container.inputData.F_SIGF.annotation)

        xmlRoot = etree.Element('Patterson')
        with open ( self.makeFileName('PROGRAMXML'),'w' ) as xmlFile:
            xmlString = etree.tostring ( xmlRoot, pretty_print=True )
            CCP4Utils.writeXML(xmlFile,xmlString)

        return CPluginScript.SUCCEEDED
