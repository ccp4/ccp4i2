"""
     CCP4CootData.py: CCP4 GUI Project
     Copyright (C) 2013 STFC

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

from core import CCP4Data
from core import CCP4File
from core import CCP4ErrorHandling

class CCootHistoryDataFile(CCP4File.CDataFile):

    SUBTYPE_INITIAL = 1
    SUBTYPE_HISTORY = 2
    CONTENTS = {}
    CONTENTS.update(CCP4File.CDataFile.CONTENTS)
    CONTENTS['subType'] = {'class' : CCP4Data.CInt,
                           'qualifiers' : {'default' : SUBTYPE_HISTORY, 'enumerators' : [SUBTYPE_INITIAL, SUBTYPE_HISTORY],
                                           'onlyEnumerators':True, 'menuText' : ['Coot 0-state.scm', 'Coot history.scm']}}
    QUALIFIERS = {'mimeTypeName' : 'application/coot-script', 'mimeTypeDescription' : 'Coot history/script file',
                  'fileExtensions' : ['scm', 'py'], 'fileContentClassName' : None,
                  'guiLabel' : 'Coot history', 'fileLabel' : 'coot_history',
                  'toolTip' : "history.scm or 0-state.scm file from Coot" }

    def assertSame(self, arg, testPath=False, testChecksum=True, testSize=False, testDiff=False, diagnostic=False, fileName=None):
        #MN Now here we have an issue that assertSame on an Coot Hstory file is fraught with difficulties, and an identical checkSum is probably far too stringent.  I'm gonna suggest that we should remove testChecksum for now, with a view to putting in a more intelligent comparison later
        report = CCP4File.CDataFile.assertSame(self, arg, testPath, False, testSize, testDiff, diagnostic, fileName)
        return report


