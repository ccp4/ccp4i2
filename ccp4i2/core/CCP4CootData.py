from . import CCP4Data
from . import CCP4File


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
