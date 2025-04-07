"""
Copyright (C) 2014 Newcastle University
Author: Martin Noble
"""

import os
import shutil

from ....core.CCP4Modules import PREFERENCES, PROJECTSMANAGER
from ....qtgui.CCP4TaskWidget import CTaskWidget


def whatNext(jobId=None,childTaskName=None,childJobNumber=None,projectName=None):
    from ....core import CCP4File
    from ....core import CCP4Container
    from ....core import CCP4Data
    try:
        jobDirectory = PROJECTSMANAGER().db().jobDirectory(jobId=jobId)
        returnList = []
        
        for taskName in ['ShelxCE','ShelxCECompareHands']:
            cont = PROJECTSMANAGER().getJobParams(jobId)
            fileRoot = 'Subsequent'+taskName+'.params.xml'
            paramsPath = os.path.normpath(os.path.join(jobDirectory,fileRoot))
            if not os.path.isfile(paramsPath):
                cont_new = CCP4Container.CContainer(name=taskName)
                cont_new.__dict__['header']=CCP4File.CI2XmlHeader()
                cont_new.__dict__['header'].pluginName = taskName
                for ic in ('inputData','controlParameters'):
                    cont_new.addObject( getattr(cont,ic), reparent=True )
                cont_new.inputData.addObject( cont.outputData.XYZOUT, reparent=True )
                cont_new.inputData.renameObject('XYZOUT','HAIN')
                cont_new.addObject(CCP4Container.CContainer(parent=cont_new,name='keywords'),name='keywords')
                cont_new.keywords.addContent( name='z', cls=CCP4Data.CBoolean, value=True )
                cont_new.saveDataToXml(fileName=paramsPath)
            returnList.append([taskName, paramsPath])
    
        if cont.controlParameters.MODE.__str__() == 'SAD':
            cont = PROJECTSMANAGER().getJobParams(jobId)
            taskName='phaser_EP_AUTO'
            fileRoot = 'Subsequent'+taskName+'.params.xml'
            paramsPath = os.path.normpath(os.path.join(jobDirectory,fileRoot))
            if not os.path.isfile(paramsPath):
                cont_new = CCP4Container.CContainer(name=taskName)
                cont_new.__dict__['header']=CCP4File.CI2XmlHeader()
                cont_new.__dict__['header'].pluginName = taskName
                cont_new.addParamsSubContainers()
                cont_new.inputData.addObject( cont.outputData.XYZOUT, reparent=True )
                cont_new.inputData.renameObject('XYZOUT','XYZIN_HA')
                cont_new.inputData.addObject( cont.inputData.SAD, reparent=True )
                cont_new.inputData.renameObject('SAD','F_SIGF')
                cont_new.inputData.addContent( name='LLGC_CYCLES', cls=CCP4Data.CInt, value=10 )
                cont_new.inputData.addContent( name='ELEMENTS', cls=CCP4Data.CList, subItem={'class':CCP4Data.CString} )
                cont_new.inputData.ELEMENTS.addItem(value=cont.controlParameters.SFAC.__str__())
                cont_new.addObject(CCP4Container.CContainer(parent=cont_new,name='keywords'),name='keywords')
                cont_new.keywords.addContent( name='HAND', cls=CCP4Data.CString, value='Off' )
                cont_new.saveDataToXml(fileName=paramsPath)
            returnList.append([taskName, paramsPath])

        return returnList

    except:
        print('Exception in ShelxCD_gui.whatNext')
        return ['ShelxCECompareHands','ShelxCE','phaser_EP_AUTO']

#-------------------------------------------------------------------
class ShelxCD_gui(CTaskWidget):
    #-------------------------------------------------------------------
    
    # Subclass CTaskWidget to give specific task window
    TASKMODULE = 'test'                               # Where this plugin will appear on the gui
    TASKTITLE = 'Find HA sites - SHELXC/D'     # A short title for gui menu
    SHORTTASKTITLE = 'ShelxCD'     # No - this one really is a short title for gui menu
    TASKVERSION = 0.1
    DESCRIPTION = 'Find sites from SAD/MAD/SIR/SIRAS/RIP/RIPAS data'
    TASKNAME = 'ShelxCD'                                  # Task name - should be same as class name
    MGDISPLAYFILES = ['XYZOUT']

    def __init__(self,*args,**kw):
        super(ShelxCD_gui,self).__init__(*args, **kw)

    def drawContents(self):
        
        self.openFolder(folderFunction='inputData')

        if (not hasattr(PREFERENCES(),'SHELXDIR')) and shutil.which('shelxc') is None:
            if (not PREFERENCES().SHELXDIR.exists()) and shutil.which('shelxc') is None:
              self.createLine ( [ 'warning','The Shelx programs have not been found. They are not part of CCP4 but you can get them from\nhttp://shelx.uni-ac.gwdg.de/SHELX/download.php\nIf you already have them make sure they are on the search path\nOR specify where they are in the Preferences window - under Other Software.' ])
        
        self.createLine ( [ 'label','Experiment type','stretch','tip','What sort of data  are available for finding heavy atoms','widget','MODE' ] )
        self.openSubFrame( frame=[True])
        self.createLine(['label','Atom type to seek','stretch','widget','SFAC'])
        self.createLine(['label','Number to find','stretch','widget','FIND'])
        self.createLine(['label','Number of trys','stretch','widget','NTRY'])
        self.closeSubFrame()
        self.openSubFrame( frame=True)
        
        self.createLine(['advice','Anomalous (SAD) dataset'],toggle=['MODE','open', ['SAD']])
        self.createLine(['widget','SAD'],toggle=['MODE','open', ['SAD']])

        self.openSubFrame( toggle=['MODE','open',['MAD']])
        self.createLine(['advice','High energy remote dataset'])
        self.createLine(['widget','HREM'])
        self.createLine(['advice','Low energy remote dataset'])
        self.createLine(['widget','LREM'])
        self.createLine(['advice','Inflection point dataset'])
        self.createLine(['widget','INFL'])
        self.createLine(['advice','Peak dataset'])
        self.createLine(['widget','PEAK'])
        self.closeSubFrame()

        self.createLine(['advice','Derivative dataset'],toggle=['MODE','open',['SIR']])
        self.createLine(['widget','SIR'],toggle=['MODE','open',['SIR']])

        self.createLine(['advice','Anomalous derivative dataset'], toggle=['MODE','open',['SIRAS']])
        self.createLine(['widget','SIRA'], toggle=['MODE','open',['SIRAS']])
        
        self.createLine(['advice','Native dataset'])
        self.createLine(['widget','NAT'])

        self.createLine(['advice','Radiation damaged dataset'], toggle=['MODE','open',['RIP']])
        self.createLine(['widget','RIP'], toggle=['MODE','open',['RIP']])

        self.createLine(['advice','Radiation damaged anomalous dataset'], toggle=['MODE','open',['RIPAS']])
        self.createLine(['widget','RIPA'], toggle=['MODE','open',['RIPAS']])
        self.closeSubFrame()
    
        self.dShelxCDKeywords()

    def dShelxCDKeywords(self):
        #From here on, GUI is autogenerated from the keywords container
        self.openFolder(folderFunction='controlParameters', title='Keywords')
        self.openSubFrame()
        self.autoGenerate(container=self.container.controlParameters,selection={'excludeParameters' : ['MODE','FIND','NTRY','SFAC']})
        self.closeSubFrame()

    def isValid(self):
        inp = self.container.inputData
        requires = {
            'SAD':[inp.SAD],
            'MAD':[inp.HREM,inp.LREM,inp.INFL,inp.PEAK],
            'SIR':[inp.SIR,inp.NAT],
            'SIRAS':[inp.SIRA,inp.NAT],
            'RIP':[inp.RIP,inp.NAT],
            'RIPAS':[inp.RIPA,inp.NAT],}
        modeStr = self.container.controlParameters.MODE.__str__()
        # Need to unset the unused inputs before calling the CTaskWidget.isValid() or invalid data
        # will be flagged to user
        for fObj in [inp.SAD,inp.HREM,inp.LREM,inp.INFL,inp.PEAK,inp.SIR,inp.SIRA,inp.RIP,inp.RIPA]:
          if fObj not in requires[modeStr]: fObj.unSet()
        # Generic validity test
        invalidElements = super(ShelxCD_gui,self).isValid()

        if modeStr == 'MAD':
            nset = len([setDataset for setDataset in requires['MAD'] if setDataset.isSet()])
            if nset < 2:
                for require in requires[modeStr]:
                    # Check whether all of the datasets required by this mode are set
                    if not require.isSet() and not require in invalidElements:
                        invalidElements.append(require)        
        else:
            for require in requires[modeStr]:
                # Check whether all of the datasets required by this mode are set
                if not require.isSet() and not require in invalidElements:
                    invalidElements.append(require)
        
        return invalidElements
