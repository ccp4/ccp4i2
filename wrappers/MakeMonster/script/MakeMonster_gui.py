"""
     MakeMonster task widget
"""

from PySide6 import QtGui, QtWidgets,QtCore

from qtgui.CCP4TaskWidget import CTaskWidget
from core import CCP4Modules

#-------------------------------------------------------------------
class MakeMonster_gui(CTaskWidget):
#-------------------------------------------------------------------

# Subclass CTaskWidget to give specific task window
    TASKNAME = 'MakeMonster'
    TASKVERSION = 0.1
    TASKMODULE='developer_tools'
    TASKTITLE='Export monster mtz'
    DESCRIPTION='Build an MTZ file from data objects'
    
    dataTypesDict = {'Obs':'reflection','Phi':'phase','FWT':'Map coeffs','DELFWT':'Difmap coeffs','FREER':'FreeR flags'}
    
    def drawContents(self):
        import functools
        self.openFolder(folderFunction='inputData')
        
        for dataType in ['Obs','Phi','FWT','DELFWT','FREER']:
            countObjectName = 'N'+dataType+'Objects'
            self.createLine( [ 'subtitle', 'Input '+self.dataTypesDict[dataType]+' data objects','label','Number to include','widget', countObjectName ])
            countObject = getattr(self.container.inputData, countObjectName)
            countObject.dataChanged.connect( functools.partial(self.handleNObjectsChanged,dataType))
            
            for i in range(6):
                dataObjectName = dataType+"_"+str(i+1)
                dataObject = getattr(self.container.inputData,dataObjectName)
                outputAsName = dataObjectName + "_OutputRepresentation"
                SuffixName = dataObjectName + "_Suffix"
                toggleValues = [6-j for j in range(6-i)]
                self.openSubFrame(toggle = [countObjectName,'open', toggleValues])
                self.createLine( [ 'widget', dataObjectName])
                dataObject.dataChanged.connect( functools.partial(self.handleObjectChanged,dataObjectName))
                self.handleObjectChanged(dataObjectName)
                if hasattr(self.container.inputData, outputAsName):
                    self.createLine( ['label','Output as: ','widget',outputAsName,'label',' with column Suffix', 'widget',SuffixName] )
                else:
                    self.createLine( ['label','Output column Suffix', 'widget',SuffixName] )
                self.closeSubFrame()
            
            self.handleNObjectsChanged(dataType)
    
        self.validate()

    @QtCore.Slot(str)
    def handleNObjectsChanged(self, dataType):
        countObjectName = 'N'+dataType+'Objects'
        newNObjects = getattr(self.container.inputData,countObjectName).__int__()
        for i in range(newNObjects):
            dataObjectName = dataType+"_"+str(i+1)
            outputAsName = dataObjectName + "_OutputRepresentation"
            SuffixName = dataObjectName + "_Suffix"
            getattr(self.container.inputData,dataObjectName).setQualifiers({'allowUndefined':False})
            if newNObjects > 1:
                getattr(self.container.inputData,SuffixName).setQualifiers({'allowUndefined':False})
            if hasattr(self.container.inputData, outputAsName):
                getattr(self.container.inputData,outputAsName).setQualifiers({'allowUndefined':False})
        for i in range(getattr(self.container.inputData,countObjectName).__int__(),6):
            dataObjectName = dataType+"_"+str(i+1)
            outputAsName = dataObjectName + "_OutputRepresentation"
            SuffixName = dataObjectName + "_Suffix"
            getattr(self.container.inputData,dataObjectName).setQualifiers({'allowUndefined':True})
            getattr(self.container.inputData,SuffixName).setQualifiers({'allowUndefined':True})
            if hasattr(self.container.inputData, outputAsName):
                getattr(self.container.inputData,outputAsName).setQualifiers({'allowUndefined':True})
            getattr(self.container.inputData,dataObjectName).unSet()
        self.validate()

    @QtCore.Slot(str)
    def handleObjectChanged(self, dataObjectName):
        if dataObjectName.startswith('Obs'):
            dataObject = self.container.inputData.Obs_1#getattr(self.container.inputData, dataObjectName)._value_
            dataObject.setContentFlag(reset=True)
            outputAsName = dataObjectName + "_OutputRepresentation"
            outputAs = getattr(self.container.inputData,outputAsName)
            try:
                if dataObject.contentFlag.__int__() == 1:
                    outputAs.set(1)
                    outputAs.setQualifiers({'enumerators':[1,2,3,4],'menuText':['I+/I-','F+/F-','Imean','Fmean']})
                elif dataObject.contentFlag.__int__() == 2:
                    outputAs.set(2)
                    outputAs.setQualifiers({'enumerators':[2,4],'menuText':['F+/F-','Fmean']})
                elif dataObject.contentFlag.__int__() == 3:
                    outputAs.set(3)
                    outputAs.setQualifiers({'enumerators':[3,4],'menuText':['Imean','Fmean']})
                elif dataObject.contentFlag.__int__() == 4:
                    outputAs.set(4)
                    outputAs.setQualifiers({'enumerators':[4],'menuText':['Fmean']})
                self.getWidget(outputAsName).populateComboBox(outputAs)
                self.getWidget(outputAsName).updateViewFromModel()
            except:
                pass
