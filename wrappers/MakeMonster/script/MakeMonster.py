
from core.CCP4PluginScript import CPluginScript
from PySide2 import QtCore
import os,re,time,sys

class MakeMonster(CPluginScript):
    
    TASKNAME = 'MakeMonster'                                  # Task name - should be same as class name
    TASKMODULE='test'
    TASKVERSION= 0.0                                     # Version of this plugin
    ASYNCHRONOUS = False
    TIMEOUT_PERIOD = 9999999.9
    MAINTAINER = 'martin.noble@newcastle.ac.uk'
    RUNEXTERNALPROCESS=False

    def makeCommandAndScript(self):
        return CPluginScript.SUCCEEDED
    
    def startProcess(self, command):
        return CPluginScript.SUCCEEDED

    def processInputFiles(self):
        inputDataObjects = []
        colOutList = []
        for dataType in ['Obs','FREER','Phi','FWT','DELFWT']:
            for i in range(6):
                dataObjectName = dataType+'_'+str(i+1)
                dataObject = getattr(self.container.inputData,dataObjectName)
                if dataObject.isSet():
                    inputDataObject = [dataObjectName]
                    outputAsName = dataObjectName + "_OutputRepresentation"
                    suffixName = dataObjectName + "_Suffix"
                    suffixObject = getattr(self.container.inputData,suffixName)
                    
                    representationNumber = 1
                    try:
                        representationNumber = getattr(self.container.inputData,outputAsName).__int__()
                    except Exception as err:
                        pass
                    inputDataObject.append(representationNumber)
                    
                    inputDataObjects.append(inputDataObject)
                    
                    for columnLabel in dataObject.CONTENT_SIGNATURE_LIST[representationNumber-1]:
                        patchedColumnLabel = columnLabel
                        if dataType == 'Obs':
                            patchedColumnLabel = patchedColumnLabel.replace('plus','(+)').replace('minus','(-)')
                        elif dataType == 'FWT':
                            patchedColumnLabel = patchedColumnLabel.replace('F','FWT').replace('PHI','PHWT')
                        elif dataType == 'DELFWT':
                            patchedColumnLabel = patchedColumnLabel.replace('F','DELFWT').replace('PHI','PHDELWT')
                        if len(suffixObject.__str__().strip()) > 0:
                            patchedColumnLabel += ('_'+suffixObject.__str__().strip())
                        print({'suffixObject':suffixObject.__str__(),'columnLabel':columnLabel,'patchedColumnLabel':patchedColumnLabel,'dataType':dataType})
                        colOutList.append(patchedColumnLabel)
        colOutListAsString = ",".join(colOutList)
        self.intermediateFile1, colInListAsString, errReport = self.makeHklin0(inputDataObjects)
        print({'colOutList':colOutListAsString, 'colInList':colInListAsString})
        
        return CPluginScript.SUCCEEDED

    def processOutputFiles(self):
        os.rename(self.intermediateFile1,self.container.outputData.HKLOUT.__str__())
        with open(self.makeFileName('PROGRAMXML'),'w') as programXML:
            programXML.write('<MakeMonster></MakeMonster>')
        return CPluginScript.SUCCEEDED
