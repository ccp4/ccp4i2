import os

from ccp4i2.core.CCP4PluginScript import CPluginScript


class MakeMonster(CPluginScript):
    TASKNAME = 'MakeMonster'                                  # Task name - should be same as class name
    TASKVERSION= 0.0                                     # Version of this plugin
    TIMEOUT_PERIOD = 9999999.9
    MAINTAINER = 'martin.noble@newcastle.ac.uk'

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
                    if hasattr(self.container.inputData,outputAsName):
                        representationNumber = getattr(self.container.inputData,outputAsName).__int__()
                    inputDataObject.append(representationNumber)
                    
                    inputDataObjects.append(inputDataObject)
                    
                    for columnLabel in dataObject.CONTENT_SIGNATURE_LIST[representationNumber-1]:
                        patchedColumnLabel = columnLabel
                        if dataType == 'Obs':
                            patchedColumnLabel = patchedColumnLabel.replace('plus','(+)').replace('minus','(-)')
                            pass
                        elif dataType == 'FWT':
                            patchedColumnLabel = patchedColumnLabel.replace('F','FWT').replace('PHI','PHWT')
                        elif dataType == 'DELFWT':
                            patchedColumnLabel = patchedColumnLabel.replace('F','DELFWT').replace('PHI','PHDELWT')
                        
                        if len(suffixObject.__str__().strip()) > 0:
                            patchedColumnLabel += ('_'+suffixObject.__str__().strip())
                        colOutList.append(patchedColumnLabel)
        self.intermediateFile1, colInListAsString, errReport = self.makeHklin0(inputDataObjects)
        
        return CPluginScript.SUCCEEDED

    def processOutputFiles(self):
        os.rename(self.intermediateFile1,self.container.outputData.HKLOUT.__str__())
        with open(self.makeFileName('PROGRAMXML'),'w') as programXML:
            programXML.write('<MakeMonster></MakeMonster>')
        return CPluginScript.SUCCEEDED
