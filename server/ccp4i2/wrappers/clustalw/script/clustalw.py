
from ccp4i2.core.CCP4PluginScript import CPluginScript
from ccp4i2.baselayer import QtCore
import os,glob,re,time,sys,shutil
from ccp4i2.core import CCP4XtalData
from lxml import etree
import math
from ccp4i2.core import CCP4Modules,CCP4Utils

class clustalw(CPluginScript):
    TASKTITLE = 'clustalw'     # A short title for gui menu
    DESCRIPTION = 'Perform multiple alignment'
    TASKNAME = 'clustalw'                                  # Task name - should be same as class name
    TASKVERSION= 0.0                                     # Version of this plugin
    ASYNCHRONOUS = False
    TIMEOUT_PERIOD =3.

    ERROR_CODES = {  200 : { 'description' : 'Failed to catenate sequences' },201 : { 'description' : 'Failed to setFullPath' },}
    
    def __init__(self,*args,**kws):
        CPluginScript.__init__(self, *args,**kws)
        self.TASKCOMMAND = os.path.normpath(os.path.join(CCP4Utils.getCCP4Dir(),'libexec','clustalw2'))
        self.xmlroot = etree.Element('Clustalw')

    def makeCommandAndScript(self):
        self.appendCommandLine('-infile='+self.startFileName)
        self.appendCommandLine('-type=protein')
        self.appendCommandLine('-outfile='+self.container.outputData.ALIGNMENTOUT.__str__())
        self.dndFilepath = os.path.normpath(os.path.join(self.getWorkDirectory(),'alignment.dnd'))
        self.appendCommandLine('-newtree='+self.dndFilepath)
        self.statsFilepath = os.path.normpath(os.path.join(self.getWorkDirectory(),'alignment.stats'))
        self.appendCommandLine('-stats='+self.statsFilepath)
        self.appendCommandLine('-align')
        return CPluginScript.SUCCEEDED
    
    def processInputFiles(self):
        if True:
            self.startFileName = os.path.normpath(os.path.join(self.getWorkDirectory(),'tempAlignment.fasta'))
            if self.container.inputData.SEQUENCELISTORALIGNMENT == 'ALIGNMENT':
                shutil.copyfile(self.container.inputData.ALIGNMENTIN.__str__(),self.startFileName)
            else:
                sequenceStack = []
                sequenceStack_check = []
                for sequenceFile in self.container.inputData.SEQIN:
                    possibleNameRoot = '>'+sequenceFile.fileContent.identifier.__str__()
                    if len(possibleNameRoot.strip()) == 1:
                        possibleNameRoot = ">unk"
                    possibleName = possibleNameRoot
                    iCount = 0
                    while possibleName in sequenceStack_check:
                        possibleName = possibleNameRoot + '_' + str(iCount)
                        iCount += 1
                    sequenceStack.append([possibleName, sequenceFile.fileContent.sequence.__str__()])
                    sequenceStack_check.append(possibleName)
                if sys.version_info > (3,0):
                    with open(self.startFileName,'w') as startFile:
                        for identifier,seq in sequenceStack:
                            startFile.write(identifier+'\n')
                            startFile.write(seq+'\n')
                else:
                    with open(self.startFileName,'wb') as startFile:
                        for identifier,seq in sequenceStack:
                            startFile.write(identifier+'\n')
                            startFile.write(seq+'\n')
        else:
            self.appendErrorReport(200)
            return CPluginScript.FAILED
        return CPluginScript.SUCCEEDED

    def processOutputFiles(self):
        try:
            with open(self.dndFilepath,'r') as aFile:
                aNode = etree.SubElement(self.xmlroot, 'Dendogram')
                aNode.text = aFile.read()
            with open(self.statsFilepath,'r') as aFile:
                aNode = etree.SubElement(self.xmlroot, 'Statistics')
                aNode.text = aFile.read()
            with open(self.container.outputData.ALIGNMENTOUT.__str__(),'r') as aFile:
                aNode = etree.SubElement(self.xmlroot, 'Alignment')
                aNode.text = aFile.read()
            with open (self.makeFileName('LOG'),'r') as logFile:
                lines = logFile.readlines()
                bestPair = None
                bestScore = -99999
                for line in lines:
                    if line.startswith('Sequences ('):
                        tokens = line.split()
                        pairwiseScoreNode = etree.SubElement(self.xmlroot,'PairwiseScore')
                        scoreNode = etree.SubElement(pairwiseScoreNode,'Score')
                        scoreNode.text = tokens[-1]
                        for iPartner in range(2):
                            partnerNode = etree.SubElement(pairwiseScoreNode,'Partner')
                            partnerNode.text = tokens[1][1:-1].split(':')[iPartner]
                        if int(scoreNode.text)>bestScore:
                            bestScore = int(scoreNode.text)
                            bestPair = tokens[1][1:-1].split(':')
                if bestPair is not None:
                    bestPairNode = etree.SubElement(self.xmlroot,'BestPair')
                    scoreNode = etree.SubElement(bestPairNode,'Score')
                    scoreNode.text = str(bestScore)
                    for iPartner in range(2):
                        partnerNode = etree.SubElement(bestPairNode,'Partner')
                        partnerNode.text = bestPair[iPartner]

        
            with open (self.makeFileName('PROGRAMXML'),'w') as programXML:
                CCP4Utils.writeXML(programXML,etree.tostring(self.xmlroot,pretty_print=True))
        except:
            return CPluginScript.FAILED

        anno = 'Alignment: '
        if self.container.inputData.SEQUENCELISTORALIGNMENT == 'ALIGNMENT':
          anno = anno + re.sub('Alignment: ','',str(self.container.inputData.ALIGNMENTIN.annotation)) + ', '
        else:
          for seq in self.container.inputData.SEQIN:
            anno = anno + str(seq.annotation) + ', '
        self.container.outputData.ALIGNMENTOUT.annotation = anno[0:-2]

        return CPluginScript.SUCCEEDED



    
