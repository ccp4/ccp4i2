import os
import re
import pickle

from lxml import etree
import phaser_ext

from ......core import CCP4Data
from ......core.CCP4PluginScript import CPluginScript
from ......pimple import MGQTmatplotlib


class CallbackObject(object):
    def __init__(self, xmlroot=None, xmlResponders = []):
        super(CallbackObject,self).__init__()
        self.xmlroot = xmlroot
        self.xmlResponders = xmlResponders
    def currentActivityNode(self):
        try:
            aNode = self.xmlroot.xpath('CurrentActivity')[0]
        except:
            aNode = etree.SubElement(self.xmlroot,'CurrentActivity')
            self.currentActivityLabel = etree.SubElement(aNode,'label')
            self.currentActivityMax = etree.SubElement(aNode,'max')
            self.currentActivityValue = etree.SubElement(aNode,'value')
        return aNode
    def startProgressBar(self, arg1, arg2):
        aNode = self.currentActivityNode()
        self.currentActivityLabel.text = str(arg1)
        self.currentActivityMax.text = str(arg2)
        self.currentActivityValue.text = '0'
        self.notifyResponders()
    def incrementProgressBar(self):
        #print '**Increment progress bar called\n'
        aNode = self.currentActivityNode()
        self.currentActivityValue.text = str(int(self.currentActivityValue.text) + 1)
        self.notifyResponders()
    def endProgressBar(self):
        #print '**End progress bar called\n'
        pass
    def warn(self, warning):
        try:
            warningsElement = self.xmlroot.xpath('PhaserWarnings')[0]
        except: warningsElement = etree.SubElement(self.xmlroot,'PhaserWarnings')
        warningElement = etree.SubElement(warningsElement,'Warning')
        warningElement.text = warning
        self.notifyResponders()
    def loggraph(self, arg1, arg2):
        #print '\n**loggraph called',arg1
        try:
            tableelement = MGQTmatplotlib.CCP4LogToEtree(arg2)
            self.xmlroot.append(tableelement)
        except:
            print('\n\n\ failed reading MGQTMatPlotLib')
        self.notifyResponders()
    def call_back(self, arg1, arg2):
        #print '\n**call_back called: *Arg1:', arg1, '*Arg2:', arg2
        self.notifyResponders()
    def notifyResponders(self):
        try:
            for boundMethod in self.xmlResponders:
                boundMethod(self.xmlroot)
        except:
            print('Exception raised in notifyResponders')
    def addResponder(self, responder):
        self.xmlResponders.append(responder)

class phaser_MR(CPluginScript):

    TASKNAME = 'phaser_MR'                                  # Task name - should be same as class name
    TASKCOMMAND = ''                                     # The command to run the executable
    TASKVERSION= 0.0                                     # Version of this plugin
    COMTEMPLATE = None                                   # The program com file template
    COMTEMPLATEFILE = None                               # Name of file containing com file template

    ERROR_CODES = {101 : { 'description' : 'Failed to run MR_DAT' },
                   102 : { 'description' : 'Failed to set keywords' },
                   103 : { 'description' : 'Failed to parse Ensembles' },
                   104 : { 'description' : 'Failed to parse Asymmetric unit contents' },
                   105 : { 'description' : 'Phaser error' },
                   106 : { 'description' : 'Failed to parse Solutions' },
                   107 : { 'description' : 'Failed to parse Solutions' },
                   108 : { 'description' : 'Failed to add Fixed solutions' },}

    '''
    def __init__(self,parent=None,name=None,workDirectory=''):
      CPluginScript. __init__(self,parent=parent,name=name)
    '''
    excludedSections = ['macano', 'macmr', 'sgalternative', 'OUTP_LEVE']
    excludedKeywords = ['TITL', 'ROOT', 'HKLO', 'OUTP_LEVE', 'RESO_HIGH', 'RESO_LOW']
    requiredDefaultList = []
    
    def setKeywords(self, inputObject):
        #methodsOfInputObject = [getattr(inputObject, maybeMethod) for maybeMethod in dir(inputObject) if callable(getattr(inputObject, maybeMethod))]
        for parameterName, parameterObject in [(name, getattr(self.container.keywords, name)) for name in self.container.keywords.dataOrder() if name not in self.excludedSections]:
            #print parameterName, parameterObject
            if len(parameterObject.dataOrder()) > 0:
                for parameterName, parameterObject in [(name, getattr(parameterObject, name)) for name in parameterObject.dataOrder() if name not in self.excludedKeywords]:
                    self.setKeyword(inputObject, parameterName, parameterObject)
            else:
                self.setKeyword(inputObject, parameterName, parameterObject)
        
        rootPath = os.path.join(self.getWorkDirectory(), "PHASER")
        inputObject.setROOT(rootPath)
        return CPluginScript.SUCCEEDED
        
    def setKeyword(self, inputObject, parameterName, parameterObject):
        try:
            if parameterObject.isSet():
                setterName = 'set'+parameterName
                if hasattr(inputObject, setterName) and callable(getattr(inputObject,setterName)):
                    # don't set keywords set to their default values or 'Auto'
                    if (not parameterObject.isDefault() and not str(parameterObject) == 'Auto') or parameterName in self.requiredDefaultList:
                        if type(parameterObject) == CCP4Data.CString:
                            # the following functions expect the argument to be a tuple of floats that is converted to scitbx::vec3<floatType>
                            if parameterName in ['ROTA_EULE', 'TNCS_TRAN_VECT', 'TNCS_ROTA_ANGL', 'TRAN_STAR', 'TRAN_END']:
                                vec3 = tuple((float(i) for i in re.split(r',\s*|\s*',str(parameterObject).strip())))
                                setterEval = 'inputObject.'+setterName+'({})'.format(vec3)
                            else:
                                setterEval = 'inputObject.'+setterName+'("'+str(parameterObject)+'")'
                        else:
                            setterEval = 'inputObject.'+setterName+'('+str(parameterObject)+')'
                        print('Setting: ', setterEval)
                        eval(setterEval,{"__builtins__":None},{"inputObject":inputObject,"True":True,"False":False})
                    else:
                        print('parameterObject for {} isDefault:{} isAuto:{} in requiredDefaultList:{}'.format(setterName, parameterObject.isDefault(), str(parameterObject) == 'Auto', parameterName in self.requiredDefaultList))
                else:
                    print('Setter {} not present or not callable'.format(setterName))
        except:
            self.appendErrorReport(102)
            return CPluginScript.FAILED

    def parseContent(self, inputObject):
        inputData = self.container.inputData
        try:
            #ASymmetric unit contents
            if inputData.COMP_BY == 'DEFAULT':
                #Default is 50% solvent ?
                inputObject.setCOMP_BY("AVERAGE")
            elif inputData.COMP_BY == 'MW':
                if inputData.ASU_PROTEIN_MW.isSet():
                    inputObject.addCOMP_PROT_MW_NUM(float(inputData.ASU_PROTEIN_MW), 1.0)
                if inputData.ASU_NUCLEICACID_MW.isSet():
                    inputObject.addCOMP_NUCL_MW_NUM(float(inputData.ASU_NUCLEICACID_MW), 1.0)
            elif inputData.COMP_BY == 'ASU':
                inputObject.setCOMP_BY("ASU")
                indx = -1
                for seqObj in inputData.ASUFILE.fileContent.seqList:
                    indx += 1
                    if inputData.ASUFILE.selection[str(seqObj.name)]:
                        seqFilename = os.path.join(self.workDirectory,'ASUCONTENT_'+str(indx+1)+'.fasta')
                        inputData.ASUFILE.writeFasta(seqFilename,indx)
                        if seqObj.polymerType in ['RNA', 'DNA']:
                            inputObject.addCOMP_NUCL_SEQ_NUM(seqFilename,int(seqObj.nCopies))
                        else:
                            inputObject.addCOMP_PROT_SEQ_NUM(seqFilename,int(seqObj.nCopies))
        except Exception as e:
            self.appendErrorReport(104,details=str(e))
            return CPluginScript.FAILED
        return CPluginScript.SUCCEEDED

    def parseEnsembles(self, inputObject):
        inputData = self.container.inputData
        if True:
            for i, ensemble in enumerate(inputData.ENSEMBLES):
                for j, pdbItem in enumerate(ensemble.pdbItemList):
                    #Get selected subset of atoms if appropriate
                    atomsFile = str(pdbItem.structure)
                    #MN added following to force the selection of coordinates to make
                    #sure that the file used is a PDB (as output from getSelectedAtomsPdbFile(atomsFile)
                    #HJ reverted 25/04/19 as it does not work with multi-model PDB files
                    #if not pdbItem.structure.isSelectionSet():
                    #   pdbItem.structure.selection.text="*"
                    #End MN edit
                    if pdbItem.structure.isSelectionSet():
                        try:
                            atomsFile = os.path.join(self.workDirectory,'selectedAtomModel_'+str(i)+'_'+str(j)+'.pdb')
                            pdbItem.structure.getSelectedAtomsPdbFile(atomsFile)
                        except:
                            self.appendErrorReport(107,'Issue with getSelected')
                            return CPluginScript.FAILED
                    # CARD ON is selected by unSet pdbItem.identity_to_target and pdbItem.rms_to_target  
                    if not pdbItem.identity_to_target.isSet() and not pdbItem.rms_to_target.isSet():
                        # Check for special remark
                        cardFound = False
                        with open (atomsFile,'r') as pdbFile:
                            for line in pdbFile:
                                if "REMARK PHASER ENSEMBLE MODEL" in line:
                                    cardFound = True
                                    break
                                elif line[0:4] == "ATOM":
                                    break
                        if cardFound:
                            self.appendErrorReport(107,'So far...')
                            inputObject.addENSE_PDB_CARD(str(ensemble.label), atomsFile, True)
                            self.appendErrorReport(107,[a for a in dir(inputObject) if 'addENSE' in a].__str__())
                        else:
                            self.appendErrorReport(107,'Failed to addENSE_PDB_CARD')
                    else:
                        if pdbItem.identity_to_target.isSet() and pdbItem.identity_to_target is not None:
                            try:
                                inputObject.addENSE_PDB_ID(str(ensemble.label), atomsFile, float(pdbItem.identity_to_target))
                            except:
                                self.appendErrorReport(107,'Failed to addENSE_PDB_ID')
                                return CPluginScript.FAILED
                        elif pdbItem.rms_to_target.isSet() and pdbItem.rms_to_target is not None:
                            try:
                                inputObject.addENSE_PDB_RMS(str(ensemble.label), atomsFile, float(pdbItem.rms_to_target))
                            except:
                                self.appendErrorReport(107,'Failed to addENSE_PDB_RMS')
                                return CPluginScript.FAILED
                    
        else:
            self.appendErrorReport(103,'Failed to parse ensembles')
            return CPluginScript.FAILED
        return CPluginScript.SUCCEEDED

    def parseSolutions(self, inputObject):
        inp = self.container.inputData
        if inp.MODE_TY == "MR_FTF":
            try: 
                if inp.RFILEIN.isSet() and len(inp.USINGSOLELEMENTS) > 0:
                    with open(str(self.container.inputData.RFILEIN.fullPath),'r') as rf_file:
                        loadrf = pickle.load(rf_file)
                        filteredSolutionsObject = phaser_ext.mr_solution()
                        for solutionString in inp.USINGSOLELEMENTS:
                            iSol = int(solutionString.split(":")[0].strip())
                            trial = loadrf.__getitem__(0).RLIST[iSol]
                            inputObject.addSOLU_TRIAL_ENSE_EULER(trial.MODLID, trial.EULER, trial.RF, trial.RFZ)
            except:
                self.appendErrorReport(106)
                return CPluginScript.FAILED
        elif inp.MODE_TY in ['MR_AUTO', 'MR_PAK', 'MR_RNP']:
            try:
                #Solutions for elements of the Ensemblelist provided by the user
                if inp.SOLIN.isSet() and len(inp.USINGSOLELEMENTS) > 0:
                    with open(str(self.container.inputData.SOLIN.fullPath),'r') as file:
                        solutionsObject = pickle.load(file)
                        filteredSolutionsObject = phaser_ext.mr_solution()
                        for solutionString in inp.USINGSOLELEMENTS:
                            iSol = int(solutionString.split(":")[0].strip())
                            filteredSolutionsObject.append(solutionsObject.__getitem__(iSol))
                        inputObject.setSOLU(filteredSolutionsObject)
            except:
                self.appendErrorReport(106)
                return CPluginScript.FAILED
        try:
            if self.container.inputData.FIXENSEMBLES.isSet():
                for ensembleName in self.container.inputData.FIXENSEMBLES:
                    inputObject.addSOLU_ORIG_ENSE(ensembleName.__str__())
        except:
            self.appendErrorReport(108)
            return CPluginScript.FAILED
        return CPluginScript.SUCCEEDED

def subElementWithText(parent, type, value):
    result = etree.SubElement(parent,type)
    result.text = str(value)
    return result


def xmlFromSol(solText, parentNode):
    #print 'Soltext is ',solText
    solutionElement = None
    for line in solText.split('\n'):
        if line.strip().startswith( 'SOLU SET' ):
            solutionElement = etree.SubElement(parentNode,'Solution')
            splitLine = line.strip().split()
            if len(splitLine) > 2:
                annotationElement = etree.SubElement(solutionElement,'Annotation')
                annotationElement.text = ' '.join(splitLine[2:])
        elif line.strip().startswith( 'SOLU SPAC' ):
            spaceGroupElement = etree.SubElement(solutionElement,'spaceGroup')
            spaceGroupElement.text = line.strip()[10:]
        elif line.strip().startswith('SOLU 6DIM'):
            componentElement=etree.SubElement(solutionElement,'Component')
            nameElement = etree.SubElement(componentElement,'Name')
            nameElement.text = line.strip().split()[3]
        elif not line.strip().startswith('SOLU') and not line.strip().startswith('#') and solutionElement is not None:
            #Here if the annotation is spread over multiple lines:
            annotationElements = solutionElement.xpath('Annotation')
            if len(annotationElements) > 0: annotationElements[0].text += (' ' + line.strip())
    solutionSets = parentNode.xpath('Solution')
    for solutionSet in solutionSets:
        expandSolutionAnnotation(solutionSet)

def parseAnnotation(annotation):
    noNewLines = " ".join(annotation.split('\n'))
    words = [word.strip() for word in noNewLines.split()]
    
    solutions = []
    inAmalgamation = False
    for word in words:
        if word.startswith('RFZ='):
            solutions.append({'RFZ':word.split('=')[1],'TFZ':'-','PAK':'-','LLG':'-','refTFZ':'-','overallLLG':'-'})
        elif word.startswith('RF++'):
            solutions.append({'RFZ':'=','TFZ':'-','PAK':'-','LLG':'-','refTFZ':'-','overallLLG':'-'})
        elif word.startswith('RF*0'):
            solutions.append({'RFZ':'Zero Rot.','TFZ':'-','PAK':'-','LLG':'-','refTFZ':'-','overallLLG':'-'})
        elif word.startswith('(&'):
            inAmalgamation = True
        elif word.startswith('TFZ'):
            if inAmalgamation:
                solutions.append({'RFZ':'Amalg.','TFZ':'-','PAK':'-','LLG':'-','refTFZ':'-','overallLLG':'-'})
            if word.startswith('TFZ*0'):
                solutions[-1]['TFZ'] = 'Zero trans.'
            elif word.startswith('TFZ=*'):
                solutions[-1]['TFZ'] = 'P1 - Arbitrary'
            elif word.startswith('TFZ=='):
                solutions[-1]['refTFZ'] = word.split('==')[1]
            elif word.startswith('TFZ='):
                solutions[-1]['TFZ'] = word.split('=')[1]
            if solutions[-1]['TFZ'].endswith(')'):
                solutions[-1]['TFZ'] = solutions[-1]['TFZ'][:-1]
                inAmalgamation = False
        elif word.startswith('PAK='):
            try:
                if len(solutions) == 0:
                    solutions.append({'RFZ':'-','TFZ':'-','PAK':word.split('=')[1],'LLG':'-','refTFZ':'-','overallLLG':'-'})
                else:
                    solutions[-1]['PAK'] = word.split('=')[1]
                    if solutions[-1]['PAK'].endswith(')'):
                        solutions[-1]['PAK'] = solutions[-1]['PAK'][:-1]
                        inAmalgamation = False
            except Exception as inst:
                    print(inst)
        elif word.startswith('LLG='):
            #For phaser_rnp, we might get here with no solution...check and make one if necessary
            if len(solutions)<1:
                solutions.append({'RFZ':'Amalg.','TFZ':'-','PAK':'-','LLG':'-','refTFZ':'-','overallLLG':'-'})
            if 'LLG' in solutions[-1] and solutions[-1]['LLG'] != '-':
                solutions[-1]['overallLLG'] = word.split('=')[1]
            else:
                solutions[-1]['LLG'] = word.split('=')[1]
                if solutions[-1]['LLG'].endswith(')'):
                    solutions[-1]['LLG'] = solutions[-1]['LLG'][:-1]
                    inAmalgamation = False
    return solutions

def expandSolutionAnnotation(solutionSet):
    componentNodes = solutionSet.xpath('Component')
    annotationNodes = solutionSet.xpath('Annotation')
    if len(annotationNodes)>0:
        annotation = annotationNodes[0].text
        solnArray = parseAnnotation(annotation)
        for iComponent, solution in enumerate(solnArray):
            try:
                for property in solution:
                    result = subElementWithText(componentNodes[iComponent],property,solution[property])
            except:
                pass
