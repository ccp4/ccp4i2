from __future__ import print_function

import re

from ccp4i2.core.CCP4PluginScript import CPluginScript
from ccp4i2.core import CCP4Utils
import pathlib
import csv

class gesamt(CPluginScript):

    TASKTITLE='Gesamt - structural alignment'
    TASKNAME = 'gesamt'
    TASKMODULE= 'model_data_utility'
    TASKCOMMAND = 'gesamt'
    TASKVERSION= 0.0
    PERFORMANCECLASS = 'CSuperposePerformance'

    def makeCommandAndScript(self):

      inp = self.container.inputData
      par = self.container.controlParameters
      out = self.container.outputData

      import os

      self.csvout = os.path.join(self.getWorkDirectory(),"csv.out")

      if str(self.container.controlParameters.PAIRMULTI) == "MULTIPLE":
          for f in inp.XYZIN_LIST:
              xyzin_file = str( f.fullPath )
              f.loadFile()
              self.appendCommandLine( [ xyzin_file ] )

          """
          #I guess everythinbg moves.
          if par.OUTPUT_COORDS == "minusO":
              self.appendCommandLine( [ '-o', str( out.XYZOUT.fullPath ) ] )
          elif par.OUTPUT_COORDS == "minusOF":
              a,b = os.path.splitext ( inp.XYZIN_TARGET.fullPath )
              self.tmpOut = a + "_2" + b
              self.appendCommandLine( [ '-o-f' ] )
          """
          if par.MODE.isSet():
              self.appendCommandLine( [ '-'+str( par.MODE ) ] )

          self.appendCommandLine( [ '-csv',self.csvout ] )

          return CPluginScript.SUCCEEDED

      elif inp.XYZIN_TARGET.fullPath.isSet():
          if inp.XYZIN_TARGET.isSelectionSet():
              xyzin_target_file = os.path.join(self.getWorkDirectory(),"XYZIN_TARGET_sel.pdb")
              self.container.inputData.XYZIN_TARGET.getSelectedAtomsPdbFile(xyzin_target_file)
          else:
              xyzin_target_file = str( inp.XYZIN_TARGET.fullPath )
          self.appendCommandLine( [ xyzin_target_file ] )

      if inp.XYZIN_QUERY.isSelectionSet():
          xyzin_query_file = os.path.join(self.getWorkDirectory(),"XYZIN_QUERY_sel.pdb")
          self.container.inputData.XYZIN_QUERY.loadFile()
          if self.container.inputData.XYZIN_QUERY.isMMCIF():
            xyzin_query_file = str(pathlib.Path(xyzin_query_file).with_suffix('.cif'))
            out.XYZOUT.setFullPath(str(pathlib.Path(out.XYZOUT.fullPath.__str__()).with_suffix('.cif')))
          self.container.inputData.XYZIN_QUERY.getSelectedAtomsPdbFile(xyzin_query_file)
      else:
          xyzin_query_file = str( inp.XYZIN_QUERY.fullPath )
          self.container.inputData.XYZIN_QUERY.loadFile()
          if self.container.inputData.XYZIN_QUERY.isMMCIF():
            out.XYZOUT.setFullPath(str(pathlib.Path(out.XYZOUT.fullPath.__str__()).with_suffix('.cif')))

      self.appendCommandLine( [ xyzin_query_file ] )

      if par.MODE.isSet():
        self.appendCommandLine( [ '-'+str( par.MODE ) ] )

      if out.XYZOUT.fullPath.isSet():
        #SJM 13/12/2017 - This should always be minusO now after removing option from GUI.
        if par.OUTPUT_COORDS == "minusO":
          self.appendCommandLine( [ '-o', str( out.XYZOUT.fullPath ) ] )
        elif par.OUTPUT_COORDS == "minusOF":
          a,b = os.path.splitext ( inp.XYZIN_TARGET.fullPath )
          self.tmpOut = a + "_2" + b
          self.appendCommandLine( [ '-o-f' ] )

      self.appendCommandLine( [ '-csv',self.csvout ] )

      return CPluginScript.SUCCEEDED

    def processOutputFiles(self):

        aminoResidues = ["GLY","ALA","VAL","PRO","SER","THR","LEU","ILE","CYS","ASP","GLU","ASN","GLN","ARG","LYS","MET","MSE","HIS","PHE","TYR","TRP","HCS","ALO","PDD","UNK"]
        
        logName = self.makeFileName('LOG')
        eulerValues = None
        translationValues = None
        rmsValue = None
        qValue = None
        nResValue = None
        transformationMatrix = []
        
        from lxml import etree
        xmlRoot = etree.Element('Gesamt')
        perResidueNode = None
        iRes = 0

        import os
        out = self.container.outputData
        if hasattr(self,"tmpOut") and os.path.isfile(self.tmpOut):
          os.rename ( self.tmpOut,str(out.XYZOUT.fullPath) )

        with open (logName,'r') as logFile:
            lines = logFile.readlines()
            inPerResidue = False
            inTargetMatrix = False
            inMeatOfMatrix = False
            for line in lines:
                tokens = line.strip().split()
                if line.strip().startswith('Angle between rotation axis)'):
                    rotationAngle = [float(tokens[-1])]
                if line.strip().startswith('Polar angles (omega,phi,kappa)'):
                    polarValues = [float(tokens[-3]), float(tokens[-2]), float(tokens[-1])]
                if line.strip().startswith('Euler angles (alpha,beta,gamma)'):
                    eulerValues = [float(tokens[-3]), float(tokens[-2]), float(tokens[-1])]
                if line.strip().startswith('Orthogonal translation (Angstrom)'):
                    translationValues = [float(tokens[-3]), float(tokens[-2]), float(tokens[-1])]
                if line.strip().startswith('RMSD'):
                    rmsValue = float(tokens[-1])
                if line.strip().startswith('Q-score'):
                    qValue = float(tokens[-1])
                if line.strip().startswith('Sequence Id'):
                    seqId = float(tokens[-1])
                if line.strip().startswith('Aligned residues'):
                    nResValue = int(tokens[-1])
                if line.strip().startswith('.-------------.------------.-------------'):
                    inPerResidue = True
                    perResidueNode = etree.Element('PerResidue')
                if 'Transformation matrix for Target:' in line: inTargetMatrix = True
                if inTargetMatrix:
                    if inMeatOfMatrix:
                        if len(line.strip()) == 0:
                            inTargetMatrix = False
                            inMeatOfMatrix = False
                        else:
                            for value in line.strip().split():
                                transformationMatrix.append(float(value))
                    if 'Rx' in line:
                        inMeatOfMatrix = True
                if inPerResidue:
                    tokens = line.strip().split('|')
                    if len(tokens) > 1 and not ('Query' in tokens[1] or '-----' in tokens[1]):
                        row = etree.SubElement(perResidueNode,'Row')
                        calphaDistanceNode = etree.SubElement(row,'Distance')
                        calphaDistanceNode.text = tokens[2].strip()[3:-3]
                        equivalenceNode = etree.SubElement(row,'Equivalence')
                        equivalenceNode.text = line.strip()
                        iResNode = etree.SubElement(row,'iRes')
                        iResNode.text = str(iRes)
                        iRes += 1
                if line.strip().startswith("`-------------'------------'-------------'"):
                    inPerResidue = False

    
        if eulerValues is not None and translationValues is not None and rmsValue is not None and qValue is not None and nResValue is not None:
            from ccp4i2.core.CCP4MathsData import CTransformation
            self.container.outputData.TRANSFORMATION = CTransformation()
            self.container.outputData.TRANSFORMATION.alpha.set(eulerValues[0])
            self.container.outputData.TRANSFORMATION.beta.set(eulerValues[1])
            self.container.outputData.TRANSFORMATION.gamma.set(eulerValues[2])
            self.container.outputData.TRANSFORMATION.x.set(translationValues[0])
            self.container.outputData.TRANSFORMATION.y.set(translationValues[1])
            self.container.outputData.TRANSFORMATION.z.set(translationValues[2])
            self.container.outputData.PERFORMANCE.RMSxyz.set(rmsValue)
            #self.container.outputData.PERFORMANCE.QScore = qValue
            self.container.outputData.PERFORMANCE.nResidues.set(nResValue)
            transformationNode = etree.SubElement(xmlRoot,'Transformation',
                                                  omega=str(polarValues[0]),
                                                  phi=str(polarValues[1]),
                                                  kappa=str(polarValues[2]),
                                                  alpha=str(self.container.outputData.TRANSFORMATION.alpha),
                                                  beta=str(self.container.outputData.TRANSFORMATION.beta),
                                                  gamma=str(self.container.outputData.TRANSFORMATION.gamma),
                                                  x = str(self.container.outputData.TRANSFORMATION.x ),
                                                  y = str(self.container.outputData.TRANSFORMATION.y ),
                                                  z = str(self.container.outputData.TRANSFORMATION.z),
                                                  rms = str(rmsValue),
                                                  q = str( qValue ),
                                                  seqid = str(seqId),
                                                  nRes  = str(nResValue ) )
            matrixNode = etree.SubElement(transformationNode,'Matrix')
            
            #Storing the transform in 4x4 matrix format, we need to add the bottom row of a 4x4transformation matrix
            for element in transformationMatrix + [0.,0.,0.,1.]:
                elementNode = etree.SubElement(matrixNode,'Element')
                elementNode.text = str(element)

            if perResidueNode is not None: transformationNode.append(perResidueNode)
        
            if len(transformationMatrix) == 12:
                from ccp4mg import mmut
                from ccp4mg import pygl_coord
                inp = self.container.inputData
                xyzin_query_file = str( inp.XYZIN_QUERY.fullPath )
                molHnd = mmut.CMMANManager()
                molHnd.ReadCoorFile(xyzin_query_file)
                mat = pygl_coord.matrix(4,4,transformationMatrix + [0.,0.,0.,1.])
                molHnd.SetTransform(mat,False)
                molHnd.FinishStructEdit()
                molHnd.WritePDBASCII(str(out.XYZOUT_QUERY.fullPath))
                out.XYZOUT_QUERY.annotation.set("Transformed query model")
            else:
                print("Hmm, transformationMatrix is length",len(transformationMatrix))

        if os.path.exists(self.csvout):
            transformationCSVNode = etree.SubElement(xmlRoot,'TransformationCSV')
            perResidueNode = etree.SubElement(transformationCSVNode,'PerResidue')
            currentChain = "dummy_dummy_dummy"
            with open(self.csvout) as csvfile:
                spamreader = csv.reader(csvfile)
                inMultiDisplacement = False
                inPairwiseDistance = False
                inMultiTransformation = False
                inPairwiseTransformation = False
                inScore = False
                pairwiseTransformation = []
                multiTransformationsCombined = []
                multiTransformations = []
                
                for row in spamreader:
                    if (inMultiDisplacement or inPairwiseDistance) and len(row)>0 and len(row[0].strip())>0:
                        chain = row[-1].strip().split()[-2].split(":")[0]
                        if chain != currentChain:
                            chainNode = etree.SubElement(perResidueNode,'Chain')
                            currentChain = chain
                        name_num = (row[-1].split(":")[1])
                        resType = re.findall(r'[A-Z]+', name_num)[0]
                        iRes = re.findall(r'\d+', name_num)[0]
                        if aminoResidues.count(resType)>0:
                            rowNode = etree.SubElement(chainNode,'Row')
                            equivalenceNode = etree.SubElement(rowNode,"Equivalence")
                            distanceNode = etree.SubElement(rowNode,"Distance")
                            iResNode = etree.SubElement(rowNode,"iRes")
                            chainIDNode = etree.SubElement(rowNode,"chainID")
                            chainIDNode.text = chain
                            distanceNode.text = row[0].strip()
                            iResNode.text = iRes
                            equivalenceNode.text = ("|").join(row)

                    if inMultiTransformation:
                            if len(row)==4 and row[0].strip()!="Rx" and not row[0].strip().startswith("STRUCTURE"):
                                multiTransformationsCombined.extend([x.strip() for x in row])

                    if inPairwiseTransformation:
                            if len(row)==4 and row[0].strip()!="Rx" and not row[0].strip().startswith("Direction"):
                                pairwiseTransformation.extend([x.strip() for x in row])

                    if not inMultiDisplacement and len(row)>0 and row[0].startswith("Disp.[A]"):
                            inMultiDisplacement = True

                    if not inPairwiseDistance and len(row)>0 and row[0].startswith("Dist [A]"):
                            inPairwiseDistance = True

                    if not inMultiTransformation and len(row)>0 and row[0].strip().startswith("SUPERPOSITION MATRICES (ORTHOGONAL)"):
                            inMultiTransformation = True

                    if inMultiTransformation and len(row)>0 and row[0].strip().startswith("SCORES ACHIEVED"):
                            inMultiTransformation = False

                    if not inPairwiseTransformation and len(row)>0 and row[0].strip().startswith("Transformation matrix for Target"):
                            inPairwiseTransformation = True

                    if inPairwiseTransformation and len(row)>0 and row[0].strip().startswith("Direction cosines of the rotation axis"):
                            inPairwiseTransformation = False

                    if not inScore and len(row)>0 and (row[0].strip().startswith("SCORES ACHIEVED") or row[0].strip().startswith("SUPERPOSITION")):
                        inScore = True

                    if inScore and len(row)>0 and (row[0].strip().startswith("Transformation") or row[0].strip().startswith("Pairwise")):
                        inScore = False

                    if inScore and len(row)>0 and (row[0].strip().startswith("Q-score") or row[0].strip().startswith("quality Q")):
                        qNode = etree.SubElement(transformationCSVNode,'Q_CSV')
                        qNode.text = row[1]

                    if inScore and len(row)>0 and (row[0].strip().startswith("r.m.s.d") or row[0].strip().startswith("RMSD")):
                        rmsdNode = etree.SubElement(transformationCSVNode,'RMSD_CSV')
                        rmsdNode.text = row[1]
                        self.container.outputData.PERFORMANCE.RMSxyz.set(float(row[1]))

                    if inScore and len(row)>0 and (row[0].strip().startswith("Nalign") or row[0].strip().startswith("Aligned residues")):
                        nAlignNode = etree.SubElement(transformationCSVNode,'NALIGN_CSV')
                        nAlignNode.text = row[1]
                        self.container.outputData.PERFORMANCE.nResidues.set(int(row[1]))

                for i in range(len(multiTransformationsCombined)):
                    if i%12 == 0:
                        multiTransformations.append([])
                    multiTransformations[-1].append(multiTransformationsCombined[i])

                if len(pairwiseTransformation)>0:
                    matrixNode = etree.SubElement(transformationCSVNode,'matrixCSV')
                    matrixNode.text = ",".join(pairwiseTransformation)
                for transform in multiTransformations:
                    matrixNode = etree.SubElement(transformationCSVNode,'matrixCSV')
                    matrixNode.text = ",".join(transform)
                
        with open(self.makeFileName('PROGRAMXML'),'w') as xmlFile:
            xmlString = etree.tostring(xmlRoot, pretty_print=True)
            CCP4Utils.writeXML(xmlFile,xmlString)

        out.XYZOUT.annotation.set("Gesamt output file (selected atoms from query and target only)")
        return CPluginScript.SUCCEEDED
