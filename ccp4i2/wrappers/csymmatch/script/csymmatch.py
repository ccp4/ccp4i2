import os
import pathlib

from lxml import etree

from ....core import CCP4Utils
from ....core.CCP4PluginScript import CPluginScript


class csymmatch(CPluginScript):

    TASKTITLE='Csymmatch - move model to a reference by using symmetry'
    TASKNAME = 'csymmatch'
    TASKMODULE = 'molecular_replacement'
    TASKCOMMAND = 'csymmatch'
    TASKVERSION = 0.0
    ASYNCHRONOUS = False
    MAINTAINER = 'liz.potterton@york.ac.uk'

    def makeCommandAndScript(self):
      inp = self.container.inputData
      par = self.container.controlParameters
      out = self.container.outputData

      if inp.XYZIN_QUERY.isSelectionSet():
        xyzin_query_file = os.path.join(self.getWorkDirectory(),'XYZIN_QUERY_sel.pdb')
        self.container.inputData.XYZIN_QUERY.loadFile()
        if self.container.inputData.XYZIN_QUERY.isMMCIF():
            xyzin_query_file = str(pathlib.Path(xyzin_query_file).with_suffix('.cif'))
            out.XYZOUT.setFullPath(str(pathlib.Path(out.XYZOUT.fullPath.__str__()).with_suffix('.cif')))
        inp.XYZIN_QUERY.getSelectedAtomsPdbFile(xyzin_query_file)
      else:
        xyzin_query_file = inp.XYZIN_QUERY.fullPath.__str__()
        self.container.inputData.XYZIN_QUERY.loadFile()
        if self.container.inputData.XYZIN_QUERY.isMMCIF():
            out.XYZOUT.setFullPath(str(pathlib.Path(out.XYZOUT.fullPath.__str__()).with_suffix('.cif')))
      if inp.XYZIN_TARGET.isSelectionSet():
        xyzin_target_file = os.path.join(self.getWorkDirectory(),'XYZIN_TARGET_sel.pdb')
        self.container.inputData.XYZIN_TARGET.loadFile()
        if self.container.inputData.XYZIN_TARGET.isMMCIF():
            xyzin_target_file = str(pathlib.Path(xyzin_target_file).with_suffix('.cif'))
        inp.XYZIN_TARGET.getSelectedAtomsPdbFile(xyzin_target_file)
      else:
        xyzin_target_file = inp.XYZIN_TARGET.fullPath.__str__()
      self.appendCommandLine( [ '-pdbin', xyzin_query_file ] )
      self.appendCommandLine( [ '-pdbin-ref', xyzin_target_file ] )
      if par.ORIGIN_HAND.isSet():
        if par.ORIGIN_HAND:
          self.appendCommandLine( [ '-origin-hand'] )
      if par.CONNECTIVITY_RADIUS.isSet():
        self.appendCommandLine( [ '-connectivity-radius',str(par.CONNECTIVITY_RADIUS)] )
      self.appendCommandLine( [ '-pdbout' , out.XYZOUT.fullPath.__str__() ] )

    def processOutputFiles(self):
        logName = self.makeFileName('LOG')
        
        xmlRoot = etree.Element('Csymmatch')
        segmentNode = None
        with open (logName,'r') as logFile:
            lines = logFile.readlines()
            for line in lines:
                if line.strip().startswith('Change of hand'):
                    handNode = etree.SubElement(xmlRoot,'ChangeOfHand')
                    handNode.text = line.strip().split(':')[1]
                elif line.strip().startswith('Change of origin'):
                    originNode = etree.SubElement(xmlRoot,'ChangeOfOrigin')
                    originNode.text = line.strip().split(':')[1]
                elif line.strip().startswith('Chain'):
                    segmentNode = etree.SubElement(xmlRoot,'Segment')
                    rangeNode = etree.SubElement(segmentNode,'Range')
                    rangeNode.text = line.strip().split('will')[0]
                elif line.strip().startswith('Symmetry operator'):
                    if segmentNode is not None:
                        operatorNode = etree.SubElement(segmentNode,'Operator')
                        operatorNode.text = line.strip().split(':')[1]
                elif line.strip().startswith('Lattice shift'):
                    if segmentNode is not None:
                        shiftNode = etree.SubElement(segmentNode,'Shift')
                        shiftNode.text = line.strip().split(':')[1]
                elif line.strip().startswith('with normalised score'):
                    if segmentNode is not None:
                        scoreNode = etree.SubElement(segmentNode,'Score')
                        scoreNode.text = line.strip().split(':')[1]
    
        with open(self.makeFileName('PROGRAMXML'),'w') as xmlFile:
            xmlString = etree.tostring(xmlRoot, pretty_print=True)
            CCP4Utils.writeXML(xmlFile,xmlString)


        return CPluginScript.SUCCEEDED
