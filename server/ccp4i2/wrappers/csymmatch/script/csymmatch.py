import pathlib

from ccp4i2.core import CCP4Utils
from ccp4i2.core.CCP4PluginScript import CPluginScript


class csymmatch(CPluginScript):

    TASKNAME = 'csymmatch'
    TASKCOMMAND = 'csymmatch'

    def makeCommandAndScript(self):

      inp = self.container.inputData
      par = self.container.controlParameters
      out = self.container.outputData

      import os
      # csymmatch reads both formats: keep whichever came in, and match the
      # output's name to the query's format.
      if inp.XYZIN_QUERY.isSelectionSet():
        xyzin_query_file = inp.XYZIN_QUERY.getSelectedAtomsFile(
            'XYZIN_QUERY_sel', self.getWorkDirectory())
      else:
        xyzin_query_file = inp.XYZIN_QUERY.fullPath.__str__()
      if inp.XYZIN_QUERY.isMMCIF():
        out.XYZOUT.setFullPath(str(pathlib.Path(out.XYZOUT.fullPath.__str__()).with_suffix('.cif')))

      if inp.XYZIN_TARGET.isSelectionSet():
        xyzin_target_file = inp.XYZIN_TARGET.getSelectedAtomsFile(
            'XYZIN_TARGET_sel', self.getWorkDirectory())
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

        from lxml import etree
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
