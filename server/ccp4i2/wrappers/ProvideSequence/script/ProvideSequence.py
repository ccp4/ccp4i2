from __future__ import print_function


from ccp4i2.core.CCP4PluginScript import CPluginScript
import os
import sys

try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO
from ccp4i2.core import CCP4ModelData
from ccp4i2.core import CCP4Utils
from lxml import etree

class ProvideSequence(CPluginScript):

    TASKNAME = 'ProvideSequence'                                  # Task name - should be same as class name
    TASKCOMMAND = ''                                     # The command to run the executable
    TASKVERSION= 0.0                                     # Version of this plugin
    COMTEMPLATE = None                                   # The program com file template
    COMTEMPLATEFILE = None                               # Name of file containing com file template
    RUNEXTERNALPROCESS=False

    '''
    def __init__(self,parent=None,name=None,workDirectory=''):
      CPluginScript. __init__(self,parent=parent,name=name)
    '''
    
    def startProcess(self, command, **kw):
        import tempfile
        from ccp4i2.wrappers.ProvideAlignment.script.ProvideAlignment import importAlignment
        
        root = etree.Element('ProvideSequence')
        
        # Create a temporary file to store the sequence(s) that will be used
        tempFile = tempfile.NamedTemporaryFile(suffix='.txt',delete=False)
        if sys.version_info > (3,0):
            tempFile.file.write(self.container.controlParameters.SEQUENCETEXT.__str__().encode('utf-8'))
        else:
            tempFile.file.write(self.container.controlParameters.SEQUENCETEXT.__str__())
        tempFile.close()
        
        #Attempt to interpret that as an alignment and/or stack of sequences
        alignment, format, commentary = importAlignment(tempFile.name)
        
        commentaryNode = etree.SubElement(root,"Commentary")
        commentaryNode.text = commentary.getvalue()
        
        if alignment is None:
            with open(self.makeFileName('PROGRAMXML'),'w') as programXML:
                CCP4Utils.writeXML(programXML,etree.tostring(root, pretty_print=True))
            self.reportStatus(CPluginScript.UNSATISFACTORY)
            return
        
        formatNode = etree.SubElement(root,'Format')
        formatNode.text = format
        
        from Bio import SeqIO
        for iSeq, seq in enumerate(alignment):
            outputList = self.container.outputData.SEQUENCEFILE_LIST
            outputList.append(outputList.makeItem())
            outputFile = outputList[-1]
            outputFile.setFullPath(os.path.normpath(os.path.join(self.getWorkDirectory(),'SEQUENCE'+str(iSeq)+'.fasta')))
            outputString = StringIO()
            with open(outputFile.__str__(),'w') as outputFileHandle:
                SeqIO.write([seq],outputFileHandle,'fasta')
            outputFile.annotation = seq.id + '-' + seq.description
        
            sequenceElement = etree.SubElement(root,'Sequence')
            outputString = StringIO()
            SeqIO.write([seq],outputString,'fasta')
            sequenceElement.text = outputString.getvalue()
            for property in ['id','name','description','seq']:
                newElement = etree.SubElement(sequenceElement,property)
                newElement.text = str(getattr(seq,property,'Undefined'))

        with open(self.makeFileName('PROGRAMXML'),'w') as programXML:
            CCP4Utils.writeXML(programXML,etree.tostring(root, pretty_print=True))
        return CPluginScript.SUCCEEDED

    def processOutputFiles(self, *args, **kwargs):
        if len(self.container.outputData.SEQUENCEFILE_LIST) == 0:
            return CPluginScript.SUCCEEDED
        try:
            self.container.outputData.CASUCONTENTOUT.loadFile()
            for iFile, sequenceFile in enumerate(self.container.outputData.SEQUENCEFILE_LIST):
                sequenceFile.loadFile()
                if iFile > 0:
                    self.container.outputData.CASUCONTENTOUT.fileContent.seqList.append(CCP4ModelData.CAsuContentSeq())
                entry = self.container.outputData.CASUCONTENTOUT.fileContent.seqList[-1]
                entry.nCopies.set(1)
                entry.sequence.set(sequenceFile.fileContent.sequence)
                entry.name.set(sequenceFile.fileContent.identifier)
                entry.description.set(sequenceFile.fileContent.description)
                entry.autoSetPolymerType()
            self.container.outputData.CASUCONTENTOUT.saveFile()
        except Exception as err:
            print("Failed to create CASUCONTENTOUT with error", err)
        return CPluginScript.SUCCEEDED
