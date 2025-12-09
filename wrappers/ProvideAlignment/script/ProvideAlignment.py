from __future__ import print_function

try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO

from ccp4i2.core.CCP4PluginScript import CPluginScript
from ccp4i2.core import CCP4Utils
import os


class ProvideAlignment(CPluginScript):

    TASKNAME = 'ProvideAlignment'                                  # Task name - should be same as class name
    TASKCOMMAND = ''                                     # The command to run the executable
    TASKVERSION= 0.0                                     # Version of this plugin
    COMTEMPLATE = None                                   # The program com file template
    COMTEMPLATEFILE = None                               # Name of file containing com file template
    RUNEXTERNALPROCESS=False

    ERROR_CODES = { 201 : { 'description' : 'Failed writing standard Clustal alignment file - did input sequence lengths match?' },
                    202 : { 'description' : 'Failed reading input alignment file' },
                    203 : { 'description' : 'Failed writing alignment file extracted from HHPred file' }
                    }
    
    '''
    def __init__(self,parent=None,name=None,workDirectory=''):
      CPluginScript. __init__(self,parent=parent,name=name)
    '''
    
    def startProcess(self, command, **kw):
        import os
        from lxml import etree        
        from Bio import AlignIO, SeqIO

        status = CPluginScript.SUCCEEDED
        mode = str(self.container.controlParameters.PASTEORREAD)

        if mode in ['ALIGNIN','PASTE']:
          alignment, format, commentary = importAlignment(self.container.inputData.ALIGNIN.__str__())
          commentary = commentary.getvalue()
          if alignment is None:
            self.appendErrorReport(202,self.container.inputData.ALIGNIN.__str__())
            status =  CPluginScript.FAILED
          else:
            outputString = StringIO()
            try:
              with open(str(self.container.outputData.ALIGNMENTFILE),"w") as outputFile:
                AlignIO.write(alignment, outputFile ,'clustal')
              AlignIO.write(alignment, outputString ,'clustal')
              self.container.outputData.ALIGNMENTFILE.annotation.set(self.container.controlParameters.ANNOTATION.__str__())
            except:
              self.appendErrorReport(201,fileName)
              status = CPluginScript.FAILED
            else:
              alignmentText = outputString.getvalue()
              
        elif mode in  ['HHPREDIN','BLASTIN']:
          try:
            if mode == 'HHPREDIN':
              format = 'HHpred'
              alignmentText = self.container.inputData.HHPREDIN.fileContent.getAlignmentText(self.container.inputData.ALI_INDEX)
            else:
              format = 'Blast'
              alignmentText = self.container.inputData.BLASTIN.fileContent.getAlignmentText(self.container.inputData.ALI_INDEX)
            print('alignmentText',alignmentText)
          except CException as e:
            commentary = 'Failed extracting alignment from '+format+' file'
            self.extendErrorReport(e)
            status = CPluginScript.FAILED
          else:
            try:
              CCP4Utils.saveFile(str(self.container.outputData.ALIGNMENTFILE),alignmentText)
            except:
              self.appendErrorReport(203,str(self.container.outputData.ALIGNMENTFILE))
              status = CPluginScript.FAILED
              commentary = 'Failed writing alignment file extracted from '+format+' file'
            else:
              self.container.outputData.ALIGNMENTFILE.annotation.set(self.container.controlParameters.ANNOTATION.__str__())
              commentary = 'Succeeded extracting alignment from '+format+' file'
          
        # Write xml
        root = etree.Element('ProvideAlignment')
        commentaryNode = etree.SubElement(root,"Commentary")
        commentaryNode.text = commentary
        formatNode = etree.SubElement(root,'Format')
        formatNode.text = format
        if status == CPluginScript.SUCCEEDED:
          alignmentNode = etree.SubElement(root,'Alignment')
          alignmentNode.text = etree.CDATA(alignmentText)
          # Try putting sequences in to xml
          try:
            fileType, fileContent = self.container.outputData.ALIGNMENTFILE.identifyFile()
            for constituentSeq in fileContent:
              newSequence = etree.SubElement(root, 'Sequence')
              seqName = etree.SubElement(newSequence,'Identifier')
              seqName.text = str(constituentSeq)
          except:
            pass
            
        with open (self.makeFileName('PROGRAMXML'),'w') as programXML:
            CCP4Utils.writeXML(programXML,etree.tostring(root,pretty_print=True))
       
        return status

def getBioPyVersion():
    from Bio import __version__ as bioversion
    try:
        version=float(bioversion)
        return version
    except:
        return None

def importAlignment(filePath):
    #This tries to interpret the contents of a file in terms of a multiple sequence alignment
    #by using Biopython to read it in turn as one of a set of standard biopython sequence alignment
    #formats, and then as a bunch of catenated sequences
    #In the latter case, it does some checking (namely that the ltters are all in the IOPAC alphabet)
    #and some clean up (namely padding all of the sequences to maximum sequence length)
    from Bio import AlignIO, SeqIO
    
    outputString = StringIO()

    # Get the BioPython version number
    bioversion=getBioPyVersion()

    try:
        format = 'unknown'
        alignment = None
        for trialFormat in ['clustal','pir','fasta','stockholm','phylip']:
            outputString.write("Trying alignment format: "+trialFormat)
            try:
                with open(filePath,"r") as aliFile:
                    alignment = AlignIO.read(aliFile, trialFormat)
                    if bioversion is not None:
                        if bioversion < 1.79:
                            from Bio.Alphabet import IUPAC
                        else:
                            from Bio.Seq import IUPACData as IUPAC
                    else:
                        from Bio.Seq import IUPACData as IUPAC
                    #from Bio.Alphabet import IUPAC
                    for iseq, seq in enumerate(alignment):
                        if bioversion is not None:
                            if bioversion < 1.79:
                                letters = IUPAC.protein.letters+IUPAC.ambiguous_rna.letters+IUPAC.ambiguous_dna.letters+'-.X'
                            else:
                                letters = IUPAC.protein_letters+IUPAC.ambiguous_rna_letters+IUPAC.ambiguous_dna_letters+'-.X'
                        else:
                            letters = IUPAC.protein_letters+IUPAC.ambiguous_rna_letters+IUPAC.ambiguous_dna_letters+'-.X'
                        for iletter, letter in enumerate(seq):
                            if letter not in letters and letter not in letters.lower():
                                raise ValueError("Found unexpected character: '" + letter+"' sequence "+str(iseq)+" position "+str(iletter))
                    format= trialFormat
                outputString.write("...Success\n\n")
                break
            except Exception as e:
                alignment = None
                outputString.write("...Failed: "+str(e)+"\n\n")

        if alignment is None:
            from Bio.Align import MultipleSeqAlignment
            if bioversion is not None:
               if bioversion < 1.79:
                   from Bio.Alphabet import IUPAC, _verify_alphabet
               else:
                   from Bio.Seq import IUPACData as IUPAC
            else:
               from Bio.Seq import IUPACData as IUPAC, _verify_alphabet
            for seqFormat in ['fasta','pir','seqxml','embl','ncbi']:
                outputString.write("Trying Sequence format: "+seqFormat)
                try:
                    msa = []
                    with open(filePath,'r') as seqFile:
                        for seq in SeqIO.parse(seqFile, seqFormat):
                            msa.append(seq)
                    #Pad to the same length
                    longestSeq = max([len(seq.seq) for seq in msa])
                    paddedMSAs = []
                    for iseq, seq in enumerate(msa):
                        if bioversion is not None:
                            if bioversion < 1.79:
                                letters = IUPAC.protein.letters+IUPAC.ambiguous_rna.letters+IUPAC.ambiguous_dna.letters+'-.'
                            else:
                                letters = IUPAC.protein_letters+IUPAC.ambiguous_rna_letters+IUPAC.ambiguous_dna_letters+'-.'
                        else:
                            letters = IUPAC.protein_letters+IUPAC.ambiguous_rna_letters+IUPAC.ambiguous_dna_letters+'-.'
                        for iletter, letter in enumerate(seq):
                            if letter not in letters and letter not in letters.lower():
                                raise ValueError("Found unexpected character: '" + letter+"' sequence "+str(iseq)+" position "+str(iletter))
                        padder = ''
                        for ipad in range(longestSeq - len(seq.seq)): padder += '-'
                        new = seq+padder
                        paddedMSAs.append(new)
                
                    alignment = MultipleSeqAlignment(paddedMSAs)
                    format = seqFormat
                    outputString.write("...Success\n\n")
                    break
                except Exception as e:
                    alignment = None
                    format='unknown'
                    outputString.write("...Failed: "+str(e)+"\n\n")


        if alignment is None:
            #Last chance saloon
            with open(filePath,'r') as result_handle:
                from Bio.Blast import NCBIStandalone
                outputString.write("Trying to interpret as legacy (non XML) BLAST file")
                blast_parser = NCBIStandalone.BlastParser()
                try:
                    blast_record = blast_parser.parse(result_handle)
                    for alignment in blast_record.alignments:
                        for hsp in alignment.hsps:
                            if hsp.expect < E_VALUE_THRESH:
                                print('****Alignment****')
                                print(('sequence:', alignment.title))
                                print(('length:', alignment.length))
                                print(('e value:', hsp.expect))
                                print((hsp.query[0:75] + '...'))
                                print((hsp.match[0:75] + '...'))
                                print((hsp.sbjct[0:75] + '...'))
                except Exception as e:
                    outputString.write("...Failed: "+str(e)+"\n")
                
        return alignment, format, outputString
    except:
        return None, 'unknown', outputString




#====================================================================================================
# PLUGIN TESTS
# See Python documentation on unittest module

import unittest

class testProvideAlignment(unittest.TestCase):

   def setUp(self):
    # make all background jobs wait for completion
    # this is essential for unittest to work
    from ccp4i2.core.CCP4Modules import QTAPPLICATION,PROCESSMANAGER
    self.app = QTAPPLICATION()
    PROCESSMANAGER().setWaitForFinished(10000)

   def tearDown(self):
    from ccp4i2.core.CCP4Modules import PROCESSMANAGER
    PROCESSMANAGER().setWaitForFinished(-1)

   def test_1(self):
     from ccp4i2.core.CCP4Modules import QTAPPLICATION
     wrapper = ProvideAlignment(parent=QTAPPLICATION(),name='ProvideAlignment_test1')
     wrapper.container.loadDataFromXml()
     

def TESTSUITE():
  suite = unittest.TestLoader().loadTestsFromTestCase(testProvideAlignment)
  return suite

def testModule():
  suite = TESTSUITE()
  unittest.TextTestRunner(verbosity=2).run(suite)
