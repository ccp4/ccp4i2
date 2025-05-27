import io

from Bio import AlignIO, SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Blast import NCBIStandalone
from Bio.Seq import IUPACData as IUPAC
from lxml import etree

from ....core import CCP4Utils
from ....core.CCP4PluginScript import CPluginScript


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
        status = CPluginScript.SUCCEEDED
        mode = str(self.container.controlParameters.PASTEORREAD)

        if mode in ['ALIGNIN','PASTE']:
          alignment, format, commentary = importAlignment(self.container.inputData.ALIGNIN.__str__())
          commentary = commentary.getvalue()
          if alignment is None:
            self.appendErrorReport(202,self.container.inputData.ALIGNIN.__str__())
            status =  CPluginScript.FAILED
          else:
            outputString = io.StringIO()
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
        root = ET.Element('ProvideAlignment')
        commentaryNode = ET.SubElement(root,"Commentary")
        commentaryNode.text = commentary
        formatNode = ET.SubElement(root,'Format')
        formatNode.text = format
        if status == CPluginScript.SUCCEEDED:
          alignmentNode = ET.SubElement(root,'Alignment')
          alignmentNode.text = etree.CDATA(alignmentText)
          # Try putting sequences in to xml
          try:
            fileType, fileContent = self.container.outputData.ALIGNMENTFILE.identifyFile()
            for constituentSeq in fileContent:
              newSequence = ET.SubElement(root, 'Sequence')
              seqName = ET.SubElement(newSequence,'Identifier')
              seqName.text = str(constituentSeq)
          except:
            pass
            
        with open (self.makeFileName('PROGRAMXML'),'w') as programXML:
            CCP4Utils.writeXML(programXML,etree.tostring(root,pretty_print=True))
       
        return status

def importAlignment(filePath):
    #This tries to interpret the contents of a file in terms of a multiple sequence alignment
    #by using Biopython to read it in turn as one of a set of standard biopython sequence alignment
    #formats, and then as a bunch of catenated sequences
    #In the latter case, it does some checking (namely that the ltters are all in the IOPAC alphabet)
    #and some clean up (namely padding all of the sequences to maximum sequence length)

    outputString = io.StringIO()

    try:
        format = 'unknown'
        alignment = None
        for trialFormat in ['clustal','pir','fasta','stockholm','phylip']:
            outputString.write("Trying alignment format: "+trialFormat)
            try:
                with open(filePath,"r") as aliFile:
                    alignment = AlignIO.read(aliFile, trialFormat)
                    for iseq, seq in enumerate(alignment):
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
