from qtgui import CCP4TaskWidget

class CTaskmrparse(CCP4TaskWidget.CTaskWidget):

    TASKNAME = 'mrparse'
    TASKVERSION = 1.0
    TASKMODULE =['alpha_fold','bioinformatics']
    TASKTITLE = 'MrParse'
    SHORTTASKTITLE = "mrparse"
    DESCRIPTION = 'Search online PDB and EBI-AFDB databases to find and process search models for use in molecular replacement'
    WHATNEXT = []

    def __init__(self,parent):
        CCP4TaskWidget.CTaskWidget.__init__(self,parent)

    def drawContents(self):
        import multiprocessing
        MAXPROC=multiprocessing.cpu_count()  

        self.openFolder(folderFunction='inputData', title='Input Data')
        self.createLine(['subtitle', 'Select input data', 'Sequence is required. Observed intensities (or amplitudes) are optional but recommended'])
        self.openSubFrame(frame=[True])
        self.createLine(['tip', 'Input sequence', 'widget', 'SEQIN'])

        self.openSubFrame(frame=[True])
        self.createLine(['tip', 'Input reflections', 'widget', 'F_SIGF'])
        self.createLine(['advice',
                         'Observed intensities (or amplitudes) are optional but recommended for assessing crystal pathology and calculating eLLG scores'])
        self.closeSubFrame()

        self.openSubFrame(frame=[True])
        self.createLine( [ 'subtitle', 'Main Settings', 'Set max hits and seach database(s)' ] )
        self.createLine( [ 'advice', 'Maximum no. of search models to create for each database searched:' ] )
        self.createLine( [ 'widget', 'MAXHITS' ] )
        self.createLine( [ 'widget', 'DATABASE', 'label', 'Databases to search for models (PDB, AFDB, All)' ] )
        self.closeSubFrame()

#        self.createLine( [ 'subtitle', 'Model databases', 'Databases to search for possible search models' ] )
#        self.createLine( [ 'widget', 'SEARCH_PDB', 'label', 'Search PDB for for possible MR search models' ] )
#        self.createLine( [ 'advice', indent+'Non-redundancy level for homologue search:' ], toggle=['SEARCH_PDB', 'open', [ True ] ] )
#        self.createLine( [ 'label', indent, 'widget', 'REDUNDANCYLEVEL' ], toggle=['SEARCH_PDB', 'open', [ True ] ]  )

        self.openSubFrame(frame=[True], toggle = ['DATABASE', 'close', ['AFDB'] ] )
        self.createLine( [ 'subtitle', 'Optional PDB Settings', 'Set PDB database options' ] )
        self.createLine( [ 'widget', '-browseDb', True, 'PDBLOCAL', 'tip', 'Local PDB mirror' ], toggle=['DATABASE', 'close', ['AFDB'] ]  )# ,toggle=['PASTEORREAD','open',['HHPREDIN']])
        self.createLine( [ 'advice', 'Alternative pdb fasta sequence file (default uses internal CCP4 version)' ], toggle=['DATABASE', 'close', ['AFDB'] ]   )
        self.createLine( [ 'widget', '-browseDb', True, 'PDBSEQDB', 'tip', 'Alternative PDB sequences fasta file' ], toggle=['DATABASE', 'close', ['AFDB'] ]   )# ,toggle=['PASTEORREAD','open',['HHPREDIN']])
        self.closeSubFrame()

        self.openSubFrame(frame=[False], toggle = ['DATABASE', 'close', ['PDB'] ] )
        self.createLine( [ 'subtitle', 'Optional AFDB Settings', 'Set AFDB database options' ] )
        self.createLine( [ 'widget', 'USEAPI', 'label', 'Search EBI Alphafold database using the EBI Phmmer API' ], toggle=['DATABASE', 'close', ['PDB'] ] )
        self.createLine( [ 'advice', 'Alternative Alphafold database fasta sequence file (default uses internal CCP4 version)' ], toggle=['USEAPI', 'open', [False] ] )
        self.createLine( [ 'widget', '-browseDb', True, 'AFDBSEQDB', 'tip', 'Alternative Alphafold sequences fasta file' ], toggle=['USEAPI', 'open', [False] ] )
        self.closeSubFrame()

        self.createLine( [ 'advice', 'Number of cores for Phmmer (maximum=%d)' % MAXPROC ], toggle=['USEAPI', 'open', [False] ]  )
        self.createLine( [ 'widget', 'NPROC' ], toggle=['USEAPI', 'open', [False] ]  )

        self.createLine(['subtitle', 'Do sequence classification?', 'widget', 'DO_CLASSIFY'])
        self.closeSubFrame()
        # self.createLine(['advice',
        #                  'Coiled-coil and transmembrane classification require locally installed programs (Deepcoil & TMHMM)'])
        # self.createLine(['advice',
        #                  'DeepCoil can be downloaded from: <a href="https://github.com/labstructbioinf/DeepCoil">https://github.com/labstructbioinf/DeepCoil</a>'])
        # self.createLine(['advice',
        #                  'TMHMM can be downloaded from: <a href="https://github.com/dansondergaard/tmhmm.py">https://github.com/dansondergaard/tmhmm.py</a>'])
