from ccp4i2.baselayer import QtCore
from ccp4i2.pipelines.phaser_pipeline.wrappers.phaser_EP_AUTO.script import phaser_EP_AUTO_gui


class phaser_EP_gui(phaser_EP_AUTO_gui.phaser_EP_AUTO_gui):
    # Subclass CTaskWidget to give specific task window
    TASKNAME = 'phaser_EP'
    TASKVERSION = 0.1
    TASKMODULE='expt_phasing'
    TASKTITLE='SAD phasing - PHASER'
    DESCRIPTION = '''Complete a heavy atom model and calculate phases'''
    ERROR_CODES = {200:{'description':'FreeR not set'}}
    RANK=1
    WHATNEXT=['modelcraft']

    def __init__(self,*args,**kws):
        super(phaser_EP_gui,self).__init__(*args, **kws)
    
    def drawContents(self):
        self.setProgramHelpFile('phaser')
        self.drawPhaserFrontPage()
        self.drawPhaserKeywordsFolder()

    def drawPhaserFrontPage(self):
        self.openFolder(folderFunction='inputData')
        self.createLine( ['subtitle', 'Experimental data','tip','Phaser uses this for maximum likelihood and to suggest number of copies in the asymmetric unit' ])
        self.openSubFrame(frame=True)
        self.createLine ( ['label', '\n '])
        self.createLine ( [ 'widget', 'F_SIGF' ] )
        self.container.inputData.F_SIGF.dataChanged.connect( self.getWavelength)
        if self.container.inputData.F_SIGF.isSet(): self.getWavelength()
        self.createLine (['label','Resolution range','stretch','widget','RESOLUTION_LOW','widget','RESOLUTION_HIGH','label','Wavelength','widget','WAVELENGTH'])
        self.closeSubFrame()
        self.createLine( ['subtitle', 'Heavy atom search, known sites, phasing model, or phases','tip','known sites will be used as start point for competing the heavy atom model. Partial model or map will be used to start to find sites.' ])
        self.openSubFrame(frame=True)
        self.createLine(['label','Find sites, provide sites, partial model or phases as:','widget','PARTIALMODELORMAP'])
        self.container.inputData.PARTIALMODELORMAP.dataChanged.connect(self.handlePARTIALMODELORMAP)
        self.handlePARTIALMODELORMAP()
        self.createLine ( [ 'widget','-guiLabel','Partial protein model', 'XYZIN_PARTIAL' ], toggle=['PARTIALMODELORMAP','open',['MODEL']])
        self.createLine ( [ 'widget','-guiLabel','Partial protein model', 'MAPCOEFF_PARTIAL' ], toggle=['PARTIALMODELORMAP','open',['MAP']])
        self.createLine ( [ 'widget','-guiLabel','Partial HA model', '-browseDb','True','XYZIN_HA' ], toggle=['PARTIALMODELORMAP','open',['NONE']] )
        self.createLine(['label','Atom type to seek','stretch','widget','SFAC'], toggle=['PARTIALMODELORMAP','open',['SEARCH']] )
        self.createLine(['label','Number to find','stretch','widget','FIND'], toggle=['PARTIALMODELORMAP','open',['SEARCH']] )
        self.createLine(['label','Number of tries','stretch','widget','NTRY'], toggle=['PARTIALMODELORMAP','open',['SEARCH']] )
        self.closeSubFrame()
        self.createLine( ['subtitle', 'Composition','tip','Phaser uses this for maximum likelihood and to suggest number of copies in the asymmetric unit' ])
        self.openSubFrame(frame=True)
        self.setMenuText('COMP_BY',{ 'DEFAULT':'Use protein average solvent content',
                         'MW':'Provide an estimate of the molecular weight of protein and nucleic acid in ASU',
                         'ASU':'Provide a full specification of the ASU content by sequence'
                         })
        self.createLine ( ['tip','Phaser uses this for maximum likelihood and to suggest number of copies in the asymmetric unit','label', 'For estimating asymmetric unit contents:','widget','COMP_BY'])
        self.container.inputData.COMP_BY.dataChanged.connect(self.handleCOMP_BY)
        self.handleCOMP_BY()
        self.createLine ( ['label', '\n '], toggle=['COMP_BY','open',['DEFAULT']])
        self.createLine ( [ 'widget', 'ASUFILE' ], toggle=['COMP_BY','open',['ASU']])
        self.createLine ( [ 'label','Molecular weight (Da) of protein in the ASU','stretch','widget', 'ASU_PROTEIN_MW' ], toggle=['COMP_BY','open',['MW']])
        self.createLine ( [ 'label','Molecular weight (Da) of nucleic acid','stretch','widget', 'ASU_NUCLEICACID_MW' ], toggle=['COMP_BY','open',['MW']])
        self.closeSubFrame()
        self.createLine ( [ 'label','Cycles of heavy atom model completion','widget', 'LLGC_CYCLES' ] )
        self.createLine(['label','Search for pure anomalous scatterers','widget','PURE_ANOMALOUS'])
        self.createLine(['widget','-title','Elements to search for','ELEMENTS'])
        self.createLine(['subtitle','Run density modification and model building','Optionally run density modification with Parrot and automated model building'])
        self.openSubFrame(frame=True)
        self.createLine(['widget','RUNPARROT','label','Run density modification'])
        self.createLine(['widget','RUNBUCCANEER','label','Run ','widget','BUCCANEER_ITERATIONS','label','cycle(s) of automated model building'])
        self.connect(self.container.controlParameters.RUNBUCCANEER,QtCore.SIGNAL('dataChanged'),self.handleRUNBUCCANEER)
        self.createLine(['widget','FREERFLAG'], toggle=['RUNBUCCANEER', 'open', [True]])
        self.closeSubFrame()

    @QtCore.Slot()
    def handleRUNBUCCANEER(self):
        if self.container.controlParameters.RUNBUCCANEER:
            self.container.inputData.FREERFLAG.setQualifiers({'allowUndefined':False,'mustExist':True})
        else:
            self.container.inputData.FREERFLAG.setQualifiers({'allowUndefined':True,'mustExist':False})
            self.container.inputData.FREERFLAG.unSet()
        self.validate()

    def taskValidity(self):
        from ccp4i2.core import CCP4ErrorHandling
        rv = CCP4ErrorHandling.CErrorReport()
        # Check the space group is same in both input Mini-MTZ files 
        if self.container.controlParameters.RUNBUCCANEER:
            if not self.container.inputData.FREERFLAG.isSet():
                rv.append(self.__class__,200)
        return rv
