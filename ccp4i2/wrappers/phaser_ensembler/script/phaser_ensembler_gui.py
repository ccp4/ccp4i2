from ....qtgui.CCP4TaskWidget import CTaskWidget


#-------------------------------------------------------------------
class phaser_ensembler_gui(CTaskWidget):
#-------------------------------------------------------------------

    # Subclass CTaskWidget to give specific task window
    TASKNAME = 'phaser_ensembler'
    TASKVERSION = 0.1
    TASKMODULE='bioinformatics'
    TASKTITLE='Build an ensemble for PHASER'
    DESCRIPTION='Compile assorted related structures into an ensemble for use in PHASER'
    MGDISPLAYFILES = ['XYZOUT']
    
    def __init__(self,*args, **kws):
        super(phaser_ensembler_gui,self).__init__(*args, **kws)
        
    def drawContents(self):
        self.drawInputsPage()
        self.drawKeywordsPage()
        
    def drawInputsPage(self):
        self.openFolder(folderFunction="inputData")
        self.createLine( ['subtitle', 'Coordinates' ])
        self.openSubFrame(frame=True)
        self.createLine ( [ 'widget', 'XYZIN_LIST' ] )
        from ....qtgui.CCP4ModelWidgets import CPdbDataFileView
        for pdbDataFileView in self.findChildren(CPdbDataFileView):
            pdbDataFileView.showAtomSelection()
        self.closeSubFrame()
        self.createLine( ['subtitle', 'Sequence identity to place in header','All models in the resulting ensemble will be annotated as sharing this percentage identity to the target' ])
        self.openSubFrame(frame=True)
        self.createLine ( [ 'widget', 'OVERRIDEID' ] )
        self.closeSubFrame()
        self.createLine( ['subtitle', 'Optional multiple sequence alignment' ])
        self.openSubFrame(frame=True)
        self.createLine ( [ 'widget', 'ALIGNIN' ] )
        self.closeSubFrame()

    def drawKeywordsPage(self):
        self.openFolder(folderFunction='keywords', title='Keywords')
        #Handle nesting
        nestedWords = {}
        for attr in self.container.keywords.dataOrder():
            levels = attr.split('__')
            baseLevel = nestedWords
            for iLevel, level in enumerate(levels):
                if iLevel == (len(levels)-1):
                    baseLevel[level] = attr
                elif level not in baseLevel:
                    baseLevel[level] = {}
                baseLevel = baseLevel[level]
    
        ignoreLevels = ['input']
        ignoreParameters = ['output__gui_output_dir','output__location','output__job_title','output__root','configuration__superposition__atoms']
        
        def unpackLevels(levels, indent = 0):
            for iLevel,level in enumerate(levels):
                if isinstance(levels[level], dict):
                    if level not in ignoreLevels:
                        self.createLine(['subtitle',level])
                        if indent == 0: self.openSubFrame(frame=True)
                        unpackLevels(levels[level], indent+1)
                        if indent == 0: self.closeSubFrame()
                else:
                    if levels[level] not in ignoreParameters:
                        self.createLine(['label',level,'stretch','widget',levels[level]])
        
        unpackLevels(nestedWords)
