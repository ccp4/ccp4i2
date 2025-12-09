from ccp4i2.baselayer import QtGui, QtWidgets,QtCore

from qtgui.CCP4TaskWidget import CTaskWidget
from core.CCP4ErrorHandling import *
import functools

#-------------------------------------------------------------------
class matthews_gui(CTaskWidget):
#-------------------------------------------------------------------
    
    # Subclass CTaskWidget to give specific task window
    TASKNAME = 'matthews'
    TASKVERSION = 0.1
    TASKMODULE=[ 'data_reduction', 'expt_data_utility' ]
    TASKTITLE='Estimate AU content'
    EDITOR = True  
    DESCRIPTION='''Estimate number of molecules in the asymmetric unit and solvent content (Matthews_coeff)'''
    RANK = 2
    MAINTAINER = 'liz.potterton@york.ac.uk'
    
    def drawContents(self):
        
        self.openFolder(folderFunction='inputData',followFrom=False)
        self.createLine( ['subtitle' , 'Cell parameters taken from reflection data', 'Select a reflection file containing the cell parameters' ] )
        self.openSubFrame(frame=[True])
        self.createLine( [ 'widget', 'HKLIN' ] )
        self.closeSubFrame()

        self.createLine( ['subtitle' , 'Calculate molecular weight from..','Define the contents of the asymmetric unit'] )
        self.openSubFrame(frame=[True])
        self.group = QtWidgets.QButtonGroup(self)
        
        widget = QtWidgets.QRadioButton('.. the sequence composition',self)
        widget.setToolTip('Enter sequence to define AU content')
        widget.setChecked(True)
        self.group.addButton(widget,1)
        line = self.createLine()
        line.layout().addWidget(widget)
        self.createLine ( [ 'widget', '-title','Contents of biological unit', 'ASUIN' ] )
        #self.createLine( [ 'widget', 'SEQIN' ] )
        self.container.inputData.ASU_COMPONENTS.dataChanged.connect(functools.partial(self.updateMode,'asu_components'))
        widget.clicked.connect(self.runAnalysis)
        
        widget = QtWidgets.QRadioButton('.. number of residues',self)
        widget.setToolTip('Enter number of residues to define AU content')
        self.group.addButton(widget,2)
        line = self.createLine( ['tip','Number of residues','widget', 'NRES' ] )
        line.layout().insertWidget(0,widget)
        self.container.inputData.NRES.dataChanged.connect(functools.partial(self.updateMode,'nres'))
        widget.clicked.connect(self.runAnalysis)
        
        widget = QtWidgets.QRadioButton('.. molecular weight',self)
        widget.setToolTip('Enter molecular weight to define AU content')
        line = self.createLine( ['tip','Molecular weight','widget', 'MOLWT' ] )
        self.group.addButton(widget,3)
        line.layout().insertWidget(0,widget)
        self.container.inputData.MOLWT.dataChanged.connect(functools.partial(self.updateMode,'molwt'))
        widget.clicked.connect(self.runAnalysis)
        self.closeSubFrame()

        # KJS : remove the results widget from the interface to avoid issues with display.
        self.resultWidget = QtWidgets.QTextEdit(self)
        self.resultWidget.setReadOnly(True)
        self.resultWidget.hide()
        
        self.container.inputData.HKLIN.dataChanged.connect(self.runAnalysis)
        self.container.inputData.ASUIN.dataChanged.connect(self.runAnalysis)
        self.group.buttonClicked[int].connect(self.handleModeChange)

        self.runAnalysis()

    def handleModeChange(self,mode):
        self.container.inputData.ASUIN.setQualifiers({'allowUndefined' : (mode!=1)})
        self.getWidget('ASUIN').validate()
        
    @QtCore.Slot(str)
    def updateMode(self,mode):
        #print 'matthews_gui.updateMode',mode
        if self.container.inputData.get(mode.upper()).isSet():
            self.container.inputData.MODE = mode
            self.group.button(self.container.inputData.MODE.qualifiers('enumerators').index(mode)+1).setChecked(True)
            self.runAnalysis()
         
    @QtCore.Slot()
    def runAnalysis(self):
        if not self.isEditable():
            return
        self.container.inputData.MODE = self.container.inputData.MODE.qualifiers('enumerators')[self.group.checkedId()-1]
        #print 'runAnalysis',self.container.inputData.MODE,self.container.inputData.ASU_COMPONENTS.isSet()
        if self.container.inputData.MODE == 'asu_components' and not self.container.inputData.ASUIN.isSet():
            self.resultWidget.clear()
            return
        elif self.container.inputData.MODE == 'molwt' and not self.container.inputData.MOLWT.isSet():
            self.resultWidget.clear()
            return
        elif self.container.inputData.MODE == 'nres' and not self.container.inputData.NRES.isSet():
            self.resultWidget.clear()
            return

        if not  self.container.inputData.HKLIN.isSet():
            self.resultWidget.clear()
            return
        #print 'matthews.runAnalysis mode',self.container.inputData.MODE
       
        text = '<html><body>'
        if self.container.inputData.MODE == 'asu_components':
            text = text + '<p>Calculating molecular weight for composition\n'
            for obj in self.container.inputData.ASUIN.fileContent.seqList:
                if self.container.inputData.ASUIN.isSelected(obj):
                    text = text + '<br>'+str(obj.nCopies) + ' * '+ str(obj.name) + '</br>'
                    text = text + '<br>'+ 'Which has molecular weight: '+ str(obj.molecularWeight()) + '</br>'
            text +=  '<br>Total sequence weight: '+str(self.container.inputData.ASUIN.molecularWeight()) + '</br>'
        elif self.container.inputData.MODE == 'nres':
            text = text + '<p>Estimated molecular weight for ' + str(self.container.inputData.NRES) + ' residues: '+str(self.container.inputData.NRES*110)
        else:
            text = text + '<p>Using given molecular weight: ' + str(self.container.inputData.MOLWT)
          
        try:
            self.container.inputData.HKLIN.loadFile()
            if self.container.inputData.MODE == 'asu_components':
                molWt = self.container.inputData.ASUIN.molecularWeight()
                #print 'molWt to matthewsCoeff',molWt
                rv = self.container.inputData.HKLIN.fileContent.matthewsCoeff(molWt=molWt)
            elif self.container.inputData.MODE == 'nres':          
                rv = self.container.inputData.HKLIN.fileContent.matthewsCoeff(nRes=self.container.inputData.NRES)
            else:
                rv = self.container.inputData.HKLIN.fileContent.matthewsCoeff(molWt=self.container.inputData.MOLWT)
        except:
            rv = {}

        #print 'runAnalysis',rv
        ktxt = text + "</p>" + "</body></html>"
        vol = rv.get('cell_volume','Unkown')
        nmol=[]
        solv = []
        matt=[]
        prob=[]
        if vol == 'Unkown':
            text = text + '<p>Cell volume = Unknown</p>'
        else:
            text = text + '<p>Cell volume = {0:.1f}<p>\n'.format(float(vol)) 
        text = text +    '<table><tr><th>  Nmol  </th><th>  %solvent  </th><th>  Matthews  </th><th>  prob(Matthews)  </th></tr>\n'
        for result in rv.get('results',[]):
            text = text + '<tr><td>  {0}  </td><td>  {1:.2f}  </td><td>  {2:.2f}  </td><td>  {3:.2f}  </td></tr>'.format(result.get('nmol_in_asu'),result.get('percent_solvent'),result.get('matth_coef'),result.get('prob_matth')) + '\n'
            nmol.append(result.get('nmol_in_asu'))
            solv.append(result.get('percent_solvent'))
            matt.append(result.get('matth_coef'))
            prob.append(result.get('prob_matth'))
        text = text + '</table></body></html>'
        self.resultWidget.setReadOnly(False)
        self.resultWidget.setHtml(text)
        self.resultWidget.setReadOnly(True)


        from report import CCP4ReportParser
        from core import CCP4Modules
        report = CCP4ReportParser.Report()
       
        # Fixed CCP4ReportParser so this should actually work now
        report.append("<h1>Matthews Coefficient Report</h1>")
        report.append(ktxt)
        if vol == 'Unkown':
            report.addText(text='Cell volume = Unknown')
        else:
            report.addText(text='Cell volume = {0:.1f}\n'.format(float(vol)))
        tab = report.addTable()
        tab.addData(title='Nmol',data=nmol)
        tab.addData(title='%solvent',data=solv)
        tab.addData(title='Matthews',data=matt)
        tab.addData(title='prob(Matthews)',data=prob)
        report.as_html_file(fileName=CCP4Modules.PROJECTSMANAGER().makeFileName(jobId=self.jobId(),mode='REPORT'))


