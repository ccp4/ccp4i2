from ccp4i2.report import Report
from ccp4i2.wrappers.ShelxCDE.script import ShelxCE_report


class ShelxCECompareHands_report(Report):
    TASKNAME = 'ShelxCECompareHands'
    RUNNING = True
    SEPARATEDATA=True
    
    def __init__(self, *args, **kws):
        super(ShelxCECompareHands_report,self).__init__(*args, **kws)
        if kws.get('jobStatus','nooutput').lower() =='nooutput': return
        self.addDiv(style='clear:both')
        self.defaultReport(parent=self)

    def defaultReport(self, parent=None):
        if parent is None: parent = self
        
        choiceNodes = self.xmlnode.findall('HandChosen')
        for choiceNode in choiceNodes:
            result = choiceNode.text
            if result == 'Inverted':
                parent.addText(text='The inverted hand structure was chosen. - necessitates a change of space group enantiomorph in some space groups',style='color:red; font-color:red; font-size:125%;')
            else:
                parent.addText(text='The original hand structure was chosen.',style='color:blue; font-color:red; font-size:125%;')
        parent.append('<br/>')
        firstHandDiv = parent.addDiv()
        firstHandDiv.addText(text='Original HA structure',style='color:blue; font-color:red; font-size:125%;')
        firstHandNodes = self.xmlnode.findall('FirstHand/ShelxCE')
        for firstHandNode in firstHandNodes:
            shelxCEReport = ShelxCE_report.ShelxCE_report(xmlnode=firstHandNode,jobStatus=self.jobStatus)
            shelxCEReport.shelXEReport(parent=firstHandDiv, initiallyOpen=True,idRoot="FirstHand")
        
        secondHandDiv = parent.addDiv()
        secondHandDiv.addText(text='Inverted HA structure',style='color:blue; font-color:red; font-size:125%;')
        secondHandNodes = self.xmlnode.findall('SecondHand/ShelxCE')
        for secondHandNode in secondHandNodes:
            shelxCEReport = ShelxCE_report.ShelxCE_report(xmlnode=secondHandNode,jobStatus=self.jobStatus)
            shelxCEReport.shelXEReport(parent=secondHandDiv, initiallyOpen=True,idRoot="SecondHand")
