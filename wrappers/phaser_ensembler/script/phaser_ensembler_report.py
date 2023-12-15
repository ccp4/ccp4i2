from report import CCP4ReportParser
import base64

class phaser_ensembler_report(CCP4ReportParser.Report):
    # Specify which gui task and/or pluginscript this applies to
    TASKNAME = 'phaser_ensembler'
    RUNNING = False

    def __init__(self,*args,**kws):
        super(phaser_ensembler_report,self).__init__(*args, **kws)
        self.addDiv(style='clear:both;')
        logFold = self.addFold(label='Log stream from phaser.ensembler', initiallyOpen = False)
        logPre = logFold.addPre()
        logPre.text = base64.b64decode(self.xmlnode.findall('LOGTEXT')[0].text)
