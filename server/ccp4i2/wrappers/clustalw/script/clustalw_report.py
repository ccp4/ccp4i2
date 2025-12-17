from ccp4i2.report import Report


class clustalw_report(Report):
    TASKNAME = 'clustalw'
    RUNNING = False
    def __init__(self,xmlnode=None,jobInfo={},jobStatus=None,**kw):
        Report. __init__(self,xmlnode=xmlnode,jobInfo=jobInfo,jobStatus=jobStatus, **kw)
        
        if jobStatus is not None and jobStatus.lower() == 'nooutput': return
        
        self.defaultReport(parent=self)

    def defaultReport(self, parent= None):
        self.addDiv(style='clear:both;')
        if parent is None: parent= self

        for alignmentNode in self.xmlnode.findall('Alignment'):
                parent.addPre(text = alignmentNode.text)

        for node in self.xmlnode.findall('Statistics'):
            parent.addPre(text = node.text)
