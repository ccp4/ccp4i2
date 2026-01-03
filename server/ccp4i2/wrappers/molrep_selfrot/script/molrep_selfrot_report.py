from ccp4i2.report import Report


class molrep_selfrot_report(Report):
    TASKNAME = 'molrep_selfrot'
    RUNNING = False
    def __init__(self,xmlnode=None,jobInfo={},jobStatus=None,**kw):
        Report.__init__(self,xmlnode=xmlnode,jobInfo=jobInfo,**kw)
        
        if jobStatus is None or jobStatus.lower() == 'nooutput': return
        self.defaultReport()

    def defaultReport(self, parent=None):
        if parent is None: parent=self
        
        for structureFactorNode in self.xmlnode.findall('.//StructureFactors'):
            parent.append('<br/>')
            parent.append('<br/>')
            parent.addText(text = 'Conclusion of search for anisotropy:',style='font-size:120%;')
            parent.addText(xmlnode=structureFactorNode, select='INFO',style='font-size:120%;')
            dataTable = parent.addTable(xmlnode = self.xmlnode, select='StructureFactors')
            headings = {'LowResProvided':'Dmax',
                'HighResProvided':'Dmin',
                'Completeness':'Completeness',
                'BOverall':'B-factor',
                'OpticalHighRes':'Optical Dmin',
                'EigenValueRatioH':'Aniso H',
                'EigenValueRatioK':'Aniso K',
                'EigenValueRatioL':'Aniso L'}
            for heading in headings:
                dataTable.addData(select=heading, title=headings[heading])
        
        for pattersonNode in self.xmlnode.findall ('.//Patterson'):
            #Pull out conclusion of search for translational symmetry
            parent.append('<br/>')
            parent.append('<br/>')
            parent.addText(text = 'Conclusion of search for translational symmetry:',style='font-size:120%;')
            parent.addText(xmlnode=pattersonNode, select='INFO',style='font-size:120%;')
            
            #Provide list of self Patterson peaks as a table
            parent.append('<br/>')
            pattersonFold=parent.addFold(label='Self Patterson peaks',initiallyOpen=False,style='font-size:110%;')
            peakTable = pattersonFold.addTable(xmlnode = pattersonNode, select='Peak',title='Self-Patterson peaks')
            headings = 'No IX  IY  IZ  Xfrac  Yfrac  Zfrac  Xort   Yort    Zort   Dens Dens_sigma'.split()
            for heading in headings:
                peakTable.addData(select=heading, title=heading)

        for rotationNode in self.xmlnode.findall ('.//SelfRotation'):
            #Provide list of symmetry-expanded rotation funciton peaks as a table
            parent.append('<br/>')
            selfRotationFold=parent.addFold(label='Self Rotation peaks',initiallyOpen=False,style='font-size:110%;')
            peakTable = selfRotationFold.addTable(xmlnode = rotationNode, select='Peak',title='Self-Patterson peaks')
            headings=  'No theta    phi     chi    alpha    beta   gamma      Rf    Rf_sigma'.split()
            for heading in headings:
                peakTable.addData(select=heading, title=heading)

        parent.append('<br/>')
        parent.addText(text='To view the rotation function, right click on the icon of the "Self-rotation function" data file below, and select "View1"',style='font-size:120%;')
