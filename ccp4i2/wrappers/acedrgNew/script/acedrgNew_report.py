from report.CCP4ReportParser import *
import sys
import math

class acedrgNew_report(Report):
    # Specify which gui task and/or pluginscript this applies to
    TASKNAME = 'Acedrg'
    RUNNING = True
    
    def __init__(self, *args, **kws):
        '''
            kws['additionalJsFiles'] = ['jsme/jsme.nocache.js']
            kws['additionalScript'] = 'function jsmeOnLoad() { jsmeApplet = new JSApplet.JSME("jsme_container", "380px", "340px");'
            '''
        Report.__init__(self, *args, **kws)
        if self.jobStatus is None or self.jobStatus.lower() == 'nooutput': return
        self.defaultReport()
    
    def defaultReport(self, parent=None):
        if parent is None: parent = self
        parent.addText(text='2D representations', style='font-size:110%;')
        div2D = parent.addDiv()
        self.analyseGeometry()
    
    def analyseGeometry(self,parent=None):
        if parent is None: parent=self
        categoryNameMap = {'_chem_comp_bond':'Bonds','_chem_comp_angle':'Angles','_chem_comp_tor':'Torsions','_chem_comp_plane_atom':'Planes', '_chem_comp_chir':'Chirals'}
        columnNameMap={'atom_id':'Atom','atom_id_1':'Atom 1','atom_id_2':'Atom 2','atom_id_3':'Atom 3','atom_id_4':'Atom 4','value_dist':'Dist', 'value_dist_esd':'Sigma','value_angle':'Angle','value_angle_esd':'Sigma','period':'Period','type':'Type','plane_id':'Plane', 'atom_id_centre':'Centre','volume_sign':'Sign'}
        for categoryName in categoryNameMap:
            examples = self.xmlnode.findall('.//Acedrg/Geometry/'+categoryName)
            if len(examples) > 0:
                folder = parent.addFold(label=categoryNameMap[categoryName])
                table = folder.addTable()
                for child in examples[0]:
                    dataNodes = self.xmlnode.findall('.//Acedrg/Geometry/'+categoryName+'/'+child.tag)
                    data = [dataNode.text for dataNode in dataNodes]
                    table.addData(title=columnNameMap[child.tag],data=data)

