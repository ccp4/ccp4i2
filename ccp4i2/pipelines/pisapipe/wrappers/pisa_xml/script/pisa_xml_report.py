from report.CCP4ReportParser import *
import sys
import math

class pisa_xml_report(Report):
    # Specify which gui task and/or pluginscript this applies to
    TASKNAME = 'pisa_xml'
    RUNNING = True

    def __init__(self,*args,**kws):
        Report.__init__(self, *args, **kws)
        if self.jobStatus is None or self.jobStatus.lower() == 'nooutput': return
        self.defaultReport()

    def defaultReport(self,parent=None):
        if parent is None: parent=self
        self.addAssemblies(parent=parent)
        self.addInterfaces(parent=parent)

    def addAssemblies(self, parent=None):
        if parent is None: parent = self
        assemblies = self.xmlnode.findall('.//assembly')
        if len(assemblies) > 0:
            assemblyFold = parent.addFold(label='Assemblies', initiallyOpen=True)
            assemblyTable = assemblyFold.addTable(select = './/assembly',downloadable=True)
            for tag,title,expr in [('id','Id',"'{:4.0f}'.format(float(x))"),
                                   ('formula','Formula',None),
                                   ('composition','Composition',None),
#                              ('size','Size',"'{:4.0f}'.format(float(x))"),
#                              ('mmsize','mmsize',"'{:4.0f}'.format(float(x))"),
#                             ('freesize','freesize',"'{:4.0f}'.format(float(x))"),
                              ('diss_energy','Diss E.',"'{:10.1f}'.format(float(x))"),
                              ('diss_energy_0','Diss E.(0)',"'{:10.1f}'.format(float(x))"),
                              ('asa','ASA',"'{:10.1f}'.format(float(x))"),
                              ('bsa','BSA',"'{:10.1f}'.format(float(x))"),
                              ('entropy','Entropy',"'{:10.1f}'.format(float(x))"),
                              ('entropy_0','Entropy(0)',"'{:10.1f}'.format(float(x))"),
                              ('diss_area','Diss. area',"'{:10.1f}'.format(float(x))")]:
                assemblyTable.addData(title=title,select=tag,expr=expr)

    def addInterfaces(self, parent=None):
        if parent is None: parent = self
        interfaces = self.xmlnode.findall('.//Interfaces/pdb_entry/interface')
        if len(interfaces) > 0:
            interfaceFold = parent.addFold(label='Interfaces', initiallyOpen=True)
            interfaceTable = interfaceFold.addTable(select = './/Interfaces/pdb_entry/interface',downloadable=True)
            mol1IdData = []
            mol2IdData = []
            for interface in interfaces:
                molIDNodes = interface.findall('molecule/chain_id')
                mol1IdData += [molIDNodes[0].text]
                mol2IdData += [molIDNodes[1].text]
            interfaceTable.addData(title='Chain 1', data=mol1IdData)
            interfaceTable.addData(title='Chain 2', data=mol2IdData)
            for tag,title,expr in [('id','Id',"'{:4.0f}'.format(float(x))"),
                              ('type','Type',"'{:4.0f}'.format(float(x))"),
                              ('n_occ','n_occ',"'{:4.0f}'.format(float(x))"),
                              ('salt-bridges/n_bonds','Salt<br/>bridges',"'{:4.0f}'.format(float(x))"),
                              ('h-bonds/n_bonds','H-<br/>bonds',"'{:4.0f}'.format(float(x))"),
                              ('ss-bonds/n_bonds','SS-<br/>bonds',"'{:4.0f}'.format(float(x))"),
                               ('cov-bonds/n_bonds','Cov.<br/>bonds',"'{:4.0f}'.format(float(x))"),
                              ('int_area','Area',"'{:10.1f}'.format(float(x))"),
                             ('int_solv_en','Solvation<br/>Energy',"'{:10.1f}'.format(float(x))"),
                                   ('css','css',"'{:10.1f}'.format(float(x))"),
                              ('pvalue','p-value',"'{:10.1f}'.format(float(x))"),
                                   ('overlap','overlap',None)]:
                interfaceTable.addData(title=title,select=tag,expr=expr)

            for interface in interfaces:
                idNodes = interface.findall('id')
                typeNodes = interface.findall('type')
                if len(idNodes)>0 and len(typeNodes)>0:
                    folderLabel = '\tType: ' + typeNodes[0].text + ' id: ' + idNodes[0].text
                    specificInterfaceFold = interfaceFold.addFold(label = folderLabel, xmlnode=interface)
                    molecules = interface.findall('molecule')
                    moleculeTable = specificInterfaceFold.addTable(xmlnode=interface,downloadable=True)
                    for tag,title,expr in [('molecule/id','Id',None),
                                           ('molecule/chain_id','Chain',None),
                                           ('molecule/symop','Sym. Op.',None),
                                           ('molecule/int_natoms','No. atoms',None),
                                           ('molecule/int_nres','No. res.',None),
                                           ('molecule/int_area','Area',"'{:10.1f}'.format(float(x))"),
                                           ('molecule/int_solv_en','Solvation<br/>Energy',"'{:10.1f}'.format(float(x))"),
                                           ('molecule/pvalue','p-value',"'{:10.1f}'.format(float(x))")
                                           ]:
                        moleculeTable.addData(title=title,select=tag,expr=expr)
'''<id>1</id>
    <type>1</type>
    <n_occ>2</n_occ>
    <int_area>1219.4288508</int_area>
    <int_solv_en>-11.132901313</int_solv_en>
    <pvalue>0.23535434822</pvalue>
    <stab_en>-19.281643582</stab_en>
    <css>0.64823734571</css>
    <overlap>No</overlap>
    <x-rel>No</x-rel>
    <fixed>No</fixed>
'''
'''
    <molecule>
    <id>1</id>
    <chain_id>C</chain_id>
    <class>Protein</class>
    <symop_no>1</symop_no>
    <symop>x,y,z</symop>
    <cell_i>0</cell_i>
    <cell_j>0</cell_j>
    <cell_k>0</cell_k>
    <rxx>1</rxx>
    <rxy>0</rxy>
    <rxz>0</rxz>
    <tx>0</tx>
    <ryx>0</ryx>
    <ryy>1</ryy>
    <ryz>0</ryz>
    <ty>0</ty>
    <rzx>0</rzx>
    <rzy>0</rzy>
    <rzz>1</rzz>
    <tz>0</tz>
    <int_natoms>75</int_natoms>
    <int_nres>22</int_nres>
    <int_area>717.5508296</int_area>
    <int_solv_en>-2.1232093624</int_solv_en>
    <pvalue>0.45669242376</pvalue>
'''

'''<serial_no>0</serial_no>
    <id>1</id>
    <size>1</size>
    <mmsize>1</mmsize>
    <freesize>1</freesize>
    <diss_energy>0</diss_energy>
    <diss_energy_0>0</diss_energy_0>
    <asa>6724.6279848</asa>
    <bsa>0</bsa>
    <entropy>0</entropy>
    <entropy_0>0</entropy_0>
    <diss_area>0</diss_area>
    <int_energy>0</int_energy>
    <n_uc>0</n_uc>
    <n_diss>1</n_diss>
    <symNumber>1</symNumber>
    <formula/>
    <composition>C</composition>
    '''
