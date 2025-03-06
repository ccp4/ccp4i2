from __future__ import print_function

"""
    ZZPipelineNameZZ_report.py: CCP4 GUI Project
    
    This library is free software: you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public License
    version 3, modified in accordance with the provisions of the
    license to address the requirements of UK law.
    
    You should have received a copy of the modified GNU Lesser General
    Public License along with this library.  If not, copies may be
    downloaded from http://www.ccp4.ac.uk/ccp4license.php
    
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.
    """

from report.CCP4ReportParser import Report
import sys
#from lxml import etree
import xml.etree.ElementTree as etree

class ZZPipelineNameZZ_report(Report):
    # Specify which gui task and/or pluginscript this applies to
    TASKNAME = 'ZZPipelineNameZZ'
    RUNNING = False
    def __init__(self,xmlnode=None,jobInfo={},jobStatus=None,**kw):
        Report. __init__(self, xmlnode=xmlnode, jobInfo=jobInfo, jobStatus=jobStatus, **kw)
        clearingDiv = self.addDiv(style="clear:both;")
        self.addDefaultReport(self)
        
    def addDefaultReport(self, parent=None):
        if parent is None: parent=self
        for ZZFirstPluginNameZZNode in self.xmlnode.findall(".//ZZFirstPluginNameZZ"):
            try:
                cycleNode = ZZFirstPluginNameZZNode.findall("Cycle")[0]
            except:
                print("Missing cycle")
            try:
                logTextNode = ZZFirstPluginNameZZNode.findall("LogText")[0]
            except:
                print("Missing logText")
            try:
                newFold = parent.addFold(label="Log text for iteration "+cycleNode.text, initiallyOpen=True)
            except:
                print("Unable to make fold")
            try:
                newFold.addPre(text = logTextNode.text)
            except:
                print("Unable to add Pre")
