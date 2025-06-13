"""
    modelASUCheck_report.py: CCP4 GUI Project

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

class modelASUCheck_report(Report):
    # Specify which gui task and/or pluginscript this applies to
    TASKNAME = 'modelASUCheck'
    RUNNING = False
    def __init__(self,xmlnode=None,jobInfo={},jobStatus=None,**kw):
        Report. __init__(self,xmlnode=xmlnode,jobInfo=jobInfo, jobStatus=jobStatus, **kw)
        clearingDiv = self.addDiv(style="clear:both;")
        xmlPath = './/SequenceAlignment/Alignment'
        xmlNodes = xmlnode.findall(xmlPath)
        self.addAlignReport(xmlNodes,True)

    def addAlignReport(self,xmlNodes,initiallyOpen=False,parent=None):
        if parent is None:
            parent = self
        
        if len(xmlNodes)>0:
          clearingDiv = parent.addDiv(style="clear:both;")
          alignFold = parent.addFold(label='Sequence alignment information',brief='Seq. Align. Info.', initiallyOpen=initiallyOpen)
          k = 60
          for alignNode in xmlNodes:
                chainID = alignNode.findall("ChainID")[0].text
                match_count = alignNode.findall("match_count")[0].text
                identity = alignNode.findall("identity")[0].text
                cigar = alignNode.findall("CIGAR")[0].text
                align_1 = alignNode.findall("align_1")[0].text
                align_2 = alignNode.findall("align_2")[0].text
                align_match = alignNode.findall("align_match")[0].text
                labelDiv = alignFold.addDiv(style='border:0px solid black; width:700px; overflow:auto;')
                labelDiv.append("<b>Chain: "+chainID+"</b>")
                tableText = "<table>\n"
                tableText += "<tr>"
                tableText += "<th>Match count</th>"
                tableText += "<td>"+match_count+"</td>"
                tableText += "</tr>"
                tableText += "<tr>"
                tableText += "<th>Identity</th>"
                tableText += "<td>"+identity+"</td>"
                tableText += "</tr>"
                tableText += "<tr>"
                tableText += "<th>CIGAR</th>"
                tableText += "<td>"+cigar+"</td>"
                tableText += "</tr>"
                tableText += "</table>\n"
                alignFold.append(tableText)
                pre_text = ""

                for i in range(0, len(align_match), k):
                     pre_text += str(i+1)
                     pre_text += (60-len(str(i)+str(i+k)))*' '
                     pre_text += str(min(i+k,len(align_match)))
                     pre_text += "\n"
                     pre_text += align_1[i:i+k] + "\n"
                     pre_text += align_match[i:i+k] + "\n"
                     pre_text += align_2[i:i+k] + "\n"
                     pre_text += "\n"

                alignFold.append('<pre>'+pre_text+'</pre>')
    
