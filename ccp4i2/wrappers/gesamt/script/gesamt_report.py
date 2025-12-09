from ccp4i2.report.CCP4ReportParser import *
import sys

class gesamt_report(Report):
    # Specify which gui task and/or pluginscript this applies to
    TASKNAME = 'gesamt'
    RUNNING = False
    def __init__(self,xmlnode=None,jobInfo={},jobStatus=None,**kw):
        Report. __init__(self,xmlnode=xmlnode,jobInfo=jobInfo,**kw)
        
        if jobStatus is not None and jobStatus.lower() == 'nooutput':
            return

        textDiv = self.addDiv(style='font-size:110%;')

        transformationNodesCSV = self.xmlnode.findall('.//TransformationCSV')

        if len(transformationNodesCSV) > 0:
            textDiv.append('<br/>')
            textDiv.append('<br/>')
            textDiv.append('<br/>')
            transformationNodes = self.xmlnode.findall('.//Transformation')
            qNode = self.xmlnode.findall('.//TransformationCSV/Q_CSV')[0]
            rmsdNode = self.xmlnode.findall('.//TransformationCSV/RMSD_CSV')[0]
            nAlignNode = self.xmlnode.findall('.//TransformationCSV/NALIGN_CSV')[0]
            if len(transformationNodes) > 0:
                messageText = 'Superimposed '+transformationNodes[0].get('nRes')+' residues (sequence identity '+transformationNodes[0].get('seqid')+') with rms='+transformationNodes[0].get('rms')+' and Q='+transformationNodes[0].get('q')
                textDiv.addText(text=messageText)
                textDiv.append('<br/>')
            elif qNode is not None and nAlignNode is not None and rmsdNode is not None:
                messageText = 'Superimposed '+nAlignNode.text.strip()+' residues with rms='+rmsdNode.text.strip()+' and Q='+qNode.text.strip()
                textDiv.addText(text=messageText)
                textDiv.append('<br/>')
            if len(transformationNodes) > 0:
                messageText = 'Polar angles omega='+transformationNodes[0].get('omega')+', phi='+transformationNodes[0].get('phi')+', kappa='+transformationNodes[0].get('kappa')
                textDiv.addText(text=messageText)
                textDiv.append('<br/>')
                messageText = 'Euler angles alpha='+transformationNodes[0].get('alpha')+', beta='+transformationNodes[0].get('beta')+', gamma='+transformationNodes[0].get('gamma')
                textDiv.addText(text=messageText)
                textDiv.append('<br/>')
                messageText = 'Cartesian shift x='+transformationNodes[0].get('x')+', y='+transformationNodes[0].get('y')+', z='+transformationNodes[0].get('z')
                textDiv.addText(text=messageText)
                textDiv.append('<br/>')
            matrixNodes = self.xmlnode.findall('.//TransformationCSV/matrixCSV')
            if len(matrixNodes)==1:
                textDiv.append('<br/>')
                transFold = self.addFold(label='Transformation matrix', initiallyOpen=True)
            else:
                textDiv.append('<br/>')
                transFold = self.addFold(label='Transformation matrices', initiallyOpen=True)
            istr = 1
            for matrixNode in matrixNodes:
                if len(matrixNodes)>1:
                    messageText = 'Structure {0}'.format(istr)
                    transFold.addText(text=messageText)
                matrixVals = matrixNode.text.split(",")
                lt = "".join(["{:.5f}".format(float(x)).rjust(12) for x in matrixVals[:4]]) + "\n"
                lt += "".join(["{:.5f}".format(float(x)).rjust(12) for x in matrixVals[4:8]]) + "\n"
                lt += "".join(["{:.5f}".format(float(x)).rjust(12) for x in matrixVals[8:]]) + "\n"
                transFold.addPre(text = lt)
                istr += 1
            transFold.append('<br/>')

            alignFold = self.addFold(label='RMSD vs. residue', initiallyOpen=True)
#TODO - handle multiple chains properly
            progressGraph = alignFold.addFlotGraph(title="Inter-residue distances",select=".//TransformationCSV/PerResidue/Chain/Row",style="height:250px; width:600px;float:left;border:0px solid white;")
            progressGraph.addData(title="Residue_number",         select="iRes")
            progressGraph.addData(title="Distance", select="Distance")
            plot = progressGraph.addPlotObject()
            plot.append('title','Distances')
            plot.append('plottype','xy')
            plot.append('xintegral','true')
            for coordinate, colour in [(2,'blue')]:
                plotLine = plot.append('plotline',xcol=1,ycol=coordinate,colour=colour)
            self.addDiv(style='clear:both')

            equivalenceFold = self.addFold(label='List of equivalenced residues')
            equivalenceNodes = self.xmlnode.findall('.//TransformationCSV/PerResidue/Chain/Row/Equivalence')
            iRow = 0
            for equivalenceNode in equivalenceNodes:
                equivalenceFold.addPre(text = equivalenceNode.text+': Position '+str(iRow))
                iRow += 1
            
            return
        
        transformationNodes = self.xmlnode.findall('.//Transformation')

        if len(transformationNodes) == 0:
            messageText = 'No superpositions found : suggest you inspect the log file'
            textDiv.addText(text=messageText)
        elif len(transformationNodes) >1:
            messageText = 'Apparently more than 1 superposition found : report to developer (martin.noble@ncl.ac.uk)'
            textDiv.addText(text=messageText)
        else:
            textDiv.append('<br/>')
            textDiv.append('<br/>')
            textDiv.append('<br/>')
            messageText = 'Superimposed '+transformationNodes[0].get('nRes')+' residues (sequence identity '+transformationNodes[0].get('seqid')+') with rms='+transformationNodes[0].get('rms')+' and Q='+transformationNodes[0].get('q')
            textDiv.addText(text=messageText)
            textDiv.append('<br/>')
            messageText = 'Polar angles omega='+transformationNodes[0].get('omega')+', phi='+transformationNodes[0].get('phi')+', kappa='+transformationNodes[0].get('kappa')
            textDiv.addText(text=messageText)
            textDiv.append('<br/>')
            messageText = 'Euler angles alpha='+transformationNodes[0].get('alpha')+', beta='+transformationNodes[0].get('beta')+', gamma='+transformationNodes[0].get('gamma')
            textDiv.addText(text=messageText)
            textDiv.append('<br/>')
            messageText = 'Cartesian shift x='+transformationNodes[0].get('x')+', y='+transformationNodes[0].get('y')+', z='+transformationNodes[0].get('z')
            textDiv.addText(text=messageText)
            textDiv.append('<br/>')
            progressGraph = self.addFlotGraph(title="Inter-residue distances",select=".//Transformation/PerResidue/Row",style="height:250px; width:600px;float:left;border:0px solid white;")
            progressGraph.addData(title="Residue_number",         select="iRes")
            progressGraph.addData(title="Distance", select="Distance")
            plot = progressGraph.addPlotObject()
            plot.append('title','Distances')
            plot.append('plottype','xy')
            plot.append('xintegral','true')
            for coordinate, colour in [(2,'blue')]:
                plotLine = plot.append('plotline',xcol=1,ycol=coordinate,colour=colour)
            self.addDiv(style='clear:both')
            equivalenceFold = self.addFold(label='List of equivalenced residues')
            equivalenceNodes = self.xmlnode.findall('.//Transformation/PerResidue/Row/Equivalence')
            iRow = 0
            for equivalenceNode in equivalenceNodes:
                equivalenceFold.addPre(text = equivalenceNode.text+': Position '+str(iRow))
                iRow += 1
