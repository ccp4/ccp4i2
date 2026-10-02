from ccp4i2.report.CCP4ReportParser import Report
from ccp4i2.wrappers.validate_protein.script import validate_protein_report


class lorestr_i2_report(Report):
    TASKNAME = 'lorestr_i2'
    RUNNING = True
    def __init__(self,xmlnode=None,jobInfo={},jobStatus=None,**kw):
        Report. __init__(self,xmlnode=xmlnode,jobInfo=jobInfo,jobStatus=jobStatus,**kw)
        self.outputXml = jobStatus is not None and jobStatus.lower().count('running')
        if jobStatus is not None and jobStatus.lower() == 'nooutput':
            return

        self.molProbity = False

        self.addDiv(style='clear:both;')

        if jobStatus.lower().count('running'):
            self.addRunningProgress(self)
        else:
            #Finished job
            self.addDefaultReport(self)

    def addRunningProgress(self, parent=None):
        if parent is None: parent=self
        self.addDefaultReport(self)



    def addDefaultReport(self, parent=None):
        if parent is None: parent=self

# Decided not to display parameters
#         if len(self.xmlnode.findall("Parameters")) > 0:
#             parametersFold = parent.addFold(label="LORESTR Parameters", initiallyOpen=True)
#             # parametersFold.addPre(text = 'Received parameters:' )
#
#             table = parametersFold.addTable(select="Parameters[last()]", transpose = True ) # ,style="width:250px;")
#
#             nh = self.xmlnode.findall('Parameters/nh')
#             if len(nh) > 0:
#               table.addData(title="Maximal number of homologues", data=[nh[0].text])
#             nc = self.xmlnode.findall('Parameters/nc')
#             if len(nc) > 0:
#               table.addData(title="Maximal number of chains", data=[nc[0].text])
#             dna = self.xmlnode.findall('Parameters/dna')
#             if len(dna) > 0:
#               table.addData(title="Use DNA/RNA restraints", data=[dna[0].text])
#             mr = self.xmlnode.findall('Parameters/mr')
#             if len(mr) > 0:
#               table.addData(title="Run 100-200 cycles of jelly body after MR", data=[mr[0].text])
#             molprob = self.xmlnode.findall('Parameters/molprobity')
#             if len(molprob) > 0:
#               table.addData(title="Use Molprobity to evaluate geometrical quality", data=[molprob[0].text])
#               if molprob[0].text == 'True':
#                 self.molProbity = True
#             auto = self.xmlnode.findall('Parameters/auto')
#             if len(auto) > 0:
#               table.addData(title="Automatically download homologues from the PDB", data=[auto[0].text])
#               minres = self.xmlnode.findall('Parameters/minRes')
#               if len(minres) > 0:
#                 table.addData(title="Lowest resolution for auto download of homologues (A)", data=[minres[0].text])
#             else:
#               table.addData(title="Automatically download homologues fromt the PDB", data=['False'])
#             cpu = self.xmlnode.findall('Parameters/nCPU')
#             if len(cpu) > 0:
#               table.addData(title="Number of CPUs to use", data=[cpu[0].text])
# #            newFold.addPre(text = 'Auto: ' + self.xmlnode.findall("Parameters/auto")[0].text)
# #            newFold.addPre(text = 'Molprobity: ' + self.xmlnode.findall("Parameters/molprobity")[0].text)
# End of Parameters section

# Global parameter - availability of Molprobity - determines additional columns in the tables
        molprob = self.xmlnode.findall('Parameters/molprobity')
        if len(molprob) > 0:
          if molprob[0].text == 'True':
            self.molProbity = True

# Displaying current status
        if len(self.xmlnode.findall("Protocols/BestProtocol")) < 1: # not finished yet
            parametersFold = parent.addFold(label="LORESTR Status", initiallyOpen=True)
            table = parametersFold.addTable(transpose = True ) # ,style="width:250px;")
            currentStatus = 'Starting...'
            if len(self.xmlnode.findall("AutoDownload")) > 0:
                currentStatus = 'Downloading homologues for restraint generation...'
            if len(self.xmlnode.findall("StartingStructure")) > 0:
                currentStatus = 'Analysing input data...'
            if len(self.xmlnode.findall("Protocols")) > 0:
                currentStatus = 'Running refinement protocols...'
            table.addData(title="Current status:", data=[currentStatus])
# End of Status

# Auto Download
        if len(self.xmlnode.findall("AutoDownload")) > 0:
            if len(self.xmlnode.findall("StartingStructure")) > 0:
                autoDownloadFold = parent.addFold(label="Automatically downloaded homologues for restraint generation", initiallyOpen=False)
            else:
                autoDownloadFold = parent.addFold(label="Automatically downloaded homologues for restraint generation", initiallyOpen=True)
            # self.printXML('AutoDownload', 'failed', 'True')
            if len(self.xmlnode.findall("AutoDownload/failed")) < 1: # Not failed
                table = autoDownloadFold.addTable(select="AutoDownload[last()]" ) # ,style="width:250px;")
                table.addData(title="PDB code", select="pdb")
                table.addData(title="Target chain", select="targetChain")
                table.addData(title="% of residues aligned", select="alignmentLength")
                table.addData(title="Sequence Identity (%)", select="identity")
            else: # Failed BLAST (most likely internet connection problems)
                table = autoDownloadFold.addTable(transpose=True) # ,style="width:250px;")
                table.addData(title="BLAST search failed", data=['Unexpected BLAST output for this chain; trying to continue. Could be internet connection problem - please check your internet connection'])


# End of Auto Download

# Analysing starting structure
        if len(self.xmlnode.findall("StartingStructure")) > 0:
            if len(self.xmlnode.findall("Protocols")) > 0:
                startingStructureFold = parent.addFold(label="Starting Parameters Of The Structure", initiallyOpen=False)
            else:
                startingStructureFold = parent.addFold(label="Starting Parameters Of The Structure", initiallyOpen=True)
            table = startingStructureFold.addTable(select="StartingStructure[last()]", transpose = True ) # ,style="width:250px;")

            for title, select in  [[ "Twinning" ,"StartingStructure/twin" ],
                                  [ "Jelly-body run after MR" ,"StartingStructure/jellyAfterMR" ],
                                  [ "Scaling method"  , "StartingStructure/scaling" ],
                                  [ "Solvent parameters (VDW, ION, RSHR)" , "StartingStructure/solvent" ],
                                  ["Starting R-factor","StartingStructure/Rfact"],
                                  ["Starting R-free","StartingStructure/Rfree"]]:
              node = self.xmlnode.findall(select)
              if len(node) > 0:
                table.addData(title=title, data=[node[0].text])
              else:
                table.addData(title=title, data=['Calculating...'])

            if self.molProbity:
              for title, select in  [[ "Ramachandran outliers (%)" ,"StartingStructure/ramaOut" ],
                                    [ "Ramachandran favourite (%)"  , "StartingStructure/ramaFav" ],
                                    [ "ClashScore percentile" , "StartingStructure/clashPercentile" ],
                                    ["MolProbity score percentile","StartingStructure/molprobPercentile"]]:
                node = self.xmlnode.findall(select)
                if len(node) > 0:
                  table.addData(title=title, data=[node[0].text])
                else:
                  table.addData(title=title, data=['Calculating...'])

            dna = self.xmlnode.findall('StartingStructure/dna')
            if len(dna) > 0:
              table.addData(title="DNA/RNA restraints", select="dna")

# End Starting Structure

# Protocols description
        if len(self.xmlnode.findall("Protocols")) > 0:

          protocolsDescriptionFold = parent.addFold(label="Description of Refinement Protocols", initiallyOpen=False)

          protocols = self.xmlnode.findall('Protocols/*')
          if len(protocols)> 0:
              name = []
              desc = []
              for protocolNode in protocols:
                  try: name.append(protocolNode.findall('Name')[0].text)
                  except: name.append('')

                  # LORESTR separates the parts of a description with "<br>",
                  # which a table cell shows as written.
                  try: desc.append(protocolNode.findall('Description')[0].text.replace('<br>', '; '))
                  except: desc.append('')

              table = protocolsDescriptionFold.addTable(select="Protocols[last()]/*" ) # ,style="width:250px;")
              table.addData(title="Protocol", data=name)
              table.addData(title="Description", data=desc)
# End protocol description

# Protocols execution
          if len(self.xmlnode.findall("Protocols/BestProtocol")) > 0:
              protocolsFold = parent.addFold(label="Execution of Refinement Protocols", initiallyOpen=False)
          else:
              protocolsFold = parent.addFold(label="Execution of Refinement Protocols", initiallyOpen=True)
          table = protocolsFold.addTable(select="Protocols[last()]/*" ) # ,style="width:250px;")
          table.addData(title="Protocol", select="Name")
          table.addData(title="Status", select="Status")
          table.addData(title="R-factor", select="Rfact", expr="x if float(x)>0.0 else '-'")
          table.addData(title="R-free", select="Rfree", expr="x if float(x)>0.0 else '-'")

          if self.molProbity:
            table.addData(title="Ramachandran outliers (%)", select="ramaOut", expr="x if float(x)>0.0 else '-'")
            table.addData(title="Ramachandran favourite (%)", select="ramaFav", expr="x if float(x)>0.0 else '-'")
            table.addData(title="ClashScore percentile", select="clashPercentile", expr="x if float(x)>0.0 else '-'")
            table.addData(title="Molprobity score percentile", select="molprobPercentile", expr="x if float(x)>0.0 else '-'")
            table.addData(title="Q-factor", select="qfact", expr="x if float(x)>0.0 else '-'")
# End Protocols execution



        if len(self.xmlnode.findall("Protocols/BestProtocol")) > 0:

            protocolBestFold = parent.addFold(label="Best Refinement Protocol", initiallyOpen=True)

# Q-score graph
            if self.molProbity:
                COLOUR = (
                          ('ice blue',     '#9cb0fe', '156,176,254,255'),
                          ('gold',         '#b3b13d', '179,177,61,255'),
                          ('coral',        '#ff7f50', '255,127,80,255'),
                          ('grey',         '#808080', '128,128,128,255'),
                          ('pink',         '#ff92ff', '255,146,255,255'),
                          ('sea green',    '#7fbbb5', '127,187,181,255'),
                          ('pale brown',   '#a97d5e', '169,125,94,255'),
                          ('lilac',        '#ae87b9', '174,135,185,255'),
                          ('lemon',        '#ffff80', '255,255,128,255'),
                          ('lawn green',   '#459c4f', '69,156,79,255'),
                          ('pale crimson', '#d23c3e', '210,60,62,255'),
                          ('light blue',   '#419ae1', '65,154,225,255'),
                          ('tan',          '#780000', '120,0,0,255'),
                          ('light green',  '#9aff9a', '154,255,154,255'),
                          ('yellow',       '#ffff00', '255,255,0,255'),
                          ('white',        '#ffffff', '255,255,255,255'),
                          ('blue',         '#0000ff', '0,0,255,255'),
                          ('red',          '#ff0000', '255,0,0,255'),
                          ('green',        '#00ff00', '0,255,0,255'),
                          ('magenta',      '#ff00ff', '255,0,255,255'),
                          ('cyan',         '#00ffe1', '0,255,225,255'),
                          ('purple',       '#9400ff', '148,0,255,255'),
                          ('dark purple',  '#922057', '146,32,87,255'),
                          ('dark cyan',    '#1095a6', '16,149,166,255'),
                          ('black',        '#000000', '0,0,0,255'),
                        )

                graph = protocolBestFold.addFlotGraph( title="Statistics for Protocols", select=".//Protocols",style="height:300px; width:450px; float:left; border:0px;")
                p = graph.addPlotObject()
                p.append('title','Rfree vs. MolProbity percentile')
                p.append('plottype','xy')
                p.append('xlabel','Rfree')
                p.append('ylabel','MolProbity percentile')

                protocols = self.xmlnode.findall('Protocols/*')
                rFrees = []
                for protocolNode in protocols:
                    try: rFrees.append(float(protocolNode.findall('Rfree')[0].text))
                    except: pass
                maxRfree = max(rFrees)
                minRfree = min(rFrees)
                maxX = maxRfree  + ((maxRfree-minRfree) * 0.3)

                counter = 1
                while len(self.xmlnode.findall('.//Protocols/P%d' % counter)) > 0:
                    # The legend shows these titles: name the protocol, or the
                    # four points cannot be told apart.
                    graph.addData (title="Rfree",  select="P%d/Rfree" %counter)
                    graph.addData (title="Protocol %d" % counter,  select="P%d/molprobPercentile" %counter)
                    l = p.append('plotline',xcol=2* counter - 1,ycol=2* counter)
                    l.append('label','%d' % counter)
                    l.append('colour', COLOUR[counter][1])
                    p.append('xrange', rightaxis='false', max='%0.2f' % maxX  ) # min=
                    counter += 1
# End graph

            bestProtocolNumber = int(self.xmlnode.findall("Protocols/BestProtocol")[0].text)

            table = protocolBestFold.addTable(select="Protocols[last()]/*") # ,style="width:250px;")
            for title, select1 in  [["Name", "Protocols/P%d/Name" % bestProtocolNumber],
                                  ["Description", "Protocols/P%d/Description" % bestProtocolNumber]]:
              node1 = self.xmlnode.findall(select1)
              table.addData(title=title, data=[(node1[0].text or '').replace('<br>', '; ')])


            table = protocolBestFold.addTable(select="Protocols[last()]/*" ) # ,style="width:250px;")
            table.addData(title='Structure', data=['Before Refinement', 'After Refinement'])
            for title, select1, select2 in  [["R-factor","StartingStructure/Rfact", "Protocols/P%d/Rfact" % bestProtocolNumber],
                                  ["R-free","StartingStructure/Rfree", "Protocols/P%d/Rfree" % bestProtocolNumber]]:
              node1 = self.xmlnode.findall(select1)
              node2 = self.xmlnode.findall(select2)
              table.addData(title=title, data=[node1[0].text, node2[0].text])

            # LORESTR says in its log when no protocol beat the start; say it here.
            try:
                before = float(self.xmlnode.findall("StartingStructure/Rfree")[0].text)
                after = float(self.xmlnode.findall("Protocols/P%d/Rfree" % bestProtocolNumber)[0].text)
                if after >= before:
                    protocolBestFold.addText(
                        text="No protocol improved R-free (%.3f before, %.3f for the best): "
                             "LORESTR could not improve this structure. " % (before, after),
                        style="color:orange;")
            except (IndexError, TypeError, ValueError):
                pass

            if self.molProbity:
              for title, select1, select2 in  [[ "Ramachandran outliers (%)" ,"StartingStructure/ramaOut",  "Protocols/P%d/ramaOut" % bestProtocolNumber ],
                                    [ "Ramachandran favourite (%)"  , "StartingStructure/ramaFav", "Protocols/P%d/ramaFav" % bestProtocolNumber ],
                                    [ "ClashScore percentile" , "StartingStructure/clashPercentile", "Protocols/P%d/clashPercentile" % bestProtocolNumber ],
                                    ["MolProbity score percentile","StartingStructure/molprobPercentile", "Protocols/P%d/molprobPercentile" % bestProtocolNumber]]:
                node1 = self.xmlnode.findall(select1)
                node2 = self.xmlnode.findall(select2)
                table.addData(title=title, data=[node1[0].text, node2[0].text])



            try:
               validateReport = None
               validateReportNode = self.xmlnode.findall("Validation")[0]
               if validateReportNode is not None:
                  validateReport = validate_protein_report.validate_protein_report(xmlnode=validateReportNode, jobStatus='nooutput', jobInfo=self.jobInfo)

               if validateReport is not None:
                  # The sections the validation task writes now, drawn as
                  # prosmart_refmac draws them.
                  for name, label, draw in (
                          ("B_factors", "B-factor analysis", validateReport.add_b_factors),
                          ("Ramachandran", "Ramachandran plots", validateReport.add_ramachandran),
                          ("Molprobity", "MolProbity analysis", validateReport.add_molprobity)):
                     try:
                        found = validateReportNode.findall(".//" + name)
                        if found and found[0].text != "":
                           draw(parent=self.addFold(label=label, initiallyOpen=False))
                     except Exception:
                        self.addText(text="Warning - " + label + " failed")
            except:
               pass






# End Best Protocol
