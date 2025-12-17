import sys

from ccp4i2.report.CCP4ReportParser import Report


class nautilus_build_refine_report(Report):
    TASKNAME = 'nautilus_build_refine'
    RUNNING = True
    def __init__(self,xmlnode=None,jobInfo={},jobStatus=None,**kw):
        Report. __init__(self, xmlnode=xmlnode, jobInfo=jobInfo, jobStatus=jobStatus, **kw)
        clearingDiv = self.addDiv(style="clear:both;")
        if jobStatus is not None and jobStatus.lower() == 'nooutput':
          return
        elif jobStatus in ['Running','Running remotely']:
          self.summaryTable()
          self.completenessGraph(self, 400, 265)
        else:
          self.finishedText()
          self.table1()
          self.completenessGraph(self, 450, 300)
          clearingDiv = self.addDiv(style="clear:both;")
          self.details(self)
          self.picture(self)


    def summaryTable(self, parent=None):
        if parent is None: parent=self
        xmlPath = './/BuildRefineCycle'
        xmlNodes = self.xmlnode.findall(xmlPath)
#FIXME - XML XPATH LOGIC
        if len(xmlNodes)>0:
            selectString = ".//BuildRefineCycle[1] "
            if len(xmlNodes)>2:
                selectString += " | .//BuildRefineCycle[%d]" % (len(xmlNodes)-1)
            if len(xmlNodes)>1:
                selectString += " | .//BuildRefineCycle[%d]" % len(xmlNodes)
            tableDiv = parent.addDiv(style="height:20em; width:20em; float:left; border:0px;")
            progressTable = tableDiv.addTable(select=selectString)
            progressTable.addData(title="Cycle", select="Number")
            progressTable.addData(title="Number of Nucleotides", select="NautilusResult/Final/NucletidesBuilt")
            progressTable.addData(title="Number of Fragments", select="NautilusResult/Final/FragmentsBuilt")
            progressTable.addData(title="Longest Fragment", select="NautilusResult/Final/NucletidesLongestFragment")
            progressTable.addData(title="R<sub>Work</sub>", select="RefmacResult/r_factor")
            progressTable.addData(title="R<sub>Free</sub>", select="RefmacResult/r_free",   expr="x if float(x)>=0.0 else '-'")

    def completenessGraph(self, parent=None, graph_width=450, graph_height=300):
        graph = parent.addFlotGraph( title="Progress by build-refine iteration", select=".//BuildRefineCycle",style="height:%dpx; width:%dpx; float:right; border:0px;" % (graph_height, graph_width) )
        graph.addData (title="Cycle",  select="Number" )
        graph.addData (title="Number of Nucleotides", select="NautilusResult/Final/NucletidesBuilt")

        graph.addData (title="R<sub>Work</sub>", select="RefmacResult/r_factor")
        graph.addData (title="R<sub>Free</sub>",  select="RefmacResult/r_free" )

        graph.addData (title="RMS<sub>Bonds</sub>", select="RefmacResult/rmsBONDx100", expr="x/100.0" )
        graph.addData (title="RMS<sub>Angles</sub>",  select="RefmacResult/rmsANGLE" ) 

        p = graph.addPlotObject()

        p.append('title', 'Nucleotides built')
        p.append('plottype','xy')
        p.append('xintegral','true')
        p.append('xlabel','Cycle')
        l = p.append('plotline',xcol=1,ycol=2)

        p = graph.addPlotObject()
        p.append('title','R-factors after each iteration')
        p.append('plottype','xy')
        p.append('xintegral','true')
        p.append('xlabel','Cycle')
        l = p.append('plotline',xcol=1,ycol=4)
        l.append('label','R<sub>Free</sub>')
        l.append('colour','gold')
        l = p.append('plotline',xcol=1,ycol=3)
        l.append('label','R<sub>Work</sub>')
        l.append('colour','lightblue')
        if len(p.validate())>0: print(p.validate())

        p = graph.addPlotObject()
        p.append('title','Geometry after each iteration')
        #p.append('yrange', rightaxis='true', min='auto',max='auto')
        p.append('plottype','xy')
        p.append('xintegral','true')
        p.append('xlabel','Cycle')
        l = p.append('plotline',xcol=1,ycol=5, rightaxis='true')
        l.append('label','RMS<sub>Bonds</sub>')
        l.append('colour','firebrick')
        l = p.append('plotline',xcol=1,ycol=6)
        l.append('label','RMS<sub>Angles</sub>')
        l.append('colour','green')
        if len(p.validate())>0: print(p.validate())

    def details(self,parent=None):
        if parent is None:
            parent = self
        fold = parent.addFold(label="Detailed progress by iteration")
        table = fold.addTable(select=".//BuildRefineCycle", transpose=True, downloadable=True,id='details')
        
        table.addData ( title='Iteration', select='Number', expr='int(x)' )

        for title,select in  [[ "Nucleotides built",   "NautilusResult/Final/NucletidesBuilt" ],
                              [ "Longest fragment",    "NautilusResult/Final/NucletidesLongestFragment" ],
                              [ "Number of fragments", "NautilusResult/Final/FragmentsBuilt" ] ]:
            table.addData(title=title,select=select)
        for title,select,expr in [[ "R<sub>Work</sub>",     "RefmacResult/r_factor", "x"],
                                  [ "R<sub>Free</sub>",     "RefmacResult/r_free","x if float(x)>0.0 else '-' " ],
                                  [ "RMS<sub>Bonds</sub>" , "RefmacResult/rmsBONDx100", "round(x/100,3)"],
                                  [ "RMS<sub>Angles</sub>", "RefmacResult/rmsANGLE","x" ] ]:
            table.addData(title=title,select=select,expr=expr)

    def table1(self,parent=None):
        if parent is None:
            parent = self
        tableDiv = parent.addDiv(style="height:30em;width:20em;float:left;border:0px;")
        table = tableDiv.addTable(select=".", transpose=True, id='table_1') 
        try:
          for title,select in  [[ "Nucleotides built" ,  "BuildRefineCycle[last()]/NautilusResult/Final/NucletidesBuilt" ],
                                [ "Longest fragment" ,   "BuildRefineCycle[last()]/NautilusResult/Final/NucletidesLongestFragment" ],
                                [ "Number of fragments", "BuildRefineCycle[last()]/NautilusResult/Final/FragmentsBuilt" ] ]:
              table.addData(title=title,select=select)
          for title,select,expr in [[ "R<sub>Work</sub>"  ,   "BuildRefineCycle[last()]/RefmacResult/r_factor", "x"],
                                    [ "R<sub>Free</sub>" ,    "BuildRefineCycle[last()]/RefmacResult/r_free","x if float(x)>0.0 else '-' " ],
                                    [ "RMS<sub>Bonds</sub>" , "BuildRefineCycle[last()]/RefmacResult/rmsBONDx100", "round(x/100,3)"],
                                    [ "RMS<sub>Angles</sub>", "BuildRefineCycle[last()]/RefmacResult/rmsANGLE","x" ] ]:
              table.addData(title=title,select=select,expr=expr)
        except Exception as e:
            print("ERROR nautilus_build_refine_report ",e, file=sys.stderr)

    def finishedText(self,parent=None) :
        if parent is None: parent = self
        try:
            frgb = float(self.xmlnode.findall('BuildRefineCycle/NautilusResult/Final/FragmentsBuilt')[-1].text)
            resb = float(self.xmlnode.findall('BuildRefineCycle/NautilusResult/Final/NucletidesBuilt')[-1].text)
            parent.append( "<p>%d nucleotides were built in %d fragments.</p>"%(resb,frgb) )
            #ress = float(self.xmlnode.findall('BuildRefineCycle/NautilusResult/Final/NucletidesSequenced')[-1].text)
            #parent.append( "<p>%d nucleotides were built in %d fragments. Of these, %d nucleotides were assigned to the sequence.</p>"%(resb,frgb,ress) )
        except Exception as e:
            parent.append( "<p>Model building results are not yet available.</p>" )
            print("ERROR nautilus_build_refine_report ",e, file=sys.stderr)
        try:
            rwrk = float(self.xmlnode.findall('BuildRefineCycle/RefmacResult/r_factor')[-1].text)
            rfre = float(self.xmlnode.findall('BuildRefineCycle/RefmacResult/r_free')[-1].text)
            rbnd = float(self.xmlnode.findall('BuildRefineCycle/RefmacResult/rmsBONDx100')[-1].text)*0.01
            rang = float(self.xmlnode.findall('BuildRefineCycle/RefmacResult/rmsANGLE')[-1].text)
            if rwrk > 0.5:    s = "the model is very incomplete or wrong"
            elif rwrk > 0.4:  s = "the model is substantially incomplete and may contain incorrect regions"
            elif rwrk > 0.35: s = "the model is likely to contain correct regions but requires further work"
            else:             s = "the model is approaching completion"
            parent.append( "<p>The refinement R-factor is %5.2f, and the free-R factor is %5.2f. The RMS bond deviation is %5.3f A.<br/>On the basis of the refinement statistics, %s.</p>"%(rwrk,rfre,rbnd,s) )
        except Exception as e:
            parent.append( "<p>Refinement results are not yet available.</p>" )
            print("ERROR nautilus_build_refine_report ",e, file=sys.stderr)

    def picture(self,parent=None) :
        pic = parent.addPicture(label="Autobuilt structure",sceneFile="$CCP4I2/pipelines/nautilus_build_refine/script/nautilus_1.scene.xml",id='autobuild_1')


# Temporary hard-wire to save typing..
#from ccp4i2.core import CCP4Utils
#r = nautilus_build_refine_report(xmlnode=CCP4Utils.openFileToEtree('/home/cowtan/CCP4I2_PROJECTS/test/CCP4_JOBS/job_30/program.xml'))
#r.as_html_file('/tmp/report.html')
