from ccp4i2.report import Report


class molrep_den_report(Report):
  TASKNAME = 'molrep_den'
  def __init__(self,xmlnode=None,jobInfo={},**kw):
    Report. __init__(self,xmlnode=xmlnode,jobInfo=jobInfo,**kw)
    results = self.addResults()
    # Relative to the program.xml root (MolrepResult): an absolute path on an
    # element raises, and failed every molrep report (#786). A FlotGraph, as
    # every other report uses: the base Graph writes nothing out, so this graph
    # had never been seen. The element is TF_sig; "Tf_sig" selected nothing.
    graph = results.addFlotGraph( title="Best TF peak vs RF peak No" , select="RFpeaks/RFpeak" )
    graph.addData (title="RF_peak_No"  , select="RF" )
    graph.addData (title="Score" , select="Score" )
    graph.addData (title="TF_sig" , select="TF_sig" )
    graph.addPlot ( plot = '''<plot>
<title>Best TF peak score vs RF peak No</title>
<plottype>xy</plottype>
<plotline xcol="1" ycol="2">
<linestyle>.</linestyle>
<markeredgewidth>0</markeredgewidth>
<colour>blue</colour>
</plotline>
</plot> ''' )
    graph.addPlot ( plot = '''<plot>
<title>TF/sig(TF) vs RF peak No</title>
<plottype>xy</plottype>
<plotline xcol="1" ycol="3">
<linestyle>.</linestyle>
<markeredgewidth>0</markeredgewidth>
<colour>red</colour>
</plotline>
</plot> ''' )

    pic = results.addPicture(label="Final structure",sceneFile="$CCP4I2/wrappers/molrep_den/script/molrep_den_1.scene.xml")

    fold = results.addFold(label='Show details')
    table = fold.addTable(select=".//RFpeaks/RFpeak")
    for title,select in  [[ "RF" ,"RF" ],
                          [ "TF"  , "TF" ],
                          [ "TF_sig" , "TF_sig" ],
                          [ "TFcntrst" ,"TFcntrst" ],
                          [ "PFind" , "PFind" ],
                          [ "PF" , "PF"],
                          [ "PFmod" , "PFmod" ],
                          [ "wRfac" , "wRfac" ],
                          [ "Score" , "Score" ],
                          [ "Cntrst" , "Cntrst" ],
                          [ "For" , "for"  ] ] :
        
        table.addData(title=title,select=select)
