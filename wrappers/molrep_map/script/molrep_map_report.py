from report.CCP4ReportParser import *
import sys

class molrep_map_report(Report):
    # Specify which gui task and/or pluginscript this applies to
    TASKNAME = 'molrep_map'
    RUNNING = False
    def __init__(self,xmlnode=None,jobInfo={},jobStatus=None,**kw):
        Report. __init__(self,xmlnode=xmlnode,jobInfo=jobInfo,**kw)
        results = self.addResults()
        for hand in ['Original', 'Flipped']:
            results.addText(text=hand)
            table = results.addTable(select=f"./{hand}/MolrepResult/RFpeaks/RFpeak")
            for title,select in  [[ "RF" ,"RF" ],
                                [ "TF"  , "TF" ],
                                [ "Tf_sig" , "Tf_sig" ],
                                [ "TFcntrst" ,"TFcntrst" ],
                                [ "PFind" , "PFind" ],
                                [ "PF" , "PF"],
                                [ "PFmin" , "PFmin" ],
                                [ "wRfac" , "wRfac" ],
                                [ "Score" , "Score" ],
                                [ "Cntrst" , "Cntrst" ],
                                [ "For" , "for"  ] ] :
                
                table.addData(title=title,select=select)

        for hand in ['Original', 'Flipped']:
            fold = results.addFold(label=f'{hand} hand graphs')
            graph = fold.addGraph( title="Best TF peak vs RF peak No" , select=f"./{hand}/MolrepResult/RFpeaks/RFpeak" )
            graph.addData (title="RF_peak_No"  , select="RF" )
            graph.addData (title="Score" , select="Score" )
            graph.addData (title="Tf_sig" , select="Tf_sig" )
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


        self.addTaskReferences()
    
if __name__ == "__main__":
    import sys
    molrep_mr_report(xmlFile=sys.argv[1],jobId=sys.argv[2])       


