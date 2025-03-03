from report.CCP4ReportParser import *
import sys

class cmapcoeff_report(Report):
    # Specify which gui task and/or pluginscript this applies to
    TASKNAME = 'cmapcoeff'
    
    def __init__(self,xmlnode=None,jobInfo={},**kw):
        Report.__init__(self,xmlnode=xmlnode,jobInfo=jobInfo,**kw)

        #if jobStatus is None or jobStatus.lower() is 'nooutput': return
     
        import os
        
        results = self.addResults()
        results.append ( 'Please find below the output files. If you want to do a peak search, you can select <i>Manual model rebuilding</i>. ' ) 
#        results.append ( '<h3>Map statistics:</h3>' )
#        results.append ( '&#160; &#160; &#160; &#160; Number of points: ' + self.xmlnode.xpath ( "//Cfft/NPoints" )[0].text )
#        results.append ( '&#160; &#160; &#160; &#160; 1st moment about zero (mean) =' + self.xmlnode.xpath ( "//Cfft/FirstMomZero" )[0].text +
#                ' &#160; &#160; &#160; &#160; 1st moment about mean =' + self.xmlnode.xpath ( "//Cfft/FirstMomMean" )[0].text )
#        results.append ( '&#160; &#160; &#160; &#160; 2nd moment about zero =' + self.xmlnode.xpath ( "//Cfft/SecondMomZero" )[0].text +
#                ' &#160; &#160; &#160; &#160; 2nd moment about mean =' + self.xmlnode.xpath ( "//Cfft/SecondMomZero" )[0].text )
#        results.append ( '&#160; &#160; &#160; &#160; 3rd moment about zero =' + self.xmlnode.xpath ( "//Cfft/ThirdMomZero" )[0].text +
#                ' &#160; &#160; &#160; &#160; 3rd moment about mean =' + self.xmlnode.xpath ( "//Cfft/ThirdMomZero" )[0].text )
#        results.append ( '&#160; &#160; &#160; &#160; Range: min =' + self.xmlnode.xpath ( "//Cfft/Min" )[0].text +
#                ' &#160; &#160; &#160; &#160; max =' + self.xmlnode.xpath ( "//Cfft/Max" )[0].text )

if __name__ == "__main__":
    import sys
    cmapcoeff_report(xmlFile=sys.argv[1],jobId=sys.argv[2])

