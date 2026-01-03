from ccp4i2.report import Report


class fft_report(Report):
    TASKNAME = 'fft'
    
    def __init__(self,xmlnode=None,jobInfo={},**kw):
        Report.__init__(self,xmlnode=xmlnode,jobInfo=jobInfo,**kw)

        results = self.addResults()
#        results.append ( '<h3>Map statistics:</h3>' )
        results.append ( '&#160; &#160; &#160; &#160; Number of points: ' + self.xmlnode.findall ( ".//Cfft/NPoints" )[0].text )
        results.append ( '&#160; &#160; &#160; &#160; 1st moment about zero (mean) =' + self.xmlnode.findall ( ".//Cfft/FirstMomZero" )[0].text +
                ' &#160; &#160; &#160; &#160; 1st moment about mean =' + self.xmlnode.findall ( ".//Cfft/FirstMomMean" )[0].text )
        results.append ( '&#160; &#160; &#160; &#160; 2nd moment about zero =' + self.xmlnode.findall ( ".//Cfft/SecondMomZero" )[0].text +
                ' &#160; &#160; &#160; &#160; 2nd moment about mean =' + self.xmlnode.findall ( ".//Cfft/SecondMomZero" )[0].text )
        results.append ( '&#160; &#160; &#160; &#160; 3rd moment about zero =' + self.xmlnode.findall ( ".//Cfft/ThirdMomZero" )[0].text +
                ' &#160; &#160; &#160; &#160; 3rd moment about mean =' + self.xmlnode.findall ( ".//Cfft/ThirdMomZero" )[0].text )
        results.append ( '&#160; &#160; &#160; &#160; Range: min =' + self.xmlnode.findall ( ".//Cfft/Min" )[0].text +
                ' &#160; &#160; &#160; &#160; max =' + self.xmlnode.findall ( ".//Cfft/Max" )[0].text )
