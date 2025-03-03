from __future__ import print_function

from lxml import etree
class logScraper(object):
    
    def __init__(self, *args, **kws):
        super(logScraper,self).__init__()
        
        self.xmlroot = kws.get('xmlroot',etree.Element('REFMAC'))
        self.flushXML = kws.get('flushXML',self._flushXML)
        
        self.blockHandlers = [{'triggerText':'Twin operators with Rmerge', 'finishFunction':lambda line: '---' in line,'parseFunction':self.parseTwin1},
                              {'triggerText':'Twin domains with fraction', 'finishFunction':lambda line: '---' in line,'parseFunction':self.parseTwin2},
                              {'triggerText':'Bond distance outliers', 'finishFunction':lambda line: len(line.split()) == 0 and len(self.xmlroot.xpath('//OutliersByCriteria/Bond/Outlier')) != 0,'parseFunction':self.parseBondOutlier},
                              {'triggerText':'VDW outliers', 'finishFunction':lambda line: len(line.split()) == 0 and len(self.xmlroot.xpath('//OutliersByCriteria/VDW/Outlier')) != 0,'parseFunction':self.parseVDWOutlier},
                              {'triggerText':'Alignment results', 'finishFunction':lambda line: '---' not in line and ':' not in line and '*' not in line,'parseFunction':self.parseAlignment},
                              {'triggerText':'Status     Link', 'finishFunction':lambda line: not line.startswith("Unused") and not line.startswith("Existing"),'parseFunction':self.parseLinks}
                              ]
        self.triggers = [{'triggerText':'CGMAT cycle number =','tagPath':'Cycle/number', 'parentTag':'//REFMAC', 'extractor':lambda x: len(self.xmlroot.xpath('Cycle/number'))-1,'doFlush':True},
                         {'triggerText':'***TLS refinement cycle***','tagPath':'Cycle/number', 'parentTag':'//REFMAC', 'extractor':lambda x: len(self.xmlroot.xpath('Cycle/number'))-1,'doFlush':True},
                         {'triggerText':'Rigid body cycle =','tagPath':'Cycle/number', 'parentTag':'//REFMAC', 'extractor':lambda x: len(self.xmlroot.xpath('Cycle/number'))-1,'doFlush':True},
                         {'triggerText':'WRITTEN OUTPUT MTZ FILE','tagPath':'Cycle/number', 'parentTag':'//REFMAC', 'extractor':lambda x: len(self.xmlroot.xpath('Cycle/number'))-1,'doFlush':True},
                         {'triggerText':'***TLS refinement cycle***','tagPath':'TLSMode', 'parentTag':'//REFMAC', 'doFlush':True},
                         {'triggerText':'Rigid body cycle =','tagPath':'RigidMode', 'parentTag':'//REFMAC', 'doFlush':True},
                         {'triggerText':'Weight matrix','tagPath':'WeightUsed', 'parentTag':'//REFMAC/Cycle', 'extractor':lambda x: x.split()[2],'doFlush':True},
                         {'triggerText':'Overall R factor','tagPath':'r_factor', 'parentTag':'//REFMAC/Cycle', 'extractor':lambda x: x.split()[4],'doFlush':True},
                         {'triggerText':'Free R factor','tagPath':'r_free', 'parentTag':'//REFMAC/Cycle', 'extractor':lambda x: x.split()[4],'doFlush':True},
                         {'triggerText':'Bond distances: refined atoms','tagPath':'rmsBonds', 'parentTag':'//REFMAC/Cycle', 'extractor':lambda x: x.split()[5],'doFlush':True},
                         {'triggerText':'Bond distances: refined atoms','tagPath':'rmsBondsx10', 'parentTag':'//REFMAC/Cycle', 'extractor':lambda x: 10.*float(x.split()[5]),'doFlush':True},
                         ]
        self._currentBlockHandler = None

    def scrapeTest(self):
        logFile = getLog()
        self.scrapeBuffer(logFile)
        
    def scrapeBuffer(self, logFile):
        lines = logFile.split("\n")
        for line in lines:
            self.processLine(line)
        #print etree.tostring(self.xmlroot,pretty_print=True)

    def scrapeFile(self, fileName):
        self._oldFlushXML = self.flushXML
        self.flushXML = self._noOp
        with open(fileName,'r') as logFile:
            chunk = logFile.read()
            self.processLogChunk(chunk)
        self.flushXML = self._oldFlushXML
        self.flushXML()

    def _noOp(self):
        pass

    def processLogChunk(self, logChunk):
        if not hasattr(self, "_logTail"): self._logTail = ""
        textToProcess = self._logTail + logChunk
        
        if  textToProcess.endswith("\n"): endsWithSlashN = True
        else:  endsWithSlashN = False

        lines = textToProcess.split("\n")
        nLines = len(lines)
        if not endsWithSlashN:
            self._logTail = lines[-1]
            nLines -= 1
        
        for iLine in range(nLines):
            self.processLine(lines[iLine])

    def processLine(self, line):
        if self._currentBlockHandler is not None:
            if self._currentBlockHandler['finishFunction'](line):
                #print 'Recognised closing line',line
                self._currentBlockHandler = None
            else:
                #print 'Parsing line',line
                self._currentBlockHandler['parseFunction'](line)
            return
                
        for blockHandler in self.blockHandlers:
            if blockHandler['triggerText'] in line:
                self._currentBlockHandler = blockHandler
                #print 'Using blockHandler', blockHandler
        
        for trigger in self.triggers:
            if trigger['triggerText'] in line:
                try:
                    parentNode = self.xmlroot.xpath(trigger['parentTag'])[-1]
                    for tag in trigger['tagPath'].split("/"):
                        parentNode = etree.SubElement(parentNode,tag)
                    parentNode.text = str(trigger['extractor'](line))
                except:
                    print('parent ',trigger['parentTag'],' of trigger ',trigger['triggerText'],' not found')
                if trigger['doFlush']: self.flushXML()

    def _flushXML(self):
        print(etree.tostring(self.xmlroot, pretty_print=True))

    def parseTwin1(self, line):
        tokens = line.split(':')
        if len(tokens) > 1:
            try: parentNode = self.xmlroot.xpath('Twinning')[0]
            except: parentNode = etree.SubElement(self.xmlroot,'Twinning')
            twinningBySymmetryOperatorNode = etree.SubElement(parentNode,'TwinningSymmetryOperator')
            operatorTokens = tokens[0].split()
            if len(operatorTokens) > 3:
                operatorText = operatorTokens[-3]+operatorTokens[-2]+operatorTokens[-1]
                twinningSymmetryOperatorNode = etree.SubElement(twinningBySymmetryOperatorNode,'SymmetryOperator')
                twinningSymmetryOperatorNode.text = operatorText
            otTokens = tokens[1].split('=')
            if len(otTokens) > 1:
                operatorText = otTokens[1]
                if len(operatorTokens) > 1:
                    twinningRmergeNode = etree.SubElement(twinningBySymmetryOperatorNode,'Rmerge')
                    twinningRmergeNode.text = operatorText

    def parseTwin2(self, line):
        tokens = line.split(':')
        if len(tokens) > 3:
            try: parentNode = self.xmlroot.xpath('Twinning')[0]
            except: parentNode = etree.SubElement(self.xmlroot,'Twinning')
            twinningByTwinOperatorNode = etree.SubElement(parentNode,'TwinOperator')
            twinningTwinOperatorNode = etree.SubElement(twinningByTwinOperatorNode,'SymmetryOperator')
            twinningTwinOperatorNode.text = tokens[1]
            try:
                fractionText = tokens[2].split(';')[0].split('=')[1]
                twinningFractionNode = etree.SubElement(twinningByTwinOperatorNode,'Fraction')
                twinningFractionNode.text = fractionText
            except:
                pass

    def parseLinks(self, line):
            try: parentNode = self.xmlroot.xpath('Links')[0]
            except: parentNode = etree.SubElement(self.xmlroot,'Links')
            xml_link = etree.SubElement(parentNode,"Link")
            xml_link.text = str(line.strip())

    def parseAlignment(self, line):
        if 'No of aligned' not in line and '---' not in line and 'Alignment results' not in line:
            try: parentNode = self.xmlroot.xpath('NCS')[0]
            except: parentNode = etree.SubElement(self.xmlroot,'NCS')
            tokens = [token.strip() for token in line.split(':')]
            equivalenceNode = etree.SubElement(parentNode,'equivalence',selection1 = tokens[2],selection2=tokens[3],nEquivalent=tokens[4],score=tokens[5],rms=tokens[6],ave_rmsLoc=tokens[7])

    def parseBondOutlier(self, line):
        parentNode = self.xmlroot
        if "will be monitored" in line:
            for tag in "OutliersByCriteria/Bond".split("/"):
                oldParentNode=parentNode
                try: parentNode = oldParentNode.xpath(tag)[0]
                except: parentNode = etree.SubElement(oldParentNode,tag)
            criterionNode = etree.SubElement(parentNode,"Criteria")
            criterionNode.text = line
    
        elif len(line.split()) > 0:
            for tag in "OutliersByCriteria/Bond".split("/"):
                oldParentNode=parentNode
                try: parentNode = oldParentNode.xpath(tag)[0]
                except: parentNode = etree.SubElement(oldParentNode,tag)
            criterionNode = etree.SubElement(parentNode,"Outlier")
            outlier = {}
            outlier['chainId1'] = line[0:1]
            outlier['resId1'] = line[4:8]
            outlier['at1'] = line[12:16]
            outlier['ins1'] = line[17:18]
            outlier['chainId2'] = line[21:22]
            outlier['resId2'] = line[25:29]
            outlier['at2'] = line[33:37]
            outlier['ins2'] = line[38:39]
            outlier['mod'] = line[45:51]
            outlier['ideal'] = line[56:62]
            outlier['sigma'] = line[80:86]
            for key,value in list(outlier.items()):
                criterionNode.set(key,value)

    def parseVDWOutlier(self, line):
        parentNode = self.xmlroot
        if "will be monitored" in line:
            for tag in "OutliersByCriteria/VDW".split("/"):
                oldParentNode=parentNode
                try: parentNode = oldParentNode.xpath(tag)[0]
                except: parentNode = etree.SubElement(oldParentNode,tag)
            criterionNode = etree.SubElement(parentNode,"Criteria")
            criterionNode.text = line
    
        elif len(line.split()) > 0:
            for tag in "OutliersByCriteria/VDW".split("/"):
                oldParentNode=parentNode
                try: parentNode = oldParentNode.xpath(tag)[0]
                except: parentNode = etree.SubElement(oldParentNode,tag)
            criterionNode = etree.SubElement(parentNode,"Outlier")
            outlier = {}
            outlier['chainId1'] = line[0:1]
            outlier['resId1'] = line[4:8]
            outlier['at1'] = line[12:16]
            outlier['ins1'] = line[17:18]
            outlier['chainId2'] = line[21:22]
            outlier['resId2'] = line[25:29]
            outlier['at2'] = line[33:37]
            outlier['ins2'] = line[38:39]
            outlier['mod'] = line[45:51]
            outlier['ideal'] = line[56:62]
            outlier['sigma'] = line[79:84]
            outlier['operator'] = line[90:102]
            outlier['type'] = line[111:112]
            
            for key,value in list(outlier.items()):
                criterionNode.set(key,value)



