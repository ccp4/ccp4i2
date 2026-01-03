# Some routines useful for running and analysis
#   class Tabledata
#   class AnalyseGraph
#   class AnalysisLog
#   class Makexmlgraph
# Outside classes
#   findvalueinlist
#   interpolate
#   addElement

from lxml import etree



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
class Tabledata:
    # A table of data eg score vs. resolution
    def __init__(self, tablexml):
        # initialise from XML, expect CCP4Table, headers, data
        # This should be an Element CCP4Table

        if tablexml.tag != 'CCP4Table':
            print("Tabledata: wrong tag ",tablexml.tag)
            return
        
        self.title = tablexml.attrib['title']
        self.graphid =  tablexml.attrib['id']
        self.headers = tablexml.findall('headers')[0].text
        self.data = tablexml.findall('data')[0].text

        self.collabels = self.headers.split()
        ncols = len(self.collabels)
        columns = [[] for i in range(ncols)]  # n lists

        # Parse data into columns
        lines = self.data.splitlines()
        for line in lines:
            if len(line) > 1:  # ignore blank lines
                fields = line.split()
                for i in range(ncols):
                    columns[i].append(float(fields[i]))

        self.columndata = {}
        for i in range(ncols):
            key = self.collabels[i]
            self.columndata[key] = columns[i]

    # . . . . . . . . . . . . . . .
    def getcollabels(self):
        return self.collabels
    
    # . . . . . . . . . . . . . . .
    def findlimit(self, xcol, ycol, threshold, falling=True):
        # Tabledata.findlimit
        # x values are in self.columndata[xcol]
        # y values are in self.columndata[ycol]
        
        xlist = self.columndata[xcol]
        ylist = self.columndata[ycol]
        allOK, xval = findvalueinlist(xlist, ylist, threshold, falling)
        #print("findlimit",allOK,xval,threshold)
        # returns allOK, xval
        #  allOk = -1 no valid points, +1 all valid, else = 0 
        return allOK, xval

    # . . . . . . . . . . . . . . .
    def meanValue(self, xrange, ycol, xcol):
        # estimate mean y value in "ycol" over a range of x values
        # in column "xcol". Weight contributions to average by the
        # proportion falling in the bin
        xlist = self.columndata[xcol]
        ylist = self.columndata[ycol]

        lastx = xlist[0]
        n = 0
        sd = 0.0
        for x in xlist[1:]:
            n += 1
            sd += x-lastx
            lastx = x
        dx = sd/n   # average interval (assume equal)

        nx = len(xlist)
        inrange = False
        wysum = 0.0
        wsum = 0.0

        x0 = xlist[0] - 0.5*dx  # in case first point isn't at 0
        x1 = (xrange[0] - x0)/dx   # first point
        x2 = (xrange[1] - x0)/dx   # last point
        i1 = int(x1)    # for end points
        i2 = int(x2)
        # requested range may be beyond data range
        if i1 > nx-1:
            return 0.0;
        i1 = max(i1, 0)    # for end points
        i2 = min(i2, nx-1) 
        w1 = 1.0 - (x1-i1)%1
        w2 = (x2-i2)%1

        wsum += w1 + w2
        wysum += w1*ylist[i1] + w2*ylist[i2]
        
        if (i2-1) > (i1+1):
            for i in range(i1+1, i2):
                wsum += 1.0
                wysum += ylist[i]

                meany = 0.0
        if wsum > 0.0:
            meany = wysum/wsum

        return meany

# . . . . . . . . . . . . . . .
#   Routines not in a class
# . . . . . . . . . . . . . . .
def findvalueinlist(xlist, ylist, threshold, falling=True):
    # return the x value corresponding (approximately) to the
    # point where y = threshold
    # If falling True, then look for point where y falls below threshold
    #  else rises above

    # returns allOK, xval
    #  allOk = -1 no valid points, +1 all valid, else = 0 

    nrows = len(ylist)
    if nrows == 0:
        return 0.0
    if nrows == 1:
        return xlist[0]

    # Check minimum and maximum values against the threshold
    miny = min(ylist)
    maxy = max(ylist)

    allOK = 0
    if falling:
        if miny > threshold:
            # all points are OK
            allOK = +1
            xval = xlist[-1]
        elif maxy < threshold:
            # No points are OK
            allOK = -1
            xval = xlist[0]
    else: # rising
        if maxy < threshold:
            # all points are OK
            allOK = +1
            xval = xlist[-1]
        elif miny > threshold:
            # No points are OK
            allOK = -1
            xval = xlist[0]
                
    if allOK == 0:
        # somewhere in the middle
        i2 = None # 1st point found beyond threshold
        for i in range(nrows):
            y = ylist[i]
            if (falling and  y < threshold) or \
                   ((not falling) and  y > threshold):
                i2 = i
                i1 = max(0,i-1)
                break
        if i2 is None: print("Help")
        # Interpolate
        xval = interpolate(threshold, i1,i2, xlist, ylist)

    return allOK, xval 

# . . . . . . . . . . . . . . .
def interpolate(threshold, i1,i2, x, y):
    if i1 == i2:
        # end point
        return x[i1]

    dx = x[i2] - x[i1] 
    r = (threshold - y[i1])/ (y[i2] -y[i1])
    xval = x[i1] + r * dx
    return xval
        
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# . . . . . . . . .
class AnalyseGraph:
    #  Analyse graph from XML, eg for resolution estimate
    def __init__(self, xmlblock):
        # Initialise from XML: select all CCP4Table objects
        # convert to tablegraph objects in dictionary
        ccp4tables = xmlblock.findall('CCP4Table')
        #print("AnalyseGraph ccp4tables ",len(ccp4tables))
        self.tabledata = {}
        for ccp4table in ccp4tables:
            graphid = ccp4table.attrib['id']
            self.tabledata[graphid] = Tabledata(ccp4table)
        
    def getcollabels(self, name):
        if name in self.tabledata:
            return self.tabledata[name].getcollabels()
        return None
            
    def getresolutionlimit(self, name, xcol, ycol, threshold, falling=True):
        #print("AnalyseGraph.getresolutionlimit", name)
        limit = 0.0
        allOK = -1
        if name in self.tabledata:
            allOK, limit = self.tabledata[name].findlimit(
                xcol, ycol, threshold, falling)
        return allOK, limit


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
def addElement(containerXML, elementname, elementtext,
               attributes=None):
    e2 = etree.Element(elementname, attrib=attributes)
    if elementtext is not None: e2.text = elementtext
    containerXML.append(e2)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
class AnalysisLog:
    #  Get stuff from logfile
    #  Oh, joy!
    def __init__(self, logdata, tag=None):
        # XML tag
        if tag is None:
            xtag = 'Analysis'
        else:
            xtag = tag
        self.xmlroot = etree.Element(xtag)
                
        lines = logdata.splitlines()
        rlines = reversed(lines)  # for searching backwards
        self.totalInformation(lines)
        self.resolution(rlines)
        self.tNCS = self.tncs(lines) # False if no tNCS
        self.twinning(lines)
        self.anisotropy(lines)

    # . . . . . . . . . . . . . . . . . . . . . . . . . . .
    def totalInformation(self, lines):
        # Summary of information content

        text = 'DATA INFORMATION CONTENT'
        skip = 3
        keep = 14
        infoblock = self.getBlock(lines, text, skip, keep)
        if len(infoblock) == 0:
            return

        infocontentxml = etree.Element('InformationContent')

        tag = "Information can only be estimated for intensity data"
        if tag in infoblock[0]:
            # No data
            addElement(infocontentxml, 'Message', tag)
        else:
            totals = self.findline(infoblock,
                               'Estimated total information').split()
            totalbits = totals[4]
            nref = totals[10]
            #print("Total bits {} in {} reflections".format(totalbits, nref))
            averagebits = self.findline(infoblock, 'average').split()[2]
            #print("Average {}".format(averagebits))
            
            addElement(infocontentxml, 'Text', "\n".join(infoblock))
            addElement(infocontentxml, 'TotalBits', totalbits)
            addElement(infocontentxml, 'Nreflections', nref)
            addElement(infocontentxml, 'Averagebits', averagebits)

        self.xmlroot.append(infocontentxml)

    # . . . . . . . . . . . . . . . . . . . . . . . . . . .
    def resolution(self, rlines):
        # Resolution, last occurance
        tag = 'Resolution of All Data'
        line =  self.findline(rlines, tag)
        if line is None: return
        
        self.allresolution = line.split()[5]  # Store for access
        nrefall = line.split()[7].strip('()')
        
        tag = 'Resolution of Selected Data'
        line =  self.findline(rlines, tag)
        if line is None: return

        selectedresolution = line.split()[5]
        nrefselected = line.split()[7].strip('()')

        resolutionxml = etree.Element('Resolution')
        addElement(resolutionxml, 'ResolutionAll', self.allresolution)
        addElement(resolutionxml, 'NResolutionAll', nrefall)
        addElement(resolutionxml, 'SelectedResolution', selectedresolution)
        addElement(resolutionxml, 'NSelectedResolution', nrefselected)

        self.xmlroot.append(resolutionxml)
        
    # . . . . . . . . . . . . . . . . . . . . . . . . . . .
    def tncs(self, lines):
        tNCSxml = etree.Element('tNCS')
        tag = 'No tNCS found in Patterson'
        tNCS = False
        s = 'False'
        found = False
        if self.findline(lines, tag) is not None:
            # No tNCS
            found = True
        else:
            # Patterson
            text = 'Translational NCS Analysis Table'
            skip = 4
            keep = 1
            pattstuff = self.getBlock(lines, text, skip, keep)[0].split()
            if len(pattstuff) > 0:
                tNCS = True
                s = 'True'
                found = True
                peaksize = pattstuff[2]
                # Vector and angle
                angle = self.findline(lines, 'Final angle').split()[3:6]
                vector = self.findline(lines, 'Final vector').split()[3:6]
                dvalues = self.findline(lines, 'tNCS D-values').split()[6:]
                # xml
                pattxml = etree.SubElement(tNCSxml, 'NonOriginPatterson')
                addElement(pattxml, 'PeakHeight', peaksize)
                addElement(pattxml, 'Vector', ' '.join(vector))
                addElement(pattxml, 'Angle', ' '.join(angle))
                addElement(tNCSxml, 'DvalueRange', ' '.join(dvalues))

        if found:
            tNCSxml.set('tNCS', s)
            self.xmlroot.append(tNCSxml)

        return tNCS

    # . . . . . . . . . . . . . . . . . . . . . . . . . . .
    def twinning(self, lines):
        # Summary of information content
        text = 'tNCS/Twin Detection Table'
        skip = 2
        keep = 15
        twinblock = self.getBlock(lines, text, skip, keep)

        # Ignore uncorrected errors
        # Theoretical moments including measurement error 
        # List: centric, acentric
        theoretical = self.extractnumbers(twinblock[3], 3, 2)
        theoreticaltwin = self.extractnumbers(twinblock[4], 4, 2)
        # Experimental value, line depends on whether there is tNCS
        ln = 6; i1 = 3
        if self.tNCS: ln = 7; i1 = 4
        moments = self.extractnumbers(twinblock[ln], i1, 4)
        pvalues = moments[2:]
        # split acentric value into moment & SD
        v = moments[1].split('+/-')
        momentssd = [moments[0], v[0], v[1]]  # centric, acentric sd(acentric)

        twinxml = etree.Element('Twinning')
        centricxml = etree.SubElement(twinxml, 'CentricMoments')
        addElement(centricxml, 'Theoretical', theoretical[0])
        addElement(centricxml, 'TheoreticalTwin', theoreticaltwin[0])
        addElement(centricxml, 'Observed', momentssd[0])
        twinxml.append(centricxml)

        acentricxml = etree.SubElement(twinxml, 'AcentricMoments')
        addElement(acentricxml, 'Theoretical', theoretical[1])
        addElement(acentricxml, 'TheoreticalTwin', theoreticaltwin[1])
        addElement(acentricxml, 'Observed', momentssd[1])
        addElement(acentricxml, 'SDobserved', momentssd[2])
        twinxml.append(acentricxml)

        # P-values, untwinned, twin < 5%
        pval = pvalues[0]+"  "+pvalues[1]
        addElement(twinxml, 'Pvalues', pval)

        if 'Warning' in twinblock[12]:
            text = "\n".join(twinblock[12:13])
            addElement(twinxml, 'TwinWarning', text)

        self.xmlroot.append(twinxml)
    
    # . . . . . . . . . . . . . . . . . . . . . . . . . . .
    def anisotropy(self, lines):
        tag = 'Principal components of anisotropic part'
        skip = 2
        keep = 4
        # Search from end to get last occurance
        anisotropy = self.getBlock(lines, tag, skip, keep, False)
        eigenvalues = ''
        eigenvectors = ''
        for j in range(3):  # loop lines
            v = anisotropy[j].split()
            eigenvalues += ' '+v[0]
            eigenvectors += '   '
            eigenvectors += '   '.join(v[1:])
        #print("eigenvectors:", eigenvectors)
        
        deltaB = anisotropy[3].split()[7]

        anisotropyxml = etree.Element('Anisotropy')
        addElement(anisotropyxml, 'deltaB', deltaB)
        addElement(anisotropyxml, 'Eigenvalues', eigenvalues)
        addElement(anisotropyxml, 'Eigenvectors', eigenvectors)
        
        self.xmlroot.append(anisotropyxml)
                               
    # . . . . . . . . . . . . . . . . . . . . . . . . . . .
    # . . . . . . . . . . . . . . . . . . . . . . . . . . .
    def getXML(self):
        return self.xmlroot

    # . . . . . . . . . . . . . . . . . . . . . . . . . . .
    def getmaxresolution(self):
        # Return allresolution (A) from logfile
        return self.allresolution
    # . . . . . . . . . . . . . . . . . . . . . . . . . . .
    def extractnumbers(self, line, i1, nnum):
        # extract nnum numbers starting from i1 from line
        f = line.split()
        nums = []
        for i in range(nnum):
            nums.append(f[i1+i])
        return nums

    # . . . . . . . . . . . . . . . . . . . . . . . . . . .
    def findline(self, lines, text):
        # return first or last line containing text, else None
        for line in lines:
            if text in line:
                return line
        return None

    # . . . . . . . . . . . . . . . . . . . . . . . . . . .
    def getBlock(self, lines, text, skip, keep, forwardearch=True):
        # extract text block from lines
        #  look for "text", skip "skip" lines, then read "keep" lines
        # return block as list of lines
        # Find 1st line of block
        nlines = len(lines)
        i1 = -1
        if forwardearch:
            # Forward search
            for i in range(nlines):
                if text in lines[i]:
                    i1 = i
                    break
        else:
            # Backward search
            for i in reversed(range(nlines)):
                if text in lines[i]:
                    i1 = i
                    break
        
        # Now forwards
        block = []
        i2 = i1+skip
        for i in range(nlines)[i2:i2+keep]:
            block.append(lines[i])

        return block
# End class AnalysisLog
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
class Makexmlgraph:
    def __init__(self, graphid, loggraphdata):
        # Decode one loggraph from Phaser into XML CC4Table block
        self.graphid = graphid
        parts = loggraphdata.split('$$')
        
        # First part should contain TABLE and GRAPHS,
        # but may contain random line breaks
        t1 = " ".join(parts[0].splitlines())
        
        tg = t1.split('$GRAPHS ') # TABLE then GRAPHS
        
        # Title
        t = tg[0].split(':')
        title = ''
        if 'TABLE' in t[0]:
            title = t[1]
        self.xmlgraph = etree.Element('CCP4Table',
                                      attrib={'groupID':'Graph',
                                              'id': graphid,
                                              'title': title})
            
        # Graphs
        sg = tg[1].split(':')
        graphs = []

        props = []
        for f in sg:
            if f == ' ' or f == '':
                # blank field starts new graph
                if len(props) > 0:
                    graphs.append(props)
                    props = []
            else:
                props.append(f)

        for graph in graphs:
            plot = etree.Element('plot')
            # graph title
            addElement(plot, 'title', graph[0])
            cols = graph[2].split(',')
            xcol = cols[0]
            for ycol in cols[1:]:
                addElement(plot, 'plotline', None,
                                {'xcol':xcol, 'ycol':ycol})
            self.xmlgraph.append(plot)       
        
        # Headers are in next bit
        headers = parts[1]
        addElement(self.xmlgraph, 'headers', headers)
        collabels = headers.split()
        xcolabel = collabels[int(xcol)-1]
        if xcolabel in ['1/d^2']:  # any other labels?
            # x value is inverse resolution
            addElement(plot,'xscale','oneoversqrt')
            
        # Data is in penultimate part
        data = parts[3]

        if graphid == "2ndmoment":
            # the 2nd moment graph may go beyond meaningful data,
            # so strip out lines with zeroes
            col = 'ACent_obser'
            icol = collabels.index(col)
            data = self.stripfromend(data, icol)
            
        addElement(self.xmlgraph, 'data', data)


    def stripfromend(self, data, icol):
        lines = data.splitlines()

        # Check last line first
        fields = lines[-1].split()
        if fields[icol] != '0':
            # Nothing to do
            return data
        
        # Backward search
        nlines = len(lines)
        i1 = -1
        for i in reversed(range(nlines)):
            fields = lines[i].split()
            if fields[icol] != '0':
                i1 = i
                break

        if i1 < 0:
            i1 = 1

        data = "\n".join(lines[:i1+1])+"\n"

        return data

    def label(self):
      return self.label

    def getXML(self):
      return self.xmlgraph

