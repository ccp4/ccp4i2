import sys
import getopt
import traceback
import json
from lxml import etree
import numpy

def parseGraphData(tree):

    def fixNan(arr):
       newArr = []
       for v in arr:
           if numpy.isnan(v):
               newArr.append(None)
           else:
               newArr.append(v)
       return newArr

    t = tree
    #print(etree.tostring(t,pretty_print=True))
    ttitle = t.attrib.get('title','').strip()
    data = t.xpath('data')[0]
    a = data.text.strip().split()
    a2 = []
    for i in a:
        j = i.replace("NA","Nan")
        if j == "*" or j=="-":
             j = "NaN"
        a2.append(j)

    array1 = numpy.array(a2,dtype=float)
    ncols = len(data.text.strip().split("\n")[0].split())
    nrows = len(array1) / ncols
    array2 = numpy.transpose(numpy.array_split(array1,len(array1)/ncols))
    graph_lists = []
    ydata = []
    headers = []
    for i,p in zip(list(range(len(t.xpath('plot')))),t.xpath('plot')):
        if len(p.xpath("title")) > 0:
            ptitle = p.xpath("title")[0].text.strip()
        else:
            ptitle = 'Plot %d' % (i+1)
        graph_lists.append((ptitle,""))

    headers_el = t.xpath("headers")
    if len(headers_el)>0:
        header_sep = headers_el[0].attrib.get("separator"," ")
        if header_sep == " ":
            headers = headers_el[0].text.split()
        else:
            headers = headers_el[0].text.split(header_sep)

    all_plots = []
    for p in t.xpath('plot'):
        plots = []
        for plot in p.xpath('plotline'):
            xcol = plot.attrib.get("xcol",None)
            ycol = plot.attrib.get("ycol",None)
            if xcol is not None and ycol is not None:
                x_idx, y_idx = (int(xcol)-1,int(ycol)-1)
                print(array2[y_idx].tolist())
                x_data = fixNan(array2[x_idx].tolist())
                y_data = fixNan(array2[y_idx].tolist())
                thePlot = {"x_data":x_data,"y_data":y_data}
                if(len(headers)==len(array2)):
                    thePlot["x_header"] = headers[x_idx]
                    thePlot["y_header"] = headers[y_idx]
                else:
                    thePlot["x_header"] = "x"
                    thePlot["y_header"] = "y"
                plots.append(thePlot)
        all_plots.append(plots)

    return all_plots

def openFileToEtree(fileName=None,printout=True):
    # Use this as etree.parse() seg faults on some Linux
    parser = etree.HTMLParser()
    f = open(fileName)
    s = f.read()
    f.close()
    tree = etree.fromstring(s, parser)
    return tree

def addCCP4ReportFile(fname,select=None):
    try:
        doc = openFileToEtree(fname)
    except:
        exc_type, exc_value,exc_tb = sys.exc_info()[:3]
        sys.stderr.write(str(exc_type)+'\n')
        sys.stderr.write(str(exc_value)+'\n')
        traceback.print_tb(exc_tb)
        return []
    graphList = []
    for tag in ('ccp4_data','{http://www.ccp4.ac.uk/ccp4ns}ccp4_data'):
        for tableEle in doc.iter(tag = tag):
            if select is None or tableEle.get('id',None) in select:
                try:
                    graphList.append(tableEle)
                except:
                    print('ERROR loading table',tableEle.tag)

    # This section deals with data coming from external XML files.
    for tableEle in doc.iter(tag = "div"):
        dataData = tableEle.attrib.get('data-data','')
        # We allow for the bizarre situation where the string ends with .xml but is actually a div id.
        if len(dataData.strip())>0 and dataData.endswith(".xml") and len(doc.findall(".//*[@id='"+dataData+"']"))==0:
            fnameXML = os.path.join(os.path.dirname(fname),dataData)
            try:
                docXML = openFileToEtree(fnameXML)
                tables = docXML.xpath('/CCP4ApplicationOutput/CCP4Table')
                for t in tables:
                    graphList.append(tableEle)
            except:
                print('POSSIBLE ERROR loading table file',fnameXML,".")
                print("This may not be a table file, in which case there is no problem.")
    return graphList

def extractGraphData(file_names,select,theGraphTotal):
    all_graphs = []
    for f in file_names:
        graphs = addCCP4ReportFile(f,select=select)
        for g in graphs:
            newG = parseGraphData(g)
            all_graphs.extend(newG)

    if theGraphTotal is not None and theGraphTotal<len(all_graphs):
        return [all_graphs[theGraphTotal]]
    else:
        return all_graphs

if __name__ == "__main__":
    if "-quit" in sys.argv:
        sys.exit()

    hist = False
    nomainwin = False
    testla = False
    select = []

    import getopt

    try:
        opts, file_names = getopt.getopt(sys.argv[1:], "x:s:g:G:o", ["select=","graph=","graphTotal=","output="])
    except getopt.GetoptError as err:
        # print help information and exit:
        print(str(err)) # will print something like "option -a not recognized"
        sys.exit(2)

    theGraph = None
    theGraphTotal = None
    output_file = None

    for o, a in opts:
        if o in ("-s", "--select"):
            items = a.strip('[').strip(']').split(',')
            for item in items: select.append(item.strip("'").strip('"'))
        elif o in ("-g", "--graph"):
            items = a.strip('[').strip(']').split(',')
            theGraph = items
        elif o in ("-G", "--graphTotal"):
            theGraphTotal = int(a)
        elif o in ("-o", "--output"):
            output_file = a
        else:
            assert False, "unhandled option: "+o

    graph_data = extractGraphData(file_names,select,theGraphTotal)

    if output_file:
        with open(output_file,"w+") as of:
            of.write(json.dumps(graph_data))
    else:
        print(json.dumps(graph_data))
