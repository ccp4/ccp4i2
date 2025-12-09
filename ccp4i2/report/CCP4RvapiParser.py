from __future__ import print_function

import sys


'''
Version history:
5 Mar 2018, AL:
- original graphs widget code by JA refactored
  with added support for multi-dataset widgets,
- tables are collected from separate files and
  are shown verbatim with class changed to fancy,
- customisation functions and class variables added
  to control basic layout parameters from subclasses,
- currently used by crank2, morda and arpwarp
  (no changes to crank2 code were required)
'''

import os, re
#from lxml import etree as ET
import xml.etree.ElementTree as ET
from ccp4i2.report.CCP4ReportParser import Report,CCP4NS,PARSER

class RvapiReport(Report):
  MAINTAINER    = 'andrey.lebedev@stfc.ac.uk'
  TASKNAME      = None
  RUNNING       = True
  USEPROGRAMXML = False
  WATCHED_FILE  = 'i2.xml'
  RVAPI_XML     = 'i2.xml'
  RVAPI_DIR     = '.'
  RVAPI_MERGE_TABS = True
  RVAPI_STD_HEADER = True

  def __init__(self, jobInfo={}, **kw):
    try:
      work_absdir = jobInfo['fileroot']
      data_dir = 'tables_as_xml_files'
      data_fmt = 'insert_%s.xml'
      rvapi_absdir = os.path.join(work_absdir, self.RVAPI_DIR)
      rvapi_absdir = os.path.abspath(rvapi_absdir)
      self.rvapi_absdir = rvapi_absdir
      i2xml_absfile = os.path.join(work_absdir, self.RVAPI_XML)
      i2xml_absfile = os.path.abspath(i2xml_absfile)
      data_absdir = os.path.join(work_absdir, data_dir)
      data_reldir = os.path.join('.', data_dir)
      data_absfmt = os.path.join(data_absdir, data_fmt)
      data_relfmt = os.path.join(data_reldir, data_fmt)
      if not os.path.isdir(data_absdir):
        os.mkdir(data_absdir)

      e0 = ET.parse(i2xml_absfile).getroot()
      self.clear_indents(e0)
      table_dict = self.collect_table_data(rvapi_absdir)
      l0 = self.translate_section_data(e0, table_dict)
      self.write_data_files(e0, data_absfmt, data_relfmt)
      self.SEPARATEDATA = not self.footprint_changed(l0, work_absdir)
      self.init_section_layouts()
      self.set_styles()

    except:
      exc_type, exc_value,exc_tb = sys.exc_info()[:3]
      sys.stderr.write(str(exc_type)+'\n')
      sys.stderr.write(str(exc_value)+'\n')
      import traceback
      traceback.print_tb(exc_tb)
    
      return

    if not self.RUNNING and self.SEPARATEDATA:
      return

    Report.__init__(self, jobInfo=jobInfo, **kw)
    self.inject_css()

    if self.RVAPI_STD_HEADER:
      r0 = self.addResults()

    else:
      r0 = self
      r0.addDiv(style='clear:both;')

    if self.RVAPI_MERGE_TABS:
      toplevel = 'tab/section'

    else:
      toplevel = 'tab'

    for e1 in e0.findall(toplevel):
      self.append_section(e1, r0)

  def test(self, htmlout):
    htmlstr = ET.tostring(self.as_etree().decode("utf-8"))
    htmlstr = htmlstr.replace('><', '>\n<')
    docsdir = os.path.join(os.environ['CCP4'], 'share', 'ccp4i2', 'docs')
    htmlstr = htmlstr.replace(os.path.join('..', '..'), docsdir)
    with open(htmlout, 'w') as ostream:
      print(htmlstr, file=ostream)

  # overwrite
  def show_section_fold(self, e1, count):
    return count.total > 0 and e1.get('id') != 'summary'

  # overwrite
  def show_table_legends(self, e1, count):
    return count.float > 1

  # overwrite
  def get_section_layout(self, e1, count):
    return self.twocolumn if count.graph > 0 else self.onecolumn

  # overwrite
  def graph_geometry(self):
    return 'height:310px;width:400px;'

  # overwrite
  def set_styles(self):
    border = ' border: 1px solid black;'
    border = ''
    border += ' border-spacing: 0px; padding-right: 0px; padding-top: 0px;'
    self.styles = {
      'small': 'font: 12px Arial, sans-serif; margin-top: 44px; background-color: rgb(145,203,219); padding: 5px;',
      'div.rvapi-after-fold': 'width: 100px; height: 5px;' + border,
      'table.rvapi-page': 'padding-left: 9px; padding-bottom: 8px;' + border,
      'td.rvapi-column-left': 'text-align: left; vertical-align: top; width: 200px;' + border,
      'td.rvapi-column-right': 'text-align: right; vertical-align: top; width: 425px;' + border,
      'td.rvapi-column-single': 'text-align: left; vertical-align: top; width: 625px;' + border,
    }

  def inject_css(self):
    for style_tuple in sorted(self.styles.items()):
      self.append('<style type="text/css">%s {%s}</style>' %style_tuple)

  def init_section_layouts(self):
    f0 = ET.Element('table')
    f1 = ET.SubElement(ET.SubElement(f0, 'tbody'), 'tr')
    f2l = ET.SubElement(f1, 'td')
    f2r = ET.SubElement(f1, 'td')
    f0.set('class', 'rvapi-page')
    f2l.set('class', 'rvapi-column-left')
    f2r.set('class', 'rvapi-column-right')
    self.twocolumn = f0, f2l, f2r

    f0 = ET.Element('table')
    f1 = ET.SubElement(ET.SubElement(f0, 'tbody'), 'tr')
    f2 = ET.SubElement(f1, 'td')
    f0.set('class', 'rvapi-page')
    f2.set('class', 'rvapi-column-single')
    self.onecolumn = f0, f2, f2

  def append_separator(self, r1):
    f0 = ET.Element('div')
    f0.set('class', 'rvapi-after-fold')
    r1.append(ET.tostring(f0).decode("utf-8"))

  def append_section(self, e1, r0):
    count = SectionCounts(e1)
    if self.show_section_fold(e1, count):
      title = e1.findtext('name', 'Untitled')
      state = e1.findtext('open', 'false') == 'true'
      r1 = r0.addFold(label=title, initiallyOpen=state)
      self.append_separator(r1)

    else:
      r1 = r0

    f0, f2_left, f2_right = self.get_section_layout(e1, count)
    show_legends = self.show_table_legends(e1, count)
    e2_list = list(e1)
    e2_list.append(None)
    e2 = e2_list.pop(0)
    while e2 is not None:
      cou = 0
      while e2 is not None and e2.tag == 'section':
        self.append_section(e2, r1)
        e2 = e2_list.pop(0)
        cou += 1

      if cou:
        self.append_separator(r1)

      while e2 is not None and e2.tag != 'section':
        if e2.tag == 'table':
          if show_legends:
            legend = e2.findtext('legend')
            if legend:
              ET.SubElement(f2_left, 'p').text = legend

          f2_left.append(self.table_div(e2.get('data_path'), e2.get('id')))

        elif e2.tag == 'text':
          f3 = ET.SubElement(f2_left, 'p')
          f3.append(e2)
          ET.SubElement(f3, 'br')

        elif e2.tag == 'graph':
          data_absfile = os.path.join(self.rvapi_absdir,e2.get('data_path'))
          # the above line breaks arpwarp
          # not sure if this was intentional or a bug
          # hence if
          if not os.path.isfile(data_absfile):
            data_absfile = e2.get('data_abspath')

          with open(data_absfile) as dataFile:
              theData = dataFile.read()
          theDataTree = ET.fromstring(theData)
          tables = theDataTree.findall(".//CCP4Table")

          for i in range(len(tables)):
              dataDivName = "div_data_"+e2.get('data_path').lstrip("./").rstrip(".lxml").replace("/","_")+str(i)
              dataDiv = ET.Element('ccp4_data',nsmap = {"ccp4":CCP4NS}, **{
                  'id': dataDivName,
                  'type': 'application/xml',
                  'style': 'display:none',
                  'title': tables[i].attrib["title"],
              })

              try:
                 #lxml
                 dataDiv.extend(tables[i].getchildren())
              except:
                 try:
                 #xml
                     dataDiv.extend(tables[i].findall("./"))
                 except:
                    #Give up!
                    pass

              f2_right.append(self.graph_div(dataDivName, e2.get('id')+str(i)))
              f2_right.append(dataDiv)

        e2 = e2_list.pop(0)

      if len(f2_left) or len(f2_right):
        r1.append(ET.tostring(f0).decode("utf-8"))
        f2_left[:] = ()
        f2_right[:] = ()

  def table_div(self, file_relpath, id1):
    return ET.Element('div', **{
      'data-data': file_relpath,
      'data-initially-drawn': 'True',
      'data-is-urls': 'True',
      'data-renderer': 'CCP4i2Widgets.CCP4i2HTMLChunk',
      'data-require': 'CCP4i2Widgets',
      'data-widget-type': 'CCP4i2DrawnDiv',
      'id': 'div_' + id1,
    })

  def graph_div(self, data_div, id1):
    return ET.Element('div', **{
      'data-app-font-size': '12px',
      'data-data': data_div,
      'data-initially-drawn': 'True',
      'data-is-urls': 'False',
      'data-renderer': 'CCP4i2Widgets.CCP4FlotRenderer',
      'data-require': 'CCP4i2Widgets',
      'data-widget-type': 'CCP4i2DrawnDiv',
      'id': 'div_' + id1,
      'style': self.graph_geometry() + 'display:inline-block;'
    })

  def collect_table_data(self, rvapi_absdir):
    t1_dict = dict()
    tsk_absfile = os.path.join(rvapi_absdir, 'task.tsk')
    if not os.path.isfile(tsk_absfile):
      return

    t1_list = list()
    with open(os.path.join(rvapi_absdir, 'task.tsk')) as istream:
      for s1 in re.findall('<table .+?</table>', istream.read(), re.S):
        t1 = ET.fromstring(s1)
        if len(t1):
          t1_list.append(t1)

    for f1 in os.listdir(rvapi_absdir):
      if f1.endswith('.table'):
        t1 = ET.parse(os.path.join(rvapi_absdir, f1)).getroot()
        if len(t1):
          t1_list.append(t1)

    for t1 in t1_list:
      id1 = t1.get('id', None)
      if id1 and id1.endswith('-grid'):
        id1 = id1[:-5]
        t1.set('id', id1)
        t1_dict[id1] = t1

    return t1_dict

  def translate_table_data(self, e1, t1_dict):
    id1 = e1.get('id')
    t1 = t1_dict[id1]
    tags = [t2.tag for t2 in t1]
    if tags == ['thead', 'tbody']:
      assert len(t1) == 2
      f1 = t1

    else:
      tset = set(tags)
      tag = tset.pop()
      assert not tset and tag == 'tr'
      f1 = ET.Element('table')
      f1.append(t1)
      f1.attrib.update(t1.attrib)
      t1.attrib.clear()
      t1.tag = 'tbody'

    assert id1 == f1.get('id')
    for f2 in f1.iter():
      f2.attrib.pop('class', None)

    f1.find('tbody').set('class', 'fancy')
    f1.tail = '\n'
    legend = e1.findtext('legend')
    del e1[:]
    ET.SubElement(e1, 'file_content').append(f1)
    if legend:
      ET.SubElement(e1, 'legend').text = legend

  def translate_graph_data(self, e1):
    f1 = ET.Element('CCP4ApplicationOutput')
    e2_list = list()
    e2_dict = dict()
    for e2 in e1.findall('graphdata'):
      e3_list = list()
      e2_dict[e2.get('id')] = e3_list
      e2_list.append((e2, e3_list))

    for e3 in e1.findall('plot'):
      if e3.find('plotline') is not None:
        key = e3.find('plotline').get('dataset')
        e3_list = e2_dict[key]
        e3_list.append(e3)

    for e2, e3_list in e2_list:
      id2 = e2.get('id')
      e3 = e2.find('setids')
      idc_list = e3.text.split(e3.get('separator'))
      f2 = ET.SubElement(f1, 'CCP4Table', **{
        'title': e2.findtext('title'),
        'id': e2.get('id'),
        'style': 'display:none;',
      })
      e3 = e2.find('setnames')
      f3_head = e3.text.replace(' ', '&nbsp;').replace(e3.get('separator'), ' ')
      ET.SubElement(f2, 'headers').text = f3_head
      f3_data = re.sub(' ?\n ?', '\n', re.sub(' +', ' ', e2.findtext('data')))
      ET.SubElement(f2, 'data').text = f3_data
      for e3 in e3_list:
        f3 = ET.SubElement(f2, 'plot')
        ET.SubElement(f3, 'title').text = e3.findtext('title')
        ET.SubElement(f3, 'plottype').text = 'xy'
        ET.SubElement(f3, 'xlabel').text = e3.findtext('xname')
        ET.SubElement(f3, 'ylabel').text = e3.findtext('yname')
        ET.SubElement(f3, 'xintegral').text = e3.findtext('intx')
        ET.SubElement(f3, 'yintegral').text = e3.findtext('inty')
        cxval = ','.join(e3.findtext('custom_xvalues', '').split())
        cxlab = ','.join(e3.findtext('custom_xlabels', '').split())
        if cxval and cxlab:
          ET.SubElement(f3, 'customXTicks').text = cxval
          ET.SubElement(f3, 'customXLabels').text = cxlab

        f4 = ET.Element('xrange')
        e4 = e3.find('xmax')
        if e4.get('on') == 'true':
          f4.attrib['max'] = e4.text

        e4 = e3.find('xmin')
        x_min = e4.text
        if e4.get('on') == 'true':
          f4.attrib['min'] = x_min

        if len(f4.attrib):
          f3.append(f4)

        f4 = ET.Element('yrange')
        e4 = e3.find('ymax')
        y_max = e4.text
        if e4.get('on') == 'true':
          f4.attrib['max'] = y_max

        e4 = e3.find('ymin')
        y_min = e4.text
        if e4.get('on') == 'true':
          f4.attrib['min'] = y_min

        if len(f4.attrib):
          f4.attrib['rightaxis'] = 'false'
          f3.append(f4)

        for e4 in e3.findall('plotline'):
          if e4.get('dataset') == id2:
            xind = idc_list.index(e4.get('xset'))
            yind = idc_list.index(e4.get('yset'))
            if e4.findtext('fill') == 'true':
              all_cols = list(zip(*[line.split() for line in f3_data.split('\n')]))
              xy_rows = list(zip(all_cols[xind], all_cols[yind]))
              xy_rows = [(x, y) for x, y in xy_rows if x + y != '..']
              x1, y1 = xy_rows[0]
              x2, y2 = xy_rows[-1]
              if y1 == y2:
                xpos = float(x1) + float(x2)/ 5.0
                if float(y_max) > 0.0 and float(y2) > float(y_max):
                  ypos = float(y_max)/ 2.1

                else:
                  ypos = float(y2) + float(y2)/ 45.0

                f4 = ET.SubElement(f3, 'polygon')
                f4.text = ' '.join((x1, y_min, x1, y1, x2, y2, x2, y_min))
                f4.set('linecolour', '#cccccc')
                f4.set('fillcolour', '#aaaaaa')
                f4.set('alpha', '0.2')
                f4 = ET.SubElement(f3, 'text')
                f4.text = f3_head.split()[yind].replace('&nbsp;', ' ')
                f4.set('xpos', str(xpos))
                f4.set('ypos', str(ypos))
                f4.set('font', 'bold 1.25em Helvetica')
                f4.set('colour', '#aaaaaa')

            else:
              xcol = str(xind + 1)
              ycol = str(yind + 1)
              line_style = e4.findtext('style')
              line_color = e4.findtext('color')
              if line_style == 'bars':
                f4 = ET.SubElement(f3, 'barchart', col=xcol, tcol=ycol)
                ET.SubElement(f4, 'width').text = '1'
                if line_color:
                  ET.SubElement(f4, 'colour').text = line_color

              else:
                f4 = ET.SubElement(f3, 'plotline', xcol=xcol, ycol=ycol)
                ET.SubElement(f4, 'linestyle').text = '.' if line_style == 'off' else '-'
                ET.SubElement(f4, 'symbolsize').text = str(int(float(e4.findtext('width'))))
                if line_color:
                  ET.SubElement(f4, 'colour').text = line_color
                  ET.SubElement(f4, 'fillcolour').text = line_color

    self.set_indents(f1, '\n', '  ')
    id1 = e1.get('id')
    del e1[:]
    ET.SubElement(e1, 'file_content').append(f1)

  def translate_section_data(self, e0, table_dict):
    e1_list = list()
    key_e1_list = list()
    cou0 = 0
    cou3 = 0
    for e1 in e0:
      key = None
      if e1.tag == 'tab':
        assert e1.get('id')
        cou0 += 1
        key = cou0, 0, 0, 0
        cou0 += 1

      elif e1.get('id'):
        row = e1.get('row', '_')
        col = e1.get('col', '_')
        if e1.get('id') and row.isdigit() and col.isdigit():
          cou3 -= 1
          key = cou0, int(row), int(col), cou3

      if key:
        key_e1_list.append((key, e1))

      else:
        e1_list.append(e1)

    e0[:] = e1_list
    l0 = list()
    for key, e1 in sorted(key_e1_list):
      e0.append(e1)
      l1 = None
      if e1.tag in ('tab', 'section'):
        l1 = self.translate_section_data(e1, table_dict)

      elif e1.tag == 'table':
        self.translate_table_data(e1, table_dict)

      elif e1.tag == 'graph':
        self.translate_graph_data(e1)

      l0.append((key, e1.tag, l1))

    return l0

  def clear_indents(self, e1):
    e1.text = e1.text.strip() if e1.text else None
    e1.tail = e1.tail.strip() if e1.tail else None
    for e2 in e1:
      self.clear_indents(e2)

  def set_indents(self, e1, sep1, sep_add):
    if not len(e1):
      return

    if e1.text:
      return

    for e2 in e1:
      if e2.tail:
        return

    sep2 = sep1 + sep_add
    sep_var = sep1
    for e2 in reversed(e1):
      self.set_indents(e2, sep2, sep_add)
      e2.tail = sep_var
      sep_var = sep2

    e1.text = sep_var
    e1.tail = '\n'

  def write_data_fileso(self, e0, data_absfmt, data_relfmt):
    e1 = e0.find('file_content') # Original variant.
    if e1 is None:
      for e1 in e0:
        self.write_data_files(e1, data_absfmt, data_relfmt)

    else:
      assert len(e1) == 1
      id1 = e0.get('id')
      data_absfile = data_absfmt %id1
      data_relfile = data_relfmt %id1
      ET.ElementTree(e1).write(data_absfile)  # This line is causing trouble. Changed e1[0] to e1.
      e0.remove(e1)
      del e1
      e0.set('data_path', data_relfile)
      e0.set('data_abspath', data_absfile)

  def write_data_files(self, e0, data_absfmt, data_relfmt):
    feles = e0.findall('.//file_content') # re-write of previous function.
    for ele in feles:
      pnd = [ x for x in e0.findall(".//file_content/...") if ele in x.findall("file_content")][0]
      id1 = pnd.get('id')
      data_absfile = data_absfmt %id1
      data_relfile = data_relfmt %id1
      ET.ElementTree(ele).write(data_absfile)
      pnd.remove(ele)
      del ele
      pnd.set('data_path', data_relfile)
      pnd.set('data_abspath', data_absfile)

  def footprint_changed(self, l0, work_absdir):
    fp_absfile = os.path.join(work_absdir, 'footprint')
    fp_old = None
    if os.path.isfile(fp_absfile):
      with open(fp_absfile) as istream:
        fp_old = istream.read()

    fp_changed = False
    fp_new = repr(l0)
    if fp_new != fp_old:
      fp_changed = True
      with open(fp_absfile, 'w') as ostream:
        ostream.write(fp_new)

    return fp_changed

class SectionCounts(object):

  def __init__(self, e1):
    self.section = len(e1.findall('section'))
    self.graph = len(e1.findall('graph'))
    self.table = len(e1.findall('table'))
    self.float = self.graph + self.table
    self.text = len(e1.findall('text'))
    self.total = self.section + self.float + self.text

