import pyrvapi
import os,sys,subprocess, time
from . import common

flush_time=None
output_dir=None
uri_prefix=None
verbose_err=True
init_meta={}

def call_pyrvapi(func, *args):
  # call the requested pyrvapi function, on pyrvapi error print warning and return False, otherwise return True
  # this ensures backwards compatibility with the older pyrvapi versions and assures crank2 does not stop at each tiny pyrvapi issue
  # if Andrey adds more specific exceptions to pyrvapi, this can be further improved
  try:
    getattr(pyrvapi,func)(*args)
  except Exception as e:
    if verbose_err:
      common.Warning('pyrvapi could not {0}{2} due to error: {1}'.format(func,e,args))
    return False
  return True

class rvbase(object):
  def __init__(self, ID, parent=None, param=None):
    self.ID=ID
    self.parent=parent
    self.param=param
    self.children=[]
  def Remove(self):
    if not call_pyrvapi('rvapi_remove_widget', self.ID):
      common.Warning('pyrvapi could not remove widget {0} {1}'.format(type(self).__name__,self.ID))
    elif self in self.parent.children:
      self.parent.children.remove(self)


def Document(ID, outdir, title=None, mode=1, jsuri=None, viewer=None, header=None, layout=7, \
              helpfname=None, htmlfname=None, taskfname=None, rvapidoc=None, i2=False, i2xmlfname="i2.xml"):
  if title is None:  title=ID
  if jsuri is None or viewer is None:
    if not "CCP4" in os.environ:
      common.Warning('pyrvapi could not initialize document because jsrview directory not defined.')
      return None
    if jsuri is None:
      jsuri = os.path.join(os.environ["CCP4"], "share", "jsrview")
    if viewer is None:
      if not i2:
        viewer = os.path.join(os.environ["CCP4"], "libexec", "jsrview")
      #  viewer = ['ccp4-python', os.path.join(os.environ["CCP4"], "share", "ccp4i2", "bin", "browser.py")]
      #### for testing purposes only
      #viewer = os.path.join("/bin/firefox")
  #if i2:  layout = 0
  if i2 and hasattr(i2,'rvapi_converter') and i2.rvapi_converter:
    global flush_time
    flush_time=0.0
    mode=2
  try:
    if rvapidoc:
      call_pyrvapi('rvapi_restore_document2',rvapidoc)
    else:
      call_pyrvapi('rvapi_init_document', ID, outdir, title, mode, layout, jsuri, helpfname, htmlfname, taskfname, i2xmlfname)
    # error may be returned with a delay
    time.sleep(0.1)
  # a specific exception should be caught once returned by pyrvapi
  except:
    common.Warning('pyrvapi could not initialize document, probably upgrade needed.')
    # backwards compatibility. should not be needed after some time.
    ID=ID+'_backwards'
    if not call_pyrvapi('rvapi_init_document', ID, outdir, title, mode, layout, jsuri, helpfname, htmlfname, taskfname):
      common.Warning('pyrvapi could not initialize document.')
      return None
  global output_dir
  output_dir = outdir
  global init_meta
  metadata=GetMetaData()
  if metadata:
    try:
      init_meta = eval(metadata)
    except:
      common.Warning('pyrvapi could not evaluate the initial meta.')
  if header and not call_pyrvapi('rvapi_add_header', header):
    common.Warning('pyrvapi could not add header.')
  if viewer and os.path.isfile(viewer):
    if common.is_string(viewer):  viewer=[viewer]
    with open(os.devnull, 'w') as fnull:
      proc=subprocess.Popen( viewer + [os.path.join(outdir, "index.html")], stderr=fnull, stdout=fnull)
  return rvdocument(ID)


class rvdocument(rvbase):

  def Tab(self, ID, title=None, opened=True):
    if title is None:  title=ID
    if not call_pyrvapi('rvapi_add_tab', ID, title, opened):
      common.Warning('pyrvapi could not add tab.')
      return None
    self.children.append( rvtab(ID,self) )
    return self.children[-1]

  def Grid(self, ID, filling=True, row=-1, col=0, rowspan=1, colspan=1):
    if not call_pyrvapi('rvapi_add_grid', ID, filling, self.ID, row, col, rowspan, colspan):
      common.Warning('pyrvapi could not add grid.')
      return None
    self.children.append( rvgrid(ID,self) )
    return self.children[-1]


class rvtab(rvbase):

  def Section(self, ID, title=None, row=-1, col=0, rowspan=1, colspan=1, opened=True):
    if title is None:  title=ID
    if not call_pyrvapi('rvapi_add_section', ID, title, self.ID, row, col, rowspan, colspan, opened):
      common.Warning('pyrvapi could not add section.')
      return None
    if ID not in [ch.ID for ch in self.children]:
      self.children.append( rvsection(ID,self) )
    return self.children[-1]

  def Tree(self, ID, title=None, row=-1, col=0, rowspan=1, colspan=1):
    if title is None:  title=ID
    if not call_pyrvapi('rvapi_add_tree_widget', ID, title, self.ID, row, col, rowspan, colspan):
      common.Warning('pyrvapi could not add section.')
      return None
    self.children.append( rvtree(ID,self) )
    return self.children[-1]

class rvgrid(rvtab):
  pass

class rvtree(rvtab):
  def Section(self, ID, title=None, row=-1, col=0, rowspan=1, colspan=1, opened=True):
    section = rvtab.Section(self, ID, title, row, col, rowspan, colspan, opened)
    if title is None:  title=ID
    if not call_pyrvapi('rvapi_set_tree_node', self.ID, section.ID, title, "auto", ""):
      common.Warning('pyrvapi could not add section to tree.')
      return None
    return section


class rvsection(rvbase):
  def Text(self, text, append=True, row=None, col=None, rowspan=None, colspan=None, flush=False):
    if row is not None or col is not None or rowspan is not None or colspan is not None:
      if append:
        append=False
      if row is None:  row = -1
      if col is None:  row = 0
      if rowspan is None:  rowspan = 1
      if colspan is None:  colspan = 1
    if append:
      if not call_pyrvapi('rvapi_append_text', text, self.ID):
        common.Warning('pyrvapi could not add text.')
    else:
      if not call_pyrvapi('rvapi_add_text', text, self.ID, row,col,rowspan,colspan):
        common.Warning('pyrvapi could not add text.')
    if flush:  Flush()

  def getgraph(self,graph):
    # if existing graph is not inputted then a new graph will be created, with one exception:
    #   if no ID is specified and a graph with the default ID exists already, this graph is returned
    if graph is None:
      default_graph_ID='g_'+self.ID
      if not [g for g in self.children if g.ID==default_graph_ID]:
        return self.Graph(default_graph_ID)
      else:
        return [g for g in self.children if g.ID==default_graph_ID][0]
    elif common.is_string(graph):
      return self.Graph(graph)
    else:
      return graph

  def PlotLine(self, X, Y, graph=None, block=None, plot=None ):
    graph = self.getgraph(graph)
    return graph.PlotLine(X, Y, block, plot)

  def Plot(self, ID, x_label, y_label=None, title=None, graph=None, block=None, legendloc=None, xmin=None, xmax=None, ymin=None, ymax=None, intx=None, inty=None):
    graph = self.getgraph(graph)
    return graph.Plot(ID, x_label, y_label, title, block, legendloc, xmin, xmax, ymin, ymax, intx, inty)

  def Graph(self, ID):
    if not call_pyrvapi('rvapi_append_loggraph', ID, self.ID):
      common.Warning('pyrvapi could not add graph.')
      return None
    self.children.append( rvgraph(ID,self) )
    return self.children[-1]

  def DataFile(self, fname, ftype, title=None, ID=None, folded=True, flush=False):
    if ID is None:
      ID = 'df_'+self.ID+str(len(self.children)+1)
    self.children.append( rvfileblock(ID,self) )
    self.children[-1].DataFile(fname, ftype, title, folded=folded)
    if flush:  Flush()
    return self.children[-1]

  def SetState(self, close=True):
    # returning 0 on failure
    if not call_pyrvapi('rvapi_set_section_state', self.ID, not close):
      common.Warning('pyrvapi could not set section state.')
      return 0
    return 1


# allow Coot to read in the map from Refmac automatically...  unfortunately, rvapi does not seem to support labels.
def SetMtzFWT(proc, mtz_obj):
  if not mtz_obj:
    return
  mtz_in=mtz_obj.GetFileName('mtz')
  mtz_out=mtz_obj.GetFileName('mtz').rsplit('.',1)[0]+'_deffwt.mtz'
  sft=proc.AddProg('sftools',propagate_out=False)
  sft.runname+='_set_fwt'
  sft.SetKey('read', mtz_in)
  sft.SetKey('calc', ('F','col','FWT','=','col','REFM_FWT'))
  sft.SetKey('calc', ('P','col','PHWT','=','col','REFM_PHWT'))
  sft.SetKey('write', mtz_out)
  sft.SetKey('quit\nY')
  sft.Run()
  proc.programs.remove(sft)
  return mtz_out

class rvfileblock(rvbase):
  empty = True

  def DataFile(self, fname, ftype, title=None, folded=True, flush=False):
    if uri_prefix:
      #common_prefix = os.path.commonprefix( [output_dir, fname] )
      #if common_prefix:
      #  fname = os.path.join( uri_prefix, os.path.relpath( fname, common_prefix ) )
      if output_dir and fname:
        if uri_prefix=='0':
          fname = os.path.relpath( fname, output_dir )
        else:
          fname = os.path.join( uri_prefix, os.path.relpath( fname, output_dir ) )
    # rvapi file paths should use slashes under windows since they are interpreted as URL's in case of 
    # server-client (eg jsCoFe) and in case of local jsrview, Qt functions will interpret them correctly too
    if os.name=='nt':
      fname = fname.replace('\\','/')
    if self.empty:
      if not call_pyrvapi('rvapi_add_data', self.ID, title, fname, ftype, self.parent.ID, -1, 0, 1, 1, folded):
        common.Warning('pyrvapi could not add data file.')
      else:
        self.empty=False
    else:
      if not call_pyrvapi('rvapi_append_to_data', self.ID, fname, ftype):
        common.Warning('pyrvapi could not append data file.')
    if flush:  Flush()



# the structure appears to be graph -> data block -> plot -> plotline

class rvgraph(rvbase):

  def getblock(self,block):
    if block is None:
      return self.DataBlock('db_'+self.ID)
    elif common.is_string(block):
      return self.DataBlock(block)
    else:
      return block

  def PlotLine(self, X, Y, block=None, plot=None ):
    block = self.getblock(block)
    return block.PlotLine(X, Y, plot)

  def Plot(self, ID, x_label, y_label=None, title=None, block=None, legendloc=None, xmin=None, xmax=None, ymin=None, ymax=None, intx=None, inty=None):
    if block is None:
      block = ID
    block = self.getblock(block)
    return block.Plot(ID, x_label, y_label, title, legendloc, xmin, xmax, ymin, ymax, intx, inty)

  def DataBlock(self, ID, title=None):
    if title is None:  title = ID
    if not call_pyrvapi('rvapi_add_graph_data', ID, self.ID, title):
      common.Warning('pyrvapi could not add data block.')
      return None
    self.children.append( rvdatablock(ID,self,0) )
    return self.children[-1]


class rvdatablock(rvbase):

  def PlotLine(self, X, Y, plot=None):
    if plot is None:
      plot = self.Plot(self.ID, X[0], Y[0])
    elif common.is_string(plot):
      plot = self.Plot(plot, X[0], Y[0])
    return plot.PlotLine(X[1:], Y[1:], self.ID)

  def Plot(self, ID, x_label, y_label=None, title=None, legendloc=None, xmin=None, xmax=None, ymin=None, ymax=None, intx=None, inty=None):
    if title is None:  title = ID
    if y_label is None: y_label = title
    if not call_pyrvapi('rvapi_add_graph_plot', ID, self.parent.ID, title, x_label, y_label):
      common.Warning('pyrvapi could not add graph plot.')
      return None
    self.children.append( rvplot(ID,self) )
    self.children[-1].SetProperties(legendloc, xmin, xmax, ymin, ymax, intx, inty)
    return self.children[-1]


class rvplot(rvbase):
  def_size=2.5
  def_color=["#4BB2C5","#EAA228"]
  def_color_fill="#B3B3B3"
  xmin,xmax,ymin,ymax=None,None,None,None

  # Y = (y_legend, y_header, [y_data]) ; y_header and/or [y_data] are optional (the same holds for X)
  # color,style,marker,size of plotline;  fill>0 creates filling under the curve, fill==1 filling only (plotline hidden)
  def PlotLine(self, X, Y, ID=None, color="", style="", marker="", size=0, fill=0, custom_x_tick=False, flush=False):
    # n should be ordernumber of plotline in the datablock
    n = str( self.parent.param + 1 )
    if ID is None:
      ID = self.parent.ID+n
    Xhead_omit, Yhead_omit = len(X)<2 or not common.is_string(X[1]), len(Y)<2 or not common.is_string(Y[1])
    X_header, Y_header = X[0]  if Xhead_omit  else X[1], Y[0]  if Yhead_omit  else Y[1]
    if not call_pyrvapi('rvapi_add_graph_dataset', "x"+n, self.parent.ID, self.parent.parent.ID, X[0], X_header):
      common.Warning('pyrvapi could not add X dataset.')
      return None
    if not call_pyrvapi('rvapi_add_graph_dataset', "y"+n, self.parent.ID, self.parent.parent.ID, Y[0], Y_header):
      common.Warning('pyrvapi could not add Y dataset.')
      return None
    if not call_pyrvapi('rvapi_add_plot_line', self.ID, self.parent.ID, self.parent.parent.ID, "x"+n, "y"+n):
      common.Warning('pyrvapi could not add plot line.')
      return None
    if fill:
      if not call_pyrvapi('rvapi_set_line_fill', "y"+n, self.ID, self.parent.ID, self.parent.parent.ID, bool(fill), bool(fill-1), "", 0.01):
        common.Warning('pyrvapi could not set fill options.')
      # setting grey color for automatic shadow fill without curve (showing up in the legend)
      if fill==1 and not color:
        color=self.def_color_fill
    if color or style or marker or size or (not color and len(self.children)<=1):
      if not size:  size=self.def_size
      if not color:
        if not self.children:                        color=self.def_color[0]
        elif len(self.children)<len(self.def_color): color=self.def_color[len(self.children)]
      if not call_pyrvapi('rvapi_set_line_options', "y"+n, self.ID, self.parent.ID, self.parent.parent.ID, color, style, marker, size, 1):
        common.Warning('pyrvapi could not set line options.')
    self.children.append( rvplotline(ID,self,n) )
    self.parent.param += 1
    x_data_ind, y_data_ind = 2  if not Xhead_omit  else 1, 2 if not Yhead_omit  else 1
    if len(X)>x_data_ind and len(Y)>y_data_ind:
      self.children[-1].Data(X[x_data_ind],Y[y_data_ind],custom_x_tick=custom_x_tick,flush=flush if flush is not None else True)
    return self.children[-1]

  def SetProperties(self, legendloc=None, xmin=None, xmax=None, ymin=None, ymax=None, intx=None, inty=None):
    if legendloc:
      if not call_pyrvapi('rvapi_set_plot_legend', self.ID, self.parent.parent.ID, legendloc, ''):
        common.Warning('pyrvapi could not set plot legend location.')
    # remembering min/max as rvapi resets them at each call!
    if xmin is not None:  self.xmin=xmin
    if xmax is not None:  self.xmax=xmax
    if ymin is not None:  self.ymin=ymin
    if ymax is not None:  self.ymax=ymax
    if self.xmin is not None and self.xmax is None:
      if not call_pyrvapi('rvapi_set_plot_xmin', self.ID, self.parent.parent.ID, self.xmin):
        common.Warning('pyrvapi could not set plot xmin.')
    elif self.xmax is not None and self.xmin is None:
      if not call_pyrvapi('rvapi_set_plot_xmax', self.ID, self.parent.parent.ID, self.xmax):
        common.Warning('pyrvapi could not set plot xmax.')
    elif self.xmax is not None and self.xmin is not None:
      if not call_pyrvapi('rvapi_set_plot_xrange', self.ID, self.parent.parent.ID, self.xmin, self.xmax):
        common.Warning('pyrvapi could not set plot xrange.')
    if self.ymin is not None and self.ymax is None:
      if not call_pyrvapi('rvapi_set_plot_ymin', self.ID, self.parent.parent.ID, self.ymin):
        common.Warning('pyrvapi could not set plot ymin.')
    elif self.ymax is not None and self.ymin is None:
      if not call_pyrvapi('rvapi_set_plot_ymax', self.ID, self.parent.parent.ID, self.ymax):
        common.Warning('pyrvapi could not set plot ymax.')
    elif self.ymax is not None and self.ymin is not None:
      if not call_pyrvapi('rvapi_set_plot_yrange', self.ID, self.parent.parent.ID, self.ymin, self.ymax):
        common.Warning('pyrvapi could not set plot yrange.')
    if intx is not None or inty is not None:
      if intx is None:  intx=0
      if inty is None:  inty=0
      if hasattr(pyrvapi,'rvapi_set_plot_int'):
        if not call_pyrvapi('rvapi_set_plot_int', self.ID, self.parent.parent.ID, intx, inty):
          common.Warning('pyrvapi could not set plot integer type.')
      else:
        # backwards compatibility.  This should be not needed after some time.
        common.Warning('pyrvapi function to set plot integer type does not exist. Upgrade needed.')

  def CustomTick(self,xy,val,label):
    if xy.startswith('x'):
      if not call_pyrvapi('rvapi_add_plot_xtick', self.ID, self.parent.parent.ID, val, str(label)):
        common.Warning('pyrvapi could not add x plot tick {0}.'.format(val))
    if xy.startswith('y'):
      if not call_pyrvapi('rvapi_add_plot_ytick', self.ID, self.parent.parent.ID, val, str(label)):
        common.Warning('pyrvapi could not add y plot tick {0}.'.format(val))


class rvplotline(rvbase):

  def singleData(self,xy,val,typ=float,custom_tick=False):
    n=self.param
    # implemented, although there seems to be no difference between real and int plotting?
    if typ is int:
      if not call_pyrvapi('rvapi_add_graph_int', xy+n, self.parent.parent.ID, self.parent.parent.parent.ID, int(val)):
        common.Warning('pyrvapi could not add {0} integer data {1}.'.format(xy,val))
    else:
      if not call_pyrvapi('rvapi_add_graph_real', xy+n, self.parent.parent.ID, self.parent.parent.parent.ID, float(val), "%g"):
        common.Warning('pyrvapi could not add {0} real data {1}.'.format(xy,val))
    if custom_tick is not False:
      self.CustomTick(xy,val,custom_tick)

  def CustomTick(self,xy,val,label):
    self.parent.CustomTick(xy,val,label)
    #if xy.startswith('x'):
    #  if not call_pyrvapi('rvapi_add_plot_xtick', self.parent.ID, self.parent.parent.parent.ID, val, str(label)):
    #    common.Warning('pyrvapi could not add x plot tick {0}.'.format(val))
    #if xy.startswith('y'):
    #  if not call_pyrvapi('rvapi_add_plot_ytick', self.parent.ID, self.parent.parent.parent.ID, val, str(label)):
    #    common.Warning('pyrvapi could not add y plot tick {0}.'.format(val))

  def Reset(self, reset_data=True, reset_ticks=True):
    n=self.param
    if reset_data:
      if not call_pyrvapi('rvapi_reset_graph_dataset', "x"+n, self.parent.parent.ID, self.parent.parent.parent.ID):
        common.Warning('pyrvapi could not remove dataset x{0}.'.format(n))
      if not call_pyrvapi('rvapi_reset_graph_dataset', "y"+n, self.parent.parent.ID, self.parent.parent.parent.ID):
        common.Warning('pyrvapi could not remove dataset y{0}.'.format(n))
    if reset_ticks:
      if not call_pyrvapi('rvapi_reset_plot_xticks', self.parent.ID, self.parent.parent.parent.ID):
        common.Warning('pyrvapi could not remove xticks.')
      if not call_pyrvapi('rvapi_reset_plot_yticks', self.parent.ID, self.parent.parent.parent.ID):
        common.Warning('pyrvapi could not remove yticks.')

  def Data(self,x,y, x_type=float, y_type=float, custom_x_tick=False, custom_y_tick=False, flush=False):
    try:
      length = min( len(x),len(y) )
    except TypeError:
      if custom_x_tick is False:
        self.singleData("x",x,x_type,custom_tick=False)
      else:
        self.singleData("x",custom_x_tick,x_type,custom_tick=x)
      if custom_y_tick is False:
        self.singleData("y",y,y_type,custom_tick=False)
      else:
        self.singleData("y",custom_y_tick,y_type,custom_tick=y)
    else:
      for i in range( length ):
        if y[i] is not None:
          if custom_x_tick is False:
            self.singleData("x",x[i],x_type,custom_tick=False)
          else:
            self.singleData("x",i+1,x_type,custom_tick=x[i])
          if custom_y_tick is False:
            self.singleData("y",y[i],y_type,custom_tick=False)
          else:
            self.singleData("y",i+1,y_type,custom_tick=y[i])
    if flush:  Flush()


def AppendContent(file_path, parent_id):
  if not call_pyrvapi('rvapi_append_content', file_path, True, parent_id):
    common.Warning('pyrvapi could not append content.')

def PutMetaData(metadata, file_path=None):
  if not call_pyrvapi('rvapi_put_meta', metadata):
    common.Warning('pyrvapi could not put metadata.')
  if file_path:
    if not call_pyrvapi('rvapi_store_document2', file_path):
      common.Warning('pyrvapi could not store document 2.')

def GetMetaData():
  if not call_pyrvapi('rvapi_get_meta'):
    common.Warning('pyrvapi could not put metadata.')
    return {}
  meta = pyrvapi.rvapi_get_meta()
  return(meta)

def Flush(ignore_timing_restr=False, timing_restr=3.0):
  # workaround for flush flooding issues of jsrview: flush frequency set to at least 1 second
  # feature disabled for rvapi as treated directly by jsrview/rvapi;  needed for rvapi->i2 native though
  global flush_time
  if ignore_timing_restr or flush_time is None or (time.time()-flush_time)>timing_restr:
    if not call_pyrvapi('rvapi_flush'):
      common.Warning('pyrvapi could not flush.')
    elif flush_time is not None:
      flush_time = time.time()
