import os,sys,io
import common,manager,process
import pkgutil,csv
try:
  import argparse
  no_argparse=False
except:
  import optparse
  no_argparse=True

program="crank2"


csv.register_dialect('crank_parse', delimiter=' ', skipinitialspace=1, doublequote=False, quotechar='"')


class parse:
  share_with_crank = ('hklin','xyzin','seqin','hklout','xyzout','xyzsubout','lgout','logout',
                      'disable_rvapi','rvapi_viewer','rvapi_uri_prefix','rvapi_no_tree','rvapi_separate_steps','rvapi_document','disable_mtz_label_prefix')

  def __init__(self,dummy=False):
    self.f_in=None
    if not dummy:
      self.ParseCommandLine([])

  def ParseAndRun(self, ccp4i2crank=None, ccp4i2defaults=None, rvapi_style=None):
    """Parse command line and input file (if supplied), create the input xml and run crank"""
    self.pars_arg.rvapi_separate_steps = False
    if rvapi_style and rvapi_style not in ('tree','rvapi_tree'):
      self.pars_arg.rvapi_no_tree = True
      if rvapi_style and rvapi_style=='i2_tree':
        self.pars_arg.rvapi_separate_steps = True
    # create dummy crank object and copy information from it to XML file
    crank_prep = self.ParseAll(dummy=True, skip_commline_parse=ccp4i2crank)
    if (not self.pars_arg.xmlin or self.pars_arg.gcx) and not self.pars_arg.write_defaults:
      common.WriteXML( crank_prep, self.pars_arg.xmlout )
    if self.pars_arg.gcx:
      common.Info('Input XML file created.')
      return None
    # create crank run object 
    crank_run=self.ParseAll(skip_commline_parse=ccp4i2crank)
    crank_prep.run=crank_run
    # write complete object incl. defaults to XML if asked for it
    if self.pars_arg.write_defaults or ccp4i2defaults:
      crank_prep.disable_rvapi = True
      crank_prep.Run(emulate=True, set_all_par=True)
      if self.pars_arg.write_defaults:
        common.WriteXML( crank_prep.run, self.pars_arg.xmlout )
        common.Info('An input XML file incl. all defaults and passed objects created.')
      if ccp4i2defaults:
        return crank_run
    else:
      # connect to ccp4i2 (if asked to) and run crank (from the dummy prep object)
      if ccp4i2crank:
        crank_run.ccp4i2=ccp4i2crank
      # catch sys.exit(0)
      try:
        crank_prep.Run()
      except SystemExit as e:
        if e.code!=0:  raise
      except common.Unsuccessful as e:
        if not self.pars_arg.graceful_preliminary_stop:  raise
      return crank_run
    return None

  def ParseAll(self, dummy=False, skip_commline_parse=False):
    """Parse command line and input file (if supplied) and return crank process on success"""
    if not skip_commline_parse:
      self.ParseCommandLine()
    if self.pars_arg.process:
      self.PrintParams(self.pars_arg.process)
      sys.exit()
    if self.pars_arg.list_processes:
      self.PrintProcesses()
      sys.exit()
    if self.pars_arg.list_data_objects:
      self.PrintDataObjects()
      sys.exit()
    rundir = 'crank2'  if not self.pars_arg.dirout  else self.pars_arg.dirout
    rundir = os.path.realpath(rundir)
    if self.pars_arg.xmlin:
      crank = manager.crank.from_xml(self.pars_arg.xmlin, rundir=rundir, dummy=dummy)
    else:
      crank = manager.crank(rundir=rundir)
    # saving the input files/params as globals for crank so that they can be easily accessed there
    for strng in self.share_with_crank:
      if getattr(self.pars_arg,strng):
        if len(strng)>3 and strng[-2:]=='in' or strng[-3:]=='out':
          setattr(crank, strng, os.path.abspath(getattr(self.pars_arg,strng)))
        else:
          setattr(crank, strng, getattr(self.pars_arg,strng))
    if not self.pars_arg.xmlin:
      self.ParseInputFile(crank,dummy)
    return crank

  def PrintParams(self, proc):
    """Print all the parameters of the process and their description"""
    if proc=='processes':
      self.PrintProcesses()
    elif proc not in manager.crank.supported_procs:
      common.Warning('No such process {0}.'.format(proc))
      self.PrintProcesses()
      common.Error('No such process {0} supported.'.format(proc), debug=False)
    else:
      from process import process
      process.from_name(proc,None).PrintParams()

  def PrintProcesses(self):
    from process import process
    common.Info('The following processes are supported:')
    for proc in manager.crank.supported_procs:
      proc_obj = process.from_name(proc,None)
      common.Info('\n{0:15} {1}'.format(proc, proc_obj.name))
      if proc_obj.supported_procs:
        common.Info('{0:15}    supported subprocesses:  {1}'.format('', ', '.join(proc_obj.supported_procs)))
      if proc_obj.supported_progs:
        common.Info('{0:15}    supported programs:  {1}'.format('', ', '.join(proc_obj.supported_progs)))

  def PrintDataObjects(self):
    from data import data_container
    common.Info('The following data objects are supported:')
    for data_obj in data_container.__subclasses__():
      common.Info('\n{0:10} {1}'.format(data_obj.__name__, data_obj.description))
      common.Info('{0:10}    supported attributes:  {1}'.format( '', 
        ', '.join(data_obj.supported_attributes) ))
      if len(data_obj.types)>1:
        common.Info('{0:10}    possible types:  {1} (default: {2})'.format( '', 
          ', '.join(data_obj.types), data_obj.typ ))
      if data_obj.col_list:
        common.Info('{0:10}    available MTZ columns:  {1}'.format( '', ', '.join(data_obj.col_list) ))
      #common.Info('{0:10}    supported file types:  {1} (default: {2})'.format( '', 
      #  ', '.join(data_obj.supported_filetypes), data_obj.deffiletype ))


  @staticmethod
  def CheckFile(opt_str,fil):
    if not os.path.isfile(fil):
      common.Error('{0} file {1} does not exist.'.format(opt_str.strip('-').upper(),fil))
    if not os.path.isabs(fil):
      fil = os.path.join(os.getcwd(), fil)
  if no_argparse:
    class CheckFileAction():
      pass
    class OptionParserMod(optparse.OptionParser):
      def CheckFileAction(self,opt,fil):
        parse.CheckFile(opt,fil)
      def OpenFileAction(self, option, opt_str, value, parser):
        parse.CheckFile(opt_str,value)
        setattr(parser.values, option.dest, open(value))
      def add_argument(self,*args,**kwargs):
        if 'action' in kwargs and kwargs['action'] is parse.CheckFileAction:
          kwargs['action']="callback"
          kwargs['callback']=self.CheckFileAction
          kwargs['type']="str"
        if 'type' in kwargs and kwargs['type'] is file:
          kwargs['action']="callback"
          kwargs['callback']=self.OpenFileAction
          kwargs['type']="str"
        self.add_option(*args,**kwargs)
  else:
    class CheckFileAction(argparse.Action):
      def __call__(self, parser, namespace, values, option_string=None):
        parse.CheckFile(option_string,values)
        setattr(namespace, self.dest, values)

  def ParseCommandLine(self,argv=None):
    descr='CRANK2 software for automatic protein structure solution from experimental phases'
    epilog='All the *in and *out options are case insensitive and can be used without the initial dashes. '+\
           'If XMLIN or KEYIN is specified, it is used as input; otherwise standard input of keywords is assumed.'
    if no_argparse:
      parser = parse.OptionParserMod(description=descr, version=" ")
    else:
      parser = argparse.ArgumentParser(prog=program,description=descr,epilog=epilog)
      parser.add_argument('-v','-i','--version', action='version', version='')
    parser.add_argument('--xmlin', type=argparse.FileType('r'), help='specify input XML file to be executed')
    parser.add_argument('--keyin', type=argparse.FileType('r'), help='read keyword input from this file')
    parser.add_argument('--hklin', action=parse.CheckFileAction, help='specify input data file in MTZ format')
    parser.add_argument('--xyzin', action=parse.CheckFileAction, help='specify input coordinates file in PDB format')
    parser.add_argument('--seqin', action=parse.CheckFileAction, help='specify input sequence file')
    parser.add_argument('--xmlout', help='specify XML file name to which input will be written')
    parser.add_argument('--hklout', help='specify output MTZ file name')
    parser.add_argument('--xyzout', help='specify output PDB file name')
    parser.add_argument('--xyzsubout', help='specify output substructure PDB file name')
    parser.add_argument('--dirout', help='specify output directory for the CRANK job')
    parser.add_argument('--lgout', help='specify output loggraph file name')
    parser.add_argument('--logout', help='specify output log file name')
    parser.add_argument('-x','--gcx','--write-xml-only', action='store_true', help='write the crank2 XML file and exit')
    parser.add_argument('-d','--write-defaults', action='store_true', help='write presumable defaults in crank2 XML file and exit')
    parser.add_argument('-l','--list-params', dest='process', help='list parameters of the specified crank2 process and exit')
    parser.add_argument('-p','--list-processes', action='store_true', help='list all crank2 processes and exit')
    parser.add_argument('-o','--list-data-objects', action='store_true', help='list all crank2 data objects and exit')
    parser.add_argument('--disable-rvapi', action='store_true', help='disables the rvapi output')
    parser.add_argument('--rvapi-viewer', help='specify browser binary displaying rvapi output (0 to disable)')
    parser.add_argument('--rvapi-uri-prefix', help='specify rvapi uri prefix replacing dirout in filepaths (0 to use no prefix, ie relative paths)')
    parser.add_argument('--rvapi-no-tree', action='store_true', help='do not use tree widget for rvapi presentation')
    parser.add_argument('--rvapi-document', help='rvapi document for jsCoFE')
    parser.add_argument('--disable-mtz-label-prefix', action='store_true', help='disable program prefixes in output mtz labels')
    parser.add_argument('--graceful-preliminary-stop', action='store_true', help='quit gracefully (no exception) if pipeline stops preliminary due to no hope to succeed')
    # preprocessing of the passed arguments to make them case insensitive and the form without -- accepted
    if not no_argparse:
      inouts = parser.parse_args('').__dict__
    else:
      inouts, trash = parser.parse_args()
      inouts=inouts.__dict__
    if argv is None:
      argv=sys.argv[1:]
    argv2 = ['--'+a.strip('-').lower()  if a.strip('-').lower() in inouts else a  for a in argv]
    if not no_argparse:
      self.pars_arg = parser.parse_args(argv2)
    else:
      self.pars_arg, trash = parser.parse_args(argv2)
    # postprocessing
    if self.pars_arg.gcx or self.pars_arg.write_defaults:
      self.pars_arg.disable_rvapi = True
    self.pars_arg.rvapi_separate_steps = False


  class ReadlineIterator:
    """A file handle iterator that calls readline() to get its next value."""
    def __init__(self, f): self.f = f
    def __iter__(self): return self
    def __next__(self):
      line = self.f.readline()
      if line: 
        return line.strip()
      else:
        self.f.seek(0)
        raise StopIteration
    def next(self): return self.__next__()


  def ParseInputFile(self,crank,dummy=False):
    # may be changed to only str after py2 support is obsoleted
    self.unicod = str
    if sys.version_info[0] < 3:
      self.unicod = unicode
    # we can read from stdin or file
    if self.pars_arg.keyin:
      self.f_in=self.pars_arg.keyin
    else:
      if not self.f_in: # allow multiple reading of stdin on terminals that dont support seek in << redirection
        self.f_in=io.StringIO(self.unicod(sys.stdin.read()))
      if dummy:
        # a CCP4 CRANK header may come here
        print( 'CRANK2. Waiting for keyword input.' )
    # process the input
    row_prev=[]
    for row in csv.reader(parse.ReadlineIterator(self.f_in), dialect='crank_parse'):
      row=row_prev+row
      if row and row[-1]=='-':
        row_prev=row[:-1]
        continue
      else:
        row_prev=[]
      if row and not row[0].startswith('#'):
        if row[0].lower()=='end': 
          break
        row_orig=row[:]
        crank.InputElemInit(row,dummy)
        if row:
          common.Error('Wrong syntax or semantics: could not process keyword "{0}" in line "{1}"'.format(
            row[0], ' '.join(row_orig)))

