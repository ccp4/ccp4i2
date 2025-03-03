#!/usr/bin/python
import os,sys,multiprocessing
from process import process
import common


class mb(process):
  name="real space model building"
  short_name="model tracing"
  supported_progs=["buccaneer","shelxe"]
  supported_params={}
  supported_params['from_weighted_map'] = common.parameter(desc='Build from weigted maps (if available; Buccaneer only)',typ=bool)
  supported_params['from_mr_model'] = common.parameter(desc='Build from MR model (if available; Buccaneer only)',typ=bool)
  supported_params['low_res_par'] = common.parameter(desc='Use low resolution parameters in model building.', typ=bool )

  def TreatInOutPar(self, set_all_par=False):
    if not self.GetProg(supported=True):
      self.AddProg('buccaneer')
    self.resol=False
    bucc = self.GetProg('buccaneer')
    # changing quite a few buccaneer defaults here
    if bucc:
      # a very naive check whether anisotropy should be applied...
      abc=self.inp.Get('fsigf').GetCell(self)[:3]
      #if min(abc)*2.5<max(abc) and max(abc)>100. and not bucc.IsKeyOrArg('anisotropy-correction'):
      if not bucc.IsKeyOrArg('anisotropy-correction'):
        bucc.SetKey('anisotropy-correction',True)
      fsf = self.inp.Get('fsigf',typ='average',col='f',is_native=True)
      if not fsf:
        fsf = self.inp.Get('fsigf',typ='average',col='f')
      self.resol=fsf.GetResolution(self) if fsf else None
      if self.resol and self.resol>=3.0 and self.GetParam('low_res_par') is None:
        self.SetParam('low_res_par')
      if self.GetParam('low_res_par'):
        if not bucc.IsKey('sequence-reliability'):
          if self.resol and self.resol>4.0:
            bucc.SetKey('sequence-reliability',1.)
          else:
            bucc.SetKey('sequence-reliability',0.99)
        if not bucc.IsKeyOrArg('model-filter'):
          bucc.SetKey('model-filter')
        if not bucc.IsKeyOrArg('model-filter-sigma'):
          bucc.SetKey('model-filter-sigma',2.0)
        # removing the filtering step for very low res - disabled as not really helping for 3din
        #if self.resol and self.resol>4.0 and not bucc.IsKey('filter'):
        #  if not bucc.GetKey('find'):
        #    bucc.SetKey('find'),bucc.SetKey('grow'),bucc.SetKey('join'),bucc.SetKey('link')
        #    bucc.SetKey('sequence'),bucc.SetKey('correct'),bucc.SetKey('ncsbuild'),bucc.SetKey('prune')
        #    bucc.SetKey('rebuild'),bucc.SetKey('tidy') # should test >=1.6.3?
        #  else:
        #    bucc.SetKey('filter',False)
      if not bucc.IsKeyOrArg('resolution'):
        bucc.SetKey('resolution',2.0)
      # disabling free omission by default:  the idea is that calculated F's should be inputted
      if not bucc.IsKeyOrArg('colin-free'):
        bucc.SetKey('colin-free',False)
      if not bucc.IsKeyOrArg('jobs'):
        bucc.SetKey('jobs', int(os.getenv('OMP_NUM_THREADS',multiprocessing.cpu_count())) )
    process.TreatInOutPar(self,set_all_par)

  def RunBody(self,*args,**kwargs):
    prog=self.GetProg(supported=True)
    # run and catch if buccanner builds nothing and crashes
    try:
      prog.Run()
    except prog.ProgramRunError as e:
      if prog.nick=='buccaneer' and 'failed with error code -11' in str(e) and prog.GetStat('res_built',accept_none=True)==0:
        common.Warning(str(e))
        if hasattr(sys,'exc_clear'): sys.exc_clear()
        common.Error('Map too noisy, model could not be traced by {0}.'.format(prog.name), nosuccess=True)
      else:
        raise
    if prog.nick=='shelxe' and prog.GetStat('untraceable'):
      p = next( (p for p in self.GetProgs(supported=True) if p.nick!='shelxe'), None)
      if p:
        common.Warning('{0} could not trace the noisy map, will try {1}.'.format(prog.name,p.name))
        prog.out.ClearAll()
        p_new=self.AddProgCopy( self.programs.pop( self.programs.index(p) ), ind=0, deeper_copy=True )
        p_new.Run()
      else:
        common.Error('{0} could not trace the map, too noisy.'.format(prog.name))

  def RunPostprocess(self,restore=True,*args,**kwargs):
    prog=self.GetProg(supported=True)
    self.Info('')
    num_frag = prog.GetStat('frag_built') if prog.nick=='buccaneer' else len(list(filter(None,prog.GetStat('res_built_per_frag'))))
    self.Info("{0} residues in {1} fragments built.".format(prog.GetStat('res_built'), num_frag))
    if prog.nick=='buccaneer':
      self.Info("Completness of uniquely located residues against sequence {0}%, against built {1}%.".format(
        prog.GetStat('compl_chain'), prog.GetStat('compl_res')))
      self.Info('')
      # temporary and bad workaround for the problem of hybrid36 incompatibility between buccaneer and refmac
      # all the non-numerical seqID residues from the built models are removed!
      h36=False
      new_cont=[]
      out_pdb=self.out.Get('model',filetype='pdb').GetFileName()
      with open(out_pdb) as f:
        for line in f:
          if (line.startswith('ATOM') or line.startswith('HETATM')) and \
             len(line)>60 and line[22]!=' ' and line[22] not in [str(i) for i in range(10)]:
            h36=True
          else:
            new_cont.append(line)
      if h36:
        common.Warning('Hybrid36 extension residues found in the built PDB. These residues will be removed as Refmac does not support hybrid36.')
        new_pdb=out_pdb[:-4]+'_h36_removed.pdb'
        with open(new_pdb,'w') as g:
          for line in new_cont:
            g.write(line)
        self.out.AddFileToChild(self.out.Get('model',filetype='pdb'), new_pdb)
    elif prog.nick=='shelxe':
      prog.FixOutput()
    process.RunPostprocess(self,restore,*args,**kwargs)
