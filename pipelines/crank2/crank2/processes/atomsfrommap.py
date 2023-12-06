#!/usr/bin/python
import os,sys,re
from process import process
import common,data
par=common.parameter

class atomsfrommap(process):
  name="obtaining new atoms from (difference) map"
  short_name="atom picking"
  supported_progs=["peakmax"]
  supported_params = {}
  supported_params['rms_threshold'] = par( desc='RMS threshold for atom picking from map', typ=(float,bool) )
  supported_params['too_close'] = par( desc='distance (in A) - new atoms closer to an existing (substr.) atom won''t be added', typ=(float,bool) )
  supported_params['max_new_atoms'] = par( desc='maximal number of new atoms added', typ=(int,bool) )
  supported_params['bfactor'] = par( desc='set B factor of the new atoms to this value', typ=(float,bool) )
  supported_params['occupancy'] = par( desc='set occupancy of the new atoms to this value', typ=(float,bool) )

  #def RunPreprocess(self,*args,**kwargs):
  #  process.process.RunPreprocess(self,*args,**kwargs)

  def TreatInOutPar(self, set_all_par=False):
    # assuming peakmax as of now but the code is general and other programs can be easily added
    self.prog = self.GetOrAddProg('peakmax')
    # peakmax does not know the model type thus it is not propagated - to be set/fixed in the postprocessing
    self.prog.out.never_propagate=True
    # changing the program defaults - max_new_atoms,too_close,rms_threshold,bfactor,occupancy 
    # set to False will use the original defaults
    if not self.IsParam('max_new_atoms'):
      if self.inp.Get('model',typ=('substr','partial+substr'),has_num_atoms=True):
        self.SetParam('max_new_atoms', int(max(2,self.inp.Get('model',typ=('substr','partial+substr'),has_num_atoms=True).exp_num_atoms//2)))
      else:
        self.SetParam('max_new_atoms', 10)
    if not self.IsParam('too_close'):
      self.SetParam('too_close', 3.5)
      subs=self.inp.Get('model',typ=('substr','partial+substr'),has_atomtypes=True)
      seq_obj=self.inp.Get('sequence', filetype=data.sequence.supported_filetypes, try_convert=False)
      if seq_obj and not seq_obj.seqstr:  seq_obj.GetSequenceString()
      if subs and subs.GetAtomType() and subs.GetAtomType()=='S' and (not seq_obj or 'C' in seq_obj.seqstr):
        self.SetParam('too_close', 1.5)
    if not self.IsParam('rms_threshold'):
      self.SetParam('rms_threshold', 5.0)
    if not self.IsParam('occupancy'):
      self.SetParam('occupancy', 0.1)
    if not self.IsParam('bfactor'):
      gcx=self.AddProg('gcx', propagate_out=False)
      try:
        gcx.Run()
      except (gcx.ProgramRunError,common.CrankError):
        if hasattr(sys,'exc_clear'): sys.exc_clear()
        wilsonB = 50.
      else:
        wilsonB = gcx.GetStat('wilson_B')
      self.programs.remove(gcx)
      self.SetParam('bfactor', wilsonB)
    # just to supress warning for case of no atoms found.
    self.GetParam('too_close')
    process.TreatInOutPar(self,set_all_par)

  def RunBody(self,*args,**kwargs):
    self.prog.added_atoms=0
    # use anomalous diff maps by preference
    adm=self.prog.inp.Get('mapcoef',typ='anom-diff')
    if adm and adm is not self.prog.inp.Get('mapcoef'):
      self.prog.inp.Add(adm)
    try:
      self.prog.Run()
    except self.prog.ProgramRunError as e:
      if 'Threshold too high' in str(e):
        self.prog.out.Clear('model')
        info = 'No new atoms found in the map (no peaks above the threshold).'
        self.Info(info)
        self.LGInfo(info)
        if self.rv_report is not None:
          self.rv_report.Text('&emsp;'+info, flush=True)
      else:
        try: #python2
          raise
        except RuntimeError: #python3
          raise e #from None
    else:
      mapt,err='',''
      if self.prog.inpmap.GetType()=='anom-diff':
        mapt='anomalous difference'
        err='Please check that the data contains anomalous signal and the model contains anomalously scattering atoms with non-zero occupancies.'
      # this can happen if the (anomalous) diff map if flat/broken
      if self.prog.GetStat('numpeaks') and not self.prog.GetStat('heights'):
        common.Error('{0} failed finding atoms. {1}'.format(self.prog.name,err))
      info = '{0} peaks picked in the {1} density map.'.format( self.prog.GetStat('numpeaks'), mapt )
      self.prog.added_atoms=self.prog.GetStat('numpeaks')
      self.Info(info)
      self.Info('Peak heights (in sigma): {0}'.format( ', '.join(map(str,self.prog.GetStat('heights'))) ))
      self.LGInfo(info)
      if self.rv_report is not None:
        self.rv_report.Text('&emsp;'+info, flush=True)


  def RemoveClose(self,atomtypes,inp_pdb_obj):
    # remove atoms too close to existing atoms of the same type
    if self.GetParam('too_close') and atomtypes and inp_pdb_obj:
      disu=False
      if 'S' in atomtypes.keys():
        seq_obj=self.inp.Get('sequence', filetype=data.sequence.supported_filetypes, try_convert=False)
        if (not seq_obj or 'C' in seq_obj.GetSequenceString()):
          disu=True
      deleted_list=[0,0,0]
      try:
        import gemmi
        struct=gemmi.read_structure(self.out.Get('model').GetFileName())
        struct_inp=gemmi.read_structure(inp_pdb_obj.GetFileName())
        for ch1 in struct[0]:
          totres1=len(ch1)
          for i,r1 in enumerate(reversed(ch1)): # we need to loop excplicitly and backwards to remove residues with the current gemmi
            for s2 in struct_inp[0].all():
              if (struct.cell.find_nearest_image(r1[0].pos,s2.atom.pos).dist()<=self.GetParam('too_close')):
                del ch1[totres1-i-1]
                deleted_list[0]+=1
                break
        if disu:
          disulist=[]
          for i,s1 in enumerate(struct_inp[0].all()):
            for j,s2 in enumerate(struct_inp[0].all()):
              if j>i and s1.atom.element.name=='S' and s2.atom.element.name=='S' and struct.cell.find_nearest_image(s1.atom.pos,s2.atom.pos).dist()<=2.75:
                if s1.atom.pos not in disulist:  disulist.append(s1.atom.pos)
                if s2.atom.pos not in disulist:  disulist.append(s2.atom.pos)
          struct=gemmi.read_structure(self.out.Get('model').GetFileName())
          for ch1 in struct[0]:
            totres1=len(ch1)
            for i,r1 in enumerate(reversed(ch1)): # we need to loop excplicitly and backwards to remove residues with the current gemmi
              for dpos in disulist:
                if (struct.cell.find_nearest_image(r1[0].pos,dpos).dist()<=3): # remove peaks close to disuphides
                  del ch1[totres1-i-1]
                  deleted_list[1]+=1
                  break
              else:
                for s1 in struct[0].all():
                  if r1[0].pos!=s1.atom.pos and struct.cell.find_nearest_image(r1[0].pos,s1.atom.pos).dist()<=5: # remove peaks close to each other
                    del ch1[totres1-i-1]
                    deleted_list[2]+=1
                    break
        if any(deleted_list):
          delfilename=self.out.Get('model').GetFileName()[:-4]+'_del.pdb'
          struct.write_pdb(delfilename)
          self.out.AddFileToChild(self.out.Get('model'),delfilename,'pdb')
        #print deleted_list
      except ImportError:
      # deprecated.  kept as fallback if gemmi not available.
          # first generate the entire unit cell - pdbcur's deldist does not take symmetry into account
          pdbcur=self.AddProg('pdbcur',propagate_out=False)
          pdbcur.inp.Set(inp_pdb_obj)
          pdbcur.SetKey('genunit')
          pdbcur.SetKey('symop','X,Y,Z'),pdbcur.SetKey('symop','X+1,Y+1,Z+1'),pdbcur.SetKey('symop','X-1,Y-1,Z-1')
          pdbcur.SetKey('symop','X+1,Y,Z'),pdbcur.SetKey('symop','X,Y+1,Z'),pdbcur.SetKey('symop','X,Y,Z+1')
          pdbcur.SetKey('symop','X+1,Y+1,Z'),pdbcur.SetKey('symop','X+1,Y,Z+1'),pdbcur.SetKey('symop','X,Y+1,Z+1')
          pdbcur.SetKey('symop','X-1,Y,Z'),pdbcur.SetKey('symop','X,Y-1,Z'),pdbcur.SetKey('symop','X,Y,Z-1')
          pdbcur.SetKey('symop','X-1,Y-1,Z'),pdbcur.SetKey('symop','X-1,Y,Z-1'),pdbcur.SetKey('symop','X,Y-1,Z-1')
          pdbcur.SetKey('symop','X+1,Y-1,Z'),pdbcur.SetKey('symop','X+1,Y,Z-1'),pdbcur.SetKey('symop',' X+1,Y-1,Z-1')
          pdbcur.SetKey('symop','X+1,Y-1,Z+1'),pdbcur.SetKey('symop','X+1,Y+1,Z-1'),pdbcur.SetKey('symop','X,Y+1,Z-1'),pdbcur.SetKey('symop','X,Y-1,Z+1')
          pdbcur.SetKey('symop','X-1,Y+1,Z'),pdbcur.SetKey('symop','X-1,Y,Z+1'),pdbcur.SetKey('symop',' X-1,Y+1,Z+1')
          pdbcur.SetKey('symop','X-1,Y+1,Z-1'),pdbcur.SetKey('symop','X-1,Y-1,Z+1'),pdbcur.SetKey('symop','X,Y-1,Z+1'),pdbcur.SetKey('symop','X,Y+1,Z-1')
          pdbcur.SetKey('symcommit')
          pdbcur.Run()
          # now remove
          pdbcur.ClearAnyParams()
          pdbcur.out.never_propagate=False
          pdbcur.inp.Set(self.out.Get('model'))
          #with open(inp_pdb_obj.GetFileName()) as f:
          atom_regexp = r'((ATOM  )|(HETATM)).{6}(.{4}).{14}\s*(-?\d+.?\d*)\s*(-?\d+.?\d*)\s*(-?\d+.?\d*)'
          with open(pdbcur.out.Get('model').GetFileName()) as f:
            for line in f:
              re_res=re.match(atom_regexp,line)
              if re_res and re_res.group(4).strip() in atomtypes:
                pdbcur.SetKey('deldist',(re_res.group(5),re_res.group(6),re_res.group(7),self.GetParam('too_close')))
          pdbcur.Run()
          deleted_list=pdbcur.GetStat('atoms_deleted')
          if disu:
            pdbcur.SetRunDir(os.path.join(pdbcur.rundir,'S1'))
            pdbcur.out.never_propagate=True
            pdbcur.inp.Set(inp_pdb_obj)
            lstdel=[]
            with open(inp_pdb_obj.GetFileName()) as f:
              for line in f:
                re_res=re.match(atom_regexp,line)
                if re_res and re_res.group(4).strip() in atomtypes:
                  pdbcur.ClearAnyParams()
                  pdbcur.SetKey('deldist',(re_res.group(5),re_res.group(6),re_res.group(7),2.5))
                  pdbcur.Run()
                  deleted_list2=pdbcur.GetStat('atoms_deleted')
                  #print("dellist2",deleted_list2)
                  if deleted_list2 and deleted_list2[0]>=2: # mark peaks close to disuphides for deletion
                    lstdel.append((re_res.group(5),re_res.group(6),re_res.group(7),3.0))
            with open(self.out.Get('model').GetFileName()) as f:
              for line in f:
                re_res=re.match(atom_regexp,line)
                if re_res:
                  pdbcur.inp.Set(self.out.Get('model'))
                  pdbcur.ClearAnyParams()
                  pdbcur.SetKey('lvdist',(re_res.group(5),re_res.group(6),re_res.group(7),5.0))
                  pdbcur.Run()
                  deleted_list2=pdbcur.GetStat('atoms_deleted_left')
                  #print("dellist2",deleted_list2)
                  if (deleted_list2 and deleted_list2[0][1]>=2) or (not deleted_list2 and self.prog.GetStat('numpeaks')-sum(deleted_list)==2): # mark peaks close to each other for deletion
                    first=True
                    with open(pdbcur.out.Get('model').GetFileName()) as g:
                      for line in g:
                        re_res2=re.match(atom_regexp,line)
                        if re_res2:
                          if not first:
                            lstdel.append((re_res2.group(5),re_res2.group(6),re_res2.group(7),0.1))
                          first=False
            pdbcur.inp.Set(self.out.Get('model')) #delete marked
            pdbcur.out.never_propagate=False
            pdbcur.ClearAnyParams()
            for delkey in lstdel:
              pdbcur.SetKey('deldist',delkey)
            pdbcur.SetRunDir(os.path.join(pdbcur.rundir,'S2'))
            pdbcur.Run()
            #print("dellist3",pdbcur.GetStat('atoms_deleted'))
            deleted_list.extend(pdbcur.GetStat('atoms_deleted'))

      if any(deleted_list):
        self.prog.added_atoms=self.prog.GetStat('numpeaks')-sum(deleted_list)
        info = '{0} of the peaks too close to existing atoms, thus {1} atoms added.'.format( 
          sum(deleted_list), self.prog.added_atoms )
        self.Info(info)
        self.LGInfo(info)
        if self.rv_report is not None:
          self.rv_report.Text('&emsp;'+info, flush=True)


  def RunPostprocess(self,restore=True,*args,**kwargs):
    # in case of substructure, merge and "fix" the resulting substructure pdb file
    inp_pdb_obj=self.inp.Get('model',typ=('substr','partial+substr'),filetype='pdb')
    if not inp_pdb_obj:
      inp_substr_nopdb=self.inp.Get('model',typ=('substr','partial+substr'),has_atomtypes=True)
    if not inp_pdb_obj:
      inp_substr_nopdb=self.inp.Get('model',has_atomtypes=True)
    if inp_pdb_obj or inp_substr_nopdb:
      if self.prog.out.Get('model'):
        self.out.Add(self.prog.out.Get('model'),propagate=False)
        # retrieve atomtypes
        if inp_pdb_obj:
          atomtypes, atomtype1 = inp_pdb_obj.GetAtomTypes(), inp_pdb_obj.GetAtomType()
        else:
          atomtypes, atomtype1 = inp_substr_nopdb.GetAtomTypes(), inp_substr_nopdb.GetAtomType()
        if not atomtypes: 
          at_obj=self.inp.Get('model',has_atomtypes=True)
          if at_obj:
            atomtypes,atomtype1=at_obj.GetAtomTypes(),at_obj.GetAtomType()
        self.RemoveClose(atomtypes,inp_pdb_obj)
        # fix atom and residue names
        fixsubstr=self.AddProcess('fixsubstrpdb')
        fixsubstr.inp.Set(self.out.Get('model'))
        if not fixsubstr.inp.Get('model',has_atomtypes=True):
          fixsubstr.inp.Get('model').SetAtomTypes(atomtypes,atomtype1)
        if fixsubstr.inp.Get('model').GetType()=='unknown':
          fixsubstr.inp.Get('model').SetType('substr')
        fixsubstr.Run()
        # merge pdbs
        if inp_pdb_obj:
          mergepdb=self.AddProg('pdbmerge',propagate_inp=False)
          mergepdb.SetKey('nomerge')
          mergepdb.inp.Add(self.out.Get('model'))
          mergepdb.inp.Add(inp_pdb_obj)
          mergepdb.Run()
        elif inp_substr_nopdb.GetType()=='substr':
          self.inp.SetFileToChild(inp_substr_nopdb,fixsubstr.out.Get('model').GetFileName('pdb'))
          self.out.Add(inp_substr_nopdb)
        # clean up so that the "unknown" model is not propagated to the next step from the output
        self.out.Delete(self.prog.out.Get('model'),propagate=False)
      elif inp_pdb_obj:
        self.out.Add(inp_pdb_obj)
    elif self.prog.inpmap.GetType()=='anom-diff':
      common.Warning('Type of found anomalous atoms not known. They will be ignored.')
    process.RunPostprocess(self,restore,*args,**kwargs)
