#!/usr/bin/python
import os,sys
from program import program
import common

class parrot(program):
  name="Parrot"
  binary="cparrot"
  labelout_prefix="PARR_"
  modif='-'
  always_merge=True
  stat={}
  stat['ncs_operator'] = common.stats(name='NCS operator', multiple=True,
            regexp=r"-ncs-operator\s+(\S+,\S+,\S+,\S+,\S+,\S+,\S+,\S+,\S+)")
  stat['ncs_operator_correl'] = common.stats(name='NCS operator correlation', regexp=r'NCS operator statistics:\s+Operator_number  Mask_volume\/ASU  Correlation'+r'(?:\s+\d+\s+\S+\s+(\S+)\s)?'*99)
  stat['radius_auto'] = common.stats(name='Automatically determined filtering radius for solvent flattening',
            regexp=r"Suggested radius for solvent mask determination:\s*(\S+)")
  references = ( "Cowtan K (2010) Recent developments in classical density modification. Acta Cryst. D66, 470-478.", )

  def Init(self):
    self.SetArg('stdin')

  def TreatInput(self):
    # F,sigF
    #self.fsf = self.inp.Get('fsigf',typ='average',col='f',is_native=True)
    #if not self.fsf:
    self.fsf = self.inp.Get('fsigf',typ='average',col='f')
    if self.fsf:
      self.SetKey('mtzin-wrk', self.fsf.GetFileName('mtz'))
      self.SetKey('colin-wrk-fo', self.fsf.GetFullLabel('f','sigf'))
    else:
      common.Error('No data inputted for {0}'.format(self.name))
    # Parrot seems to always require phase distribution (even if not asked for phase combination and map is inputted...)
    for mapc in self.inp.GetAll('mapcoef',typ=('best','combined')):
      if mapc.GetLabel('hla'):
        self.SetKey('colin-hl', mapc.GetFullLabel('hla','hlb','hlc','hld'))
        break
      if mapc.GetLabel('ph') and mapc.GetLabel('fom'):
        self.SetKey('colin-phifom', mapc.GetFullLabel('ph','fom'))
        break
    else:
      common.Error('No phase distribution inputted for {0}'.format(self.name))
    # input map (optional - phase distr. and Fo used if not inputted)
    mapc=self.inp.Get('mapcoef',typ=('weighted','mask'),col=('f','ph'))
    if mapc:
      if 'dmmsk' in mapc.custom:
        self.SetKey('colin-wrk-fmsk', mapc.GetFullLabel('f','ph'))
        #self.SetKey('mapin', '../segmentmap/out.ccp4')
        mapc=self.inp.Get('mapcoef',typ='weighted',col=('f','ph'),exclude_cont=mapc) # used for unt models feedback recycling
        if mapc:
          self.SetKey('colin-wrk-fc', mapc.GetFullLabel('f','ph'))
      else:
        self.SetKey('colin-wrk-fc', mapc.GetFullLabel('f','ph'))
    # solvent content - if supplied with an input container
    solv_obj = self.inp.Get(has_solvent_content=True)
    if solv_obj and not self.IsKey('solvent-content'):
      self.SetKey('solvent-content', solv_obj.solvent_content)
    if self.process.IsNonFalseVirtPar('solvent_perturb') and self.GetKey('solvent-content'):
      perturb=self.process.GetVirtPar('solvent_perturb')
      if self.GetKey('solvent-content')>0.8:  perturb=min(perturb,0.0)
      if self.GetKey('solvent-content')<0.2:  perturb=max(perturb,0.0)
      self.SetKey('solvent-content', self.GetKey('solvent-content') + perturb, keep_previous=False )
    # detect ncs
    if self.process.GetParam('ncs_det') is not False:
      mr_model=self.inp.Get('model',typ=('partial','partial+substr'),custom='ncs',filetype='pdb')
      heavy_model=self.inp.Get('model',typ='substr',filetype='pdb')
      if mr_model and self.process.GetParam('ncs_det_mr') is not False and not self.process.GetParam('ncs_det_ha'):
        self.SetKey('pdbin-wrk-mr',mr_model.GetFileName('pdb'))
      elif heavy_model and self.process.GetParam('ncs_det_ha') is not False:
        self.SetKey('pdbin-wrk-ha',heavy_model.GetFileName('pdb'))
      if self.process.GetParam('ncs_det') and not heavy_model and not mr_model:
        common.Error('Asked for NCS detection by {0} but no appropriate model to detect from was supplied.'.format(self.name))

  def TreatParams(self):
    # real space DM
    if self.process.nick=='dm':
      self.SetKey('modify-map-and-terminate')
    if self.process.IsNonFalseVirtPar('cycles'):
      self.SetKey('cycles', self.process.GetVirtPar('cycles'))
    if self.process.IsNonFalseVirtPar('solvent_content'):
      self.SetKey('solvent-content', self.process.GetVirtPar('solvent_content'))
    if self.process.GetParam('solventmask_radius'):
      self.SetKey('solvent-mask-filter-radius', self.process.GetVirtPar('solventmask_radius'), keep_previous=False)
    program.TreatParams(self)

  def DefineOutput(self):
    out_mapc = self.out.AddNew( 'mapcoef', self.nick+'.mtz' )
    if self.process.nick=='dm':
      out_mapc.SetType('densmod')
      out_mapc.SetLabel( ['f','ph'] )
    else:
      out_mapc.SetType('combined')
      out_mapc.SetLabel( ['hla','hlb','hlc','hld'] )
    out_mapc.SetCrystalName(self.fsf.GetCrystalName())
    out_mapc.SetDataName(self.fsf.GetDataName())

  def TreatOutput(self):
    if self.process.nick=='dm':
      self.SetKey( 'colout-fc', self.out.mapcoef[-1].GetFullLabel('f','ph') )
    else:
      self.SetKey( 'colout-hl', self.out.mapcoef[-1].GetFullLabel('hla','hlb','hlc','hld') )
    self.SetKey( 'mtzout', self.out.mapcoef[-1].GetFileName('mtz') )
