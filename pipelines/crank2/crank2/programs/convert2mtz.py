#!/usr/bin/python
import os,sys
from program import program
import common

class convert2mtz(program):
  name="CONVERT2MTZ"
  binary="convert2mtz"
  modif='-'
  # mapping between the filetype ('inputformat' parameter) and convert2mtz 'input' key value
  in_key = { 'hkl': 'shelx', 'sca': 'scal', 'phs': "user '(3I4,2F7.1,I4)'" }

  def Init(self):
    self.SetArg('stdin')

  def TreatInput(self):
    self.inpcont=self.inp.GetAll(filetype=self.in_key.keys())
    if self.inpcont:
      self.SetArg('hklin', self.inpcont[0].GetFileName(self.process.GetParam('inputformat')))
      for cont in self.inpcont:
        # intensities get priority over f's
#        if cont.GetType()=='plus' and cont.GetLabel('i'):
#          self.AddToKey('colin', ('I(+)='+cont.GetLabel('i'), 'SIGI(+)='+cont.GetLabel('sigi')))
#        elif cont.GetType()=='plus' and cont.GetLabel('f'):
#          self.AddToKey('colin', ('F(+)='+cont.GetLabel('f'), 'SIGF(+)='+cont.GetLabel('sigf')))
#        elif cont.GetType()=='minus' and cont.GetLabel('i'):
#          self.AddToKey('colin', ('I(-)='+cont.GetLabel('i'), 'SIGI(-)='+cont.GetLabel('sigi')))
#        elif cont.GetType()=='minus' and cont.GetLabel('f'):
#          self.AddToKey('colin', ('F(-)='+cont.GetLabel('f'), 'SIGF(-)='+cont.GetLabel('sigf')))
#        elif cont.GetType()=='average' and cont.GetLabel('i'):
#          self.AddToKey('colin', ('I='+cont.GetLabel('i'), 'SIGI='+cont.GetLabel('sigi')))
#        elif cont.GetType()=='average':
#          self.AddToKey('colin', ('F', 'SIGF'))
        if cont.GetType()=='fa':
          if cont.GetLabel('f'):
            self.AddToKey('colin', (cont.GetLabel('f'), cont.GetLabel('sigf')))
            # shelx fa's include estimation of phases
            if self.process.GetParam('inputformat')=='hkl':
              self.AddToKey('colin', 'PA')
          elif cont.GetLabel('e'):
            self.AddToKey('colin', (cont.GetLabel('f'), cont.GetLabel('sigf')))
    else:
      Error('No mtz file inputted to {0}'.format(self.name))

  def TreatParams(self):
    if not self.process.IsParam('inputformat') and not self.IsKey('input'):
      common.Error('input file type not specified for {0}'.format(self.name))
    #elif not self.IsKey('input'):
    #  self.SetKey('input', self.in_key[self.process.GetParam('inputformat')])
    #else:
    #  if self.GetKey('input') in in_key.itervalues:
    #    of = [of for of,o in in_key if o==self.GetKey('input')][0]
    #    self.process.SetParam('inputformat', of)
    #  else:
    #    common.Warning('Input "{0}" not predefined for {1}'.format(self.GetKey('input'),self.name))
    program.TreatParams(self)

  def DefineOutput(self):
    for cont in self.inpcont:
      outmtz=self.out.Add(cont)
      self.out.AddFileToChild(outmtz, self.nick+'.mtz', 'mtz')
      if cont.GetType()=='fa':
        if cont.GetLabel('f'):
          outmtz.SetLabel(('f','sigf'))
        elif cont.GetLabel('e'):
          outmtz.SetLabel(('e','sige'))
        if self.process.GetParam('inputformat')=='hkl':
          outmtz.SetLabel('alpha','PA')

  def TreatOutput(self):
    self.SetArg('mtzout', self.out.Get().GetFileName('mtz'))
