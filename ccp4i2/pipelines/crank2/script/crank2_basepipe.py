class crank2_basepipe():
  # separation of the pipeline steps setting as the script does not work with gui objects anymore

  def __init__(self):
    # definition of pipelines.
    self.shelx_steps = [ "substrdet", "phdmmb", "building", "ref" ]
    self.crank2_steps = [ "substrdet", "phas", "handdet", "dmfull", "building", "ref" ]
    self.rebuild_steps = [ "refatompick", ] + self.crank2_steps[2:]
    #self.all_steps = [ "refatompick", "substrdet", "phdmmb", "phas", "handdet", "dmfull", "building", "ref" ]
    self.base_steps = self.crank2_steps[:]
    self.BaseStepsInd()

  def SetBaseSteps(self,container):
    self.container=container
    inp=container.inputData
    self.base_steps = self.shelx_steps[:]  if inp.SHELXCDE  else self.crank2_steps[:]
    if inp.INPUT_PARTIAL:
      self.base_steps = self.rebuild_steps[:]
    if str(inp.EXPTYPE)=='SAD' and self.base_steps[1]=='phas':
      self.base_steps[1]='refatompick'
      print('adjusted to refatompick')
    #elif str(inp.EXPTYPE)!='SAD' and not inp.INPUT_PARTIAL and self.base_steps[1]=='refatompick':
    #  self.base_steps[1]='phas'
    #  print 'adjusted to phas'
    self.BaseStepsInd()

  def BaseStepsInd(self):
    self.base_steps_ind = dict( (item,i) for i,item in enumerate(self.base_steps) )

  def CheckStartEnd(self, step):
    return step in self.base_steps and \
           self.base_steps_ind[self.container.inputData.START_PIPELINE] <= self.base_steps_ind[step] and \
           self.base_steps_ind[self.container.inputData.END_PIPELINE] >= self.base_steps_ind[step]

  def ToggleDetection(self):
    return self.CheckStartEnd('substrdet')
  def ToggleNotDetection(self):
    return not self.ToggleDetection()

  def ToggleShelxCDE(self):
    return self.CheckStartEnd('phdmmb')
    #return ((self.container.inputData.END_PIPELINE  != 'substrdet') and
    #        self.container.inputData.SHELXCDE and not self.container.inputData.INPUT_PARTIAL)

  def TogglePeakSearch(self):
    return self.CheckStartEnd('refatompick')

  def TogglePhasing(self):
    toggle = False
    if (self.container.inputData.INPUT_PARTIAL):
      if (self.container.inputData.PARTIAL_AS_SUBSTR):
        toggle = True
    else:
      if self.CheckStartEnd('phas') and not self.container.inputData.SHELXCDE:
        toggle = True
    return toggle

  def ToggleHandDetermination(self):
    toggle = False
    # the user can turn hand det off separately - this is taken into account by ToggleHandDeterminationDoHand!
    if (self.container.inputData.INPUT_PARTIAL):
      if (self.container.inputData.PARTIAL_AS_SUBSTR):
        toggle = True
    else:
      if ( self.CheckStartEnd('handdet') and not self.container.inputData.SHELXCDE ): 
        toggle = True
    return toggle

  def ToggleDensityModification(self):
    toggle     = False
    if (self.container.inputData.INPUT_PARTIAL):
      if (self.container.inputData.PARTIAL_AS_SUBSTR):
        toggle = True
    else:
      if ( self.CheckStartEnd('dmfull') and not self.container.inputData.SHELXCDE):
        toggle = True
    return toggle

  def ToggleModelBuilding(self):
    return self.CheckStartEnd('building')

  def ToggleUseComb(self):
    return self.CheckStartEnd('building') and self.container.controlParameters.USE_COMB and \
           str(self.container.controlParameters.MB_PROGRAM)!='arpwarp'
    #or self.container.inputData.INPUT_PARTIAL ) 

  def ToggleNCS(self):
    return str(self.container.controlParameters.COMB_PHDMMB_DMFULL_DM_PROGRAM)=='parrot' and \
           self.container.inputData.MONOMERS_ASYM and int(self.container.inputData.MONOMERS_ASYM)>1

  def ToggleRefine(self):
    return self.CheckStartEnd('ref')
