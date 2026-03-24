from ccp4i2.core import CCP4ErrorHandling
from ccp4i2.pipelines.crank2.script import crank2_script


class shelx(crank2_script.crank2):
  TASKNAME                                  = 'shelx'

  def validity(self):
      error = super(shelx, self).validity()
      # ATOM_TYPE is required
      atom_type = getattr(self.container.inputData, 'ATOM_TYPE', None)
      if atom_type is not None:
          val = str(atom_type).strip() if atom_type.isSet() else ''
          if not val:
              error.append(
                  klass=self.TASKNAME, code=200,
                  details='Heavy atom element type is required',
                  name=f'{self.TASKNAME}.container.inputData.ATOM_TYPE',
                  severity=CCP4ErrorHandling.SEVERITY_ERROR,
              )
      # SEQIN is recommended
      seqin = getattr(self.container.inputData, 'SEQIN', None)
      if seqin is not None and not seqin.isSet():
          error.append(
              klass=self.TASKNAME, code=201,
              details='Providing sequence information is strongly recommended',
              name=f'{self.TASKNAME}.container.inputData.SEQIN',
              severity=CCP4ErrorHandling.SEVERITY_WARNING,
          )
      return error
