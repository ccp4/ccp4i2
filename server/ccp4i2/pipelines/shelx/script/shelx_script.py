from ccp4i2.core import CCP4ErrorHandling
from ccp4i2.pipelines.crank2.script import crank2_script


class shelx(crank2_script.crank2):
  TASKNAME                                  = 'shelx'

  # This task shares crank2's parameters; what makes it SHELX is running the
  # SHELXC/D/E route. That was set only by the React interface, on a pending
  # job, so i2run, the API and the tests ran crank2's own route under this
  # name. The def.xml now defaults to the SHELX route and process() insists
  # on it, since this task offers no other.
  def process(self, container=None):
      if container is None:
          self.container.inputData.SHELXCDE.set(True)
      return super(shelx, self).process(container)

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
