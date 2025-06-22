import os

from phaser.command_line import main

from ....core import CCP4Container
from ....core.CCP4PluginScript import CPluginScript


class phaser_phil(CPluginScript):
    TASKNAME = 'phaser_phil'   # Task name - should be same as class name and match pluginTitle in the .def.xml file
    TASKVERSION= 0.1               # Version of this plugin
    MAINTAINER = 'kyle.stevenson@stfc.ac.uk;david.waterman@stfc.ac.uk'
    ERROR_CODES = { 201 : {'description' : 'Failed to analyse output files' },
                    202 : {'description' : 'Failed applying selection ot PDB file' }
                    }
    PURGESEARCHLIST = [ [ 'hklin.mtz' , 0 ],
                       ['log_mtzjoin.txt', 0]
                       ]
    TASKCOMMAND='ccp4-python'

    # Copied here from xia2_dials - if useful should be moved to the base
    # class
    def extract_parameters(self, container):
        """Walk through a container locating parameters that have been set
        and return a list of name, value pairs"""

        result=[]
        dataOrder = container.dataOrder()
        contents = [getattr(container,name) for name in dataOrder]
        for model in contents:
            if isinstance(model, CCP4Container.CContainer):
                result.extend(self.extract_parameters(model))
            elif model.isSet():
                name = model.objectName().replace('__','.')
                # ensure commas are converted to whitespace-separated lists. Only
                # whitespace appears to work correctly with PHIL multiple
                # choice definitions.
                val = str(model.get()).split()
                val = ' '.join([v[:-1] if v.endswith(',') else v for v in val])
                result.append((name, val))
        return result

    def makeCommandAndScript(self,**kw):
        par = self.container.controlParameters

        # PHIL parameters set by the gui
        phil_file = os.path.normpath(os.path.join(
                    self.getWorkDirectory(), 'phaser.phil'))
        with open(phil_file, 'w') as f:
            for (name, val) in self.extract_parameters(par):
                #self.appendCommandLine()
                f.write(name + '={0}\n'.format(val))

        # FIXME it is possibly better, rather than writing out a file of
        # definition=value lines, to fetch this against the Phaser master phil
        # and hence construct a fully formatted phil definition for the job.

        # Get path to phaser's main.py
        self.appendCommandLine([main.__file__])

        self.appendCommandLine([phil_file])

        return CPluginScript.SUCCEEDED
