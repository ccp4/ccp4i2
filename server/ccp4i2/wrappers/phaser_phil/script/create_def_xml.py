#
#  Copyright (C) 2016 STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David Waterman
#  Acknowledgements: based on ideas and code by Nat Echols and Martin Noble.
#

"""Create phaser_phil.def.xml from PHIL parameters"""

import os
from lxml import etree

from ccp4i2.utils.phil_handlers import PhilTaskCreator


class PhaserTaskCreator(PhilTaskCreator):

  def __init__(self, debug=False):

    # Get the PHIL information as a scope
    import phaser
    from iotbx import phil
    phaser_path = phaser.__path__[0]
    with open(os.path.join(phaser_path, 'phenix_interface',
      '__init__.params'), 'r') as f:
      phaser_phil = phil.parse(f.read())

    PhilTaskCreator.__init__(self, phaser_phil, debug)
    self.fmt_dic['PLUGINNAME'] = "phaser_phil"

    self.inputDataXML = etree.fromstring('''
<container id="inputData">
    <content id="XYZIN">
        <className>CPdbDataFile</className>
        <qualifiers>
            <ifAtomSelection>True</ifAtomSelection>
            <allowUndefined>False</allowUndefined>
            <mustExist>True</mustExist>
            <requiredSubType>4,3,2,1,0</requiredSubType>
        </qualifiers>
    </content>
    <content id="F_SIGF">
        <className>CObsDataFile</className>
        <qualifiers>
            <allowUndefined>False</allowUndefined>
            <mustExist>True</mustExist>
        </qualifiers>
    </content>
</container>
''')

    self.outputDataXML = etree.fromstring('''
<container id="outputData">
    <content id="XYZOUT">
        <className>CPdbDataFile</className>
        <qualifiers>
            <default><subType>1</subType></default>
        </qualifiers>
    </content>
</container>
''')

  def __call__(self):

    if self.debug:
      print(etree.tostring(self.phil_tree, pretty_print=True))

    super(PhaserTaskCreator, self).__call__()

if __name__ == "__main__":

  ptc = PhaserTaskCreator()
  ptc()
