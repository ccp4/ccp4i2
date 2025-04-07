"""
Copyright (C) 2011 STFC
Author: Martyn Winn

This wrapper uses Stuart's python wrapper to Phil's hklfile.
It deliberately does not use mtzdump.
You need hklfile.py and _hklfile.so on PYTHONPATH. Latter is
set in ccp4i2/utils/setup.sh
"""

from ccp4mg import hklfile

from ....core.CCP4PluginScript import CPluginScript


class mtzheader(CPluginScript):

    TASKMODULE = 'test'      # Where this plugin will appear on the gui
    TASKTITLE = 'Read MTZ header' # A short title for gui menu
    TASKNAME = 'mtzheader'   # Task name - should be same as class name
    TASKVERSION= 0.0               # Version of this plugin

    def process(self):

       unsetData = self.checkInputData()
       if len(unsetData)>0:
         self.reportStatus(CPluginScript.FAILED)
         return

       # No output files, so skip checkOutputData

       reflection_list = hklfile.ReflectionList()
       reflection_list.init(str(self.container.inputData.HKLIN))

       ftype = reflection_list.FileType()
       if ftype.FileType() == hklfile.ReflectionFileType.UNKNOWN:
         print("HKLIN is of unknown type")
         return
       else:
         print("Opening reflection file of type %s" % ftype.FileType())

       isMerged = reflection_list.Merged()
       if isMerged:
         print("Reflection file is merged")
       else:
         print("Reflection file is unmerged")

       # issue of which cell to take
       cell = hklfile.Scell(reflection_list.Cell())
       print("Cell = ", cell.UnitCell())
       spacegroup = hklfile.SpaceGroup(reflection_list.Spacegroup())
       print("Spacegroup = %s" % spacegroup.Symbol_hm())
       resolution = reflection_list.Resolution().ResHigh()
       print("Resolution = %f" % resolution)
       no_reflections = reflection_list.NumberReflections()
       print("Number of reflections = %d" % no_reflections)

       # map hklfile objects onto ccp4i2 data items
       self.container.outputData.MERGED = isMerged
       self.container.outputData.SPACEGROUPCELL.cell.set(
                                     { 'a' : cell.UnitCell()[0],
                                       'b' : cell.UnitCell()[1],
                                       'c' : cell.UnitCell()[2],
                                       'alpha' :cell.UnitCell()[3],
                                       'beta' : cell.UnitCell()[4],
                                       'gamma' : cell.UnitCell()[5] } )
       self.container.outputData.SPACEGROUPCELL.spaceGroup.set(spacegroup.Symbol_hm())
       self.container.outputData.MAXRESOLUTION = resolution
       self.container.outputData.NO_REFLECTIONS = no_reflections

       # Needed!
       self.reportStatus(CPluginScript.SUCCEEDED)
