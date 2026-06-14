"""
     mtzheader.py: CCP4 GUI Project
     Copyright (C) 2011 STFC
     Author: Martyn Winn

     Reads the header of an MTZ reflection file (cell, spacegroup, resolution,
     reflection count, merged/unmerged) and maps it onto ccp4i2 data items.

     Originally wrapped Phil Evans' hklfile (via ccp4mg); reimplemented with
     gemmi so the wrapper imports and runs without the CCP4 native toolkits
     (keeping it usable on the slim, CCP4-free API). HKLIN is a CMtzDataFile,
     so gemmi.read_mtz_file covers every field the old hklfile path provided.
"""

import gemmi

from ccp4i2.core.CCP4PluginScript import CPluginScript


class mtzheader(CPluginScript):

    TASKNAME = 'mtzheader'

    def process(self):

       unsetData = self.checkInputData()
       if len(unsetData)>0:
         self.reportStatus(CPluginScript.FAILED)
         return

       # No output files, so skip checkOutputData
       mtz = gemmi.read_mtz_file(str(self.container.inputData.HKLIN))

       # An unmerged MTZ carries batch headers; a merged one has none.
       isMerged = len(mtz.batches) == 0
       if isMerged:
         print("Reflection file is merged")
       else:
         print("Reflection file is unmerged")

       cell = mtz.cell
       print("Cell = ", (cell.a, cell.b, cell.c, cell.alpha, cell.beta, cell.gamma))
       spacegroup_hm = mtz.spacegroup.hm if mtz.spacegroup is not None else ''
       print("Spacegroup = %s" % spacegroup_hm)
       resolution = mtz.resolution_high()
       print("Resolution = %f" % resolution)
       no_reflections = mtz.nreflections
       print("Number of reflections = %d" % no_reflections)

       # map gemmi values onto ccp4i2 data items
       self.container.outputData.MERGED = isMerged
       self.container.outputData.SPACEGROUPCELL.cell.set(
                                     { 'a' : cell.a,
                                       'b' : cell.b,
                                       'c' : cell.c,
                                       'alpha' : cell.alpha,
                                       'beta' : cell.beta,
                                       'gamma' : cell.gamma } )
       self.container.outputData.SPACEGROUPCELL.spaceGroup.set(spacegroup_hm)
       self.container.outputData.MAXRESOLUTION = resolution
       self.container.outputData.NO_REFLECTIONS = no_reflections

       # Needed!
       self.reportStatus(CPluginScript.SUCCEEDED)
