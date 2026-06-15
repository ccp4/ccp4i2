import os

import gemmi
import numpy as np
from lxml import etree

from ccp4i2.core import CCP4ErrorHandling, CCP4Utils
from ccp4i2.core.CCP4PluginScript import CPluginScript

# Fixed seed so a re-run reproduces the same free set -- the reproducibility
# property of CCP4 freerflag (whose default also uses iseed=123456789). The
# native implementation reproduces freerflag's *semantics* (binning into
# 1/fraction segments, free set = flag 0, one draw per symmetry-equivalence
# class so equivalents share a flag) rather than its exact byte stream, which
# bottoms out in the compiler's RNG and is gfortran-version-dependent.
FREER_SEED = 123456789


def assign_class_flags(hkl, spacegroup, irfrac, seed=FREER_SEED):
    """Free-R flags reproducing freerflag's semantics.

    One random draw per symmetry/Friedel equivalence class (so equivalent
    reflections share a flag), values in [0, irfrac-1] with the free set = 0,
    reproducible for a fixed seed.

    Returns (flags array of length len(hkl), list of per-reflection class keys).
    """
    gops = spacegroup.operations()
    asu = gemmi.ReciprocalAsu(spacegroup)
    keys = [tuple(asu.to_asu((int(h[0]), int(h[1]), int(h[2])), gops)[0]) for h in hkl]
    rng = np.random.default_rng(seed)
    class_flag = {k: int(rng.integers(0, irfrac)) for k in sorted(set(keys))}
    flags = np.array([class_flag[k] for k in keys], dtype=int)
    return flags, keys


class freerflag(CPluginScript):
    TASKNAME = 'freerflag'
    # gemmi/numpy-native: no CCP4 binary (see generateFreeR). Removing TASKCOMMAND
    # is what lets this run inline on a CCP4-free server (ccp4_free=True).

    # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    def processInputFiles(self):
        print("freerflag wrapper, processInputFiles")

        # We have two possible resolution cutoffs
        #  1) global cutoff if RESMAX is set
        #  2) if COMPLETE and CUTRESOLUTION and FreeR is higher resolution
        #     than the data, then FreeR is cut to the data resolution
        self.freerCutoff = 0.0   # will be set >0 if Freer cutoff is done
        self.globalCutoff = 0.0
       
        self.highResD = float(self.container.inputData.F_SIGF.fileContent.resolutionRange.high)
        self.highResF = 0.0
        if self.container.inputData.FREERFLAG.isSet():
            self.highResF = float(self.container.inputData.FREERFLAG.fileContent.resolutionRange.high)

        if self.container.controlParameters.GEN_MODE == 'COMPLETE':
            if self.container.controlParameters.CUTRESOLUTION:
                #  cut the resolution of the FreeR set if it is higher than the data
                self.cutResolution()
        
            self.hklin,error = self.makeHklin(['F_SIGF','FREERFLAG'])
            print('freerflag.processInputFiles',self.hklin,error)
            if error.maxSeverity()>CCP4ErrorHandling.SEVERITY_WARNING:
                return CPluginScript.FAILED
            else:
                # Optional global cutoff
                if self.container.controlParameters.RESMAX:
                    resmax = self.container.controlParameters.RESMAX.get()
                    self.globalResolutionCutoff(resmax)

                return CPluginScript.SUCCEEDED

        else:
            self.hklin = self.container.inputData.F_SIGF.__str__()

        return CPluginScript.SUCCEEDED

    # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    def cutResolution(self):

        if self.highResD <= self.highResF:
            # Data go to higher resolution, no need to do anything
            return

        # File names
        FSIGF_file = str(self.container.inputData.F_SIGF)
        FREER_file = str(self.container.inputData.FREERFLAG)

        # Cut resolution of FreeR file to match data file
        outfile = os.path.join(self.workDirectory, 'FREEOUT.mtz')        
        highRes= self.highResD-0.001
        print("* Cutting Freer to ", highRes)
        mtz = gemmi.read_mtz_file(FREER_file)
        mtz.set_data(mtz.array[mtz.make_d_array() >= highRes])
        mtz.write_to_file(outfile)
        self.container.inputData.FREERFLAG.setFullPath(outfile)
        self.freerCutoff = highRes
    # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    def globalResolutionCutoff(self, resmax):
        # after mtzjoin, if needed. Joined file is in self.hklin
        mtz = gemmi.read_mtz_file(self.hklin)
        highRes = mtz.resolution_high()
        if highRes < resmax-0.001:
            # yes cut the data
            print(">>**>> cutting all data", highRes, resmax)
            mtz.set_data(mtz.array[mtz.make_d_array() >= resmax])
            mtz.write_to_file(self.hklin)
            self.globalCutoff = resmax

    # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    def makeCommandAndScript(self):
      # Native task: no external command. Just fix the output path and validate.
      self.hklout = os.path.join(self.workDirectory, "hklout.mtz")
      if self.container.controlParameters.GEN_MODE == 'COMPLETE':
          if not self.container.inputData.FREERFLAG.isSet():
              self.appendErrorReport(101)
      return 0

    # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    def startProcess(self, command=None, **kwargs):
      # Generate the free-R flags in-process with gemmi/numpy (no CCP4 binary).
      try:
          self.generateFreeR()
      except Exception as e:
          import traceback
          traceback.print_exc()
          self.appendErrorReport(102, str(e))
          return CPluginScript.FAILED
      return CPluginScript.SUCCEEDED

    # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    def postProcessCheck(self, processId=None):
      # No subprocess ran; success is simply that hklout was written.
      if os.path.isfile(self.hklout):
          return self.SUCCEEDED, self.SUCCEEDED, 0
      return self.FAILED, self.FAILED, 1

    # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    def _irfrac(self):
      cp = self.container.controlParameters
      rfrac = float(cp.FRAC) if cp.FRAC.isSet() else 0.05
      if not (0.0 < rfrac < 1.0):
          rfrac = 0.05
      return int(round(1.0 / rfrac))

    def _class_keyer(self, spacegroup):
      gops = spacegroup.operations()
      asu = gemmi.ReciprocalAsu(spacegroup)

      def key(h):
          # ASU representative merges symmetry- and Friedel-equivalent reflections,
          # so they all land in the same class (and thus share a flag).
          rep, _isym = asu.to_asu((int(h[0]), int(h[1]), int(h[2])), gops)
          return (rep[0], rep[1], rep[2])
      return key

    def generateFreeR(self):
      cp = self.container.controlParameters
      complete = (str(cp.GEN_MODE) == 'COMPLETE')
      uniqueify = bool(cp.UNIQUEIFY)

      mtz = gemmi.read_mtz_file(self.hklin)
      spacegroup = mtz.spacegroup or gemmi.SpaceGroup('P 1')
      key = self._class_keyer(spacegroup)
      rng = np.random.default_rng(FREER_SEED)

      data = np.array(mtz, copy=True)
      hkl = data[:, :3].astype(int)

      if complete:
          # Preserve existing flags, fill only the missing ones (per class).
          idx = self._freer_column_index(mtz)
          existing = data[:, idx].astype(float)
          present = ~np.isnan(existing)
          irfrac = (int(round(existing[present].max())) + 1) if present.any() else self._irfrac()
          keys = [key(h) for h in hkl]
          class_flag = {}
          # seed class flags from reflections that already have one
          for i, kk in enumerate(keys):
              if present[i]:
                  class_flag.setdefault(kk, int(existing[i]))
          # draw flags for the remaining classes
          for kk in sorted(set(keys)):
              class_flag.setdefault(kk, int(rng.integers(0, irfrac)))
          for i, kk in enumerate(keys):
              if not present[i]:
                  existing[i] = class_flag[kk]
          data[:, idx] = existing
          mtz.set_data(data)
      else:
          irfrac = self._irfrac()
          if uniqueify:
              hkl, data = self._complete_to_unique(mtz, data)
          flags, _keys = assign_class_flags(hkl, spacegroup, irfrac, seed=FREER_SEED)
          mtz.add_column('FreeR_flag', 'I')
          data = np.concatenate([data, flags.reshape(-1, 1).astype(np.float32)], axis=1)
          mtz.set_data(data)

      mtz.write_to_file(self.hklout)

    def _freer_column_index(self, mtz):
      for j, c in enumerate(mtz.columns):
          if c.type == 'I' and ('free' in c.label.lower() or 'flag' in c.label.lower()):
              return j
      raise RuntimeError("COMPLETE mode: no FreeR (type I) column found in the joined input")

    def _complete_to_unique(self, mtz, data):
      # UNIQUEIFY: extend the reflection list to the full unique (ASU) set within
      # the resolution range; generated reflections carry NaN data + a flag.
      cp = self.container.controlParameters
      dmin = mtz.resolution_high()
      if cp.RESMAX.isSet() and float(cp.RESMAX) > 0.0:
          dmin = float(cp.RESMAX)
      full = np.array(gemmi.make_miller_array(mtz.cell, mtz.spacegroup, dmin), dtype=int)
      ncol = data.shape[1]
      obs = {(int(h[0]), int(h[1]), int(h[2])): data[i] for i, h in enumerate(data[:, :3].astype(int))}
      rows = np.empty((len(full), ncol), dtype=np.float32)
      for i, h in enumerate(full):
          t = (int(h[0]), int(h[1]), int(h[2]))
          if t in obs:
              rows[i] = obs[t]
          else:
              rows[i] = np.nan
              rows[i, :3] = h
      return full, rows

    # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    def processOutputFiles(self):
      print("freerflag - processOutputFiles")
      self.xmlRoot = etree.Element('FREERFLAGINFO')
      annotation = ''

      if self.container.controlParameters.GEN_MODE == 'COMPLETE':
         annotation = 'Extended freeR;'
         self.addElement(self.xmlRoot, 'Mode', 'Complete')
         if self.container.controlParameters.CUTRESOLUTION:
             self.addElement(self.xmlRoot, 'CutFreerResolution', 'True')
         else:
             self.addElement(self.xmlRoot, 'CutFreerResolution', 'False')

         self.addElement(self.xmlRoot, 'ObservedDataResolution',
                         '{:6.2f}'.format(self.highResD))
         self.addElement(self.xmlRoot, 'FreeR_Resolution',
                         '{:6.2f}'.format(self.highResF))

         if self.freerCutoff > 0.0:
             self.addElement(self.xmlRoot, 'FreerCutResolution',
                             '{:6.2f}'.format(self.freerCutoff))
      else:
          annotation = 'New freeR'
          self.addElement(self.xmlRoot, 'Mode', 'New')
          fraction = 0.05  # default
          if self.container.controlParameters.FRAC.isSet():
              fraction = self.container.controlParameters.FRAC
          self.addElement(self.xmlRoot, 'Fraction', str(fraction))

      if self.container.controlParameters.UNIQUEIFY:
         self.addElement(self.xmlRoot, 'Unique', 'True')

      if self.container.controlParameters.RESMAX:
          resmax = self.globalCutoff
          if resmax > 0.0:
              self.addElement(self.xmlRoot, 'GlobalResolutionLimit',
                              '{:6.2f}'.format(resmax))

      with open ( self.makeFileName('PROGRAMXML'),'w' ) as xmlFile:
         xmlString = etree.tostring (self.xmlRoot, pretty_print=True )
         CCP4Utils.writeXML(xmlFile,xmlString)

      if self.container.controlParameters.GEN_MODE == 'COMPLETE':
          error = self.splitHklout(['FREEROUT'],['FREER'])
      else:          
          error = self.splitHklout(['FREEROUT'],['FreeR_flag'])

      if error.maxSeverity()>CCP4ErrorHandling.SEVERITY_WARNING:
        return CPluginScript.FAILED
     
      # Annotation cf aimless_pipe
      if self.container.outputData.FREEROUT.fileContent.spaceGroup.isSet():
          sgname = self.container.outputData.FREEROUT.fileContent.spaceGroup.__str__()
      else:
          sgname = 'Unk'

      highresFRformatted = "%7.2f" % float(self.container.outputData.FREEROUT.fileContent.resolutionRange.high)
      title =' Spg:'+str(sgname).strip()+';Resln:'+highresFRformatted.strip() + "A;"
      try:
          title = title + "Cell:"+self.container.outputData.FREEROUT.fileContent.cell.guiLabel()
      except Exception as e:
          print('Error writing cell parameters',e)

      annotation += title
      print("Annotation", annotation)
      self.container.outputData.FREEROUT.annotation.set(annotation)

      return CPluginScript.SUCCEEDED

    # - - - - - - - - -  - - - - - - - - -  - - - - - - - - - 
    def getXML(self):
        try:
            return self.xmlRoot
        except:
            return None
    # - - - - - - - - -  - - - - - - - - -  - - - - - - - - - 
    def addElement(self, containerXML, elementname, elementtext):
        #print 'addElement', elementname, type(elementtext), elementtext 
        e2 = etree.Element(elementname)
        e2.text = elementtext
        containerXML.append(e2)
