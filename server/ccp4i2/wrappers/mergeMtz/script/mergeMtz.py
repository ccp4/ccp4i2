import re

from ccp4i2.core.CCP4PluginScript import CPluginScript


class mergeMtz(CPluginScript):
    TASKNAME = 'mergeMtz'
    ERROR_CODES = {
        201: {'description': 'No input reflection files to merge'},
    }

    def startProcess(self):
      inFiles = []
      used = set()
      for n, miniMtz in enumerate(self.container.inputData.MINIMTZINLIST):
        if miniMtz.fileName.isSet() and miniMtz.fileName.exists():
          cls,contentFlag =  miniMtz.fileName.miniMtzType()
          if cls is not None:
            # The standard column names of this kind of mini-MTZ. (It called
            # cls().columnNames(True, contentFlag), the Qt-era signature: the
            # second argument is now asString, the flag was lost, and every
            # file gave no columns -- a merge of four files wrote H,K,L alone.)
            kind = cls()
            kind.contentFlag.set(contentFlag)
            inColumns = kind.columnNames(True)
            outColumns = inColumns
            if miniMtz.columnNames.isSet():
              userColumnNames = re.sub(' ','',miniMtz.columnNames.__str__())
              if userColumnNames.count(',') == inColumns.count(','):
                outColumns = userColumnNames
            # The tag is a prefix, as its tooltip says; without one, a name an
            # earlier file already used is prefixed with this file's position,
            # so two sets of map coefficients do not both become F,PHI.
            labels = [c for c in outColumns.split(',') if c]
            tag = str(miniMtz.columnTag).strip() if miniMtz.columnTag.isSet() else ''
            if tag:
              labels = ['%s_%s' % (tag, c) for c in labels]
            elif used.intersection(labels):
              labels = ['%d_%s' % (n + 1, c) for c in labels]
            used.update(labels)
            inFiles.append([miniMtz.fileName.__str__(), inColumns, ','.join(labels)])

      if not inFiles:
        # joinMtz with nothing to join succeeds and writes no file, so without
        # this the job reports success and produces no HKLOUT. A wrapper that
        # computes the inputs to its own step has to say when that computation
        # comes up empty.
        self.appendErrorReport(
            201,
            'No input reflection files to merge. %d were given; none had a '
            'file name that could be read.'
            % len(self.container.inputData.MINIMTZINLIST))
        return CPluginScript.FAILED

      rv = self.joinMtz(self.container.outputData.HKLOUT.fullPath.__str__(),inFiles)
      # Say what is in it (it was listed as "HKLOUT.mtz").
      self.container.outputData.HKLOUT.annotation = 'Merged from %d files: columns %s' % (
        len(inFiles), ','.join(f[2] for f in inFiles))
      return rv
