
from ccp4i2.core.CCP4PluginScript import CPluginScript
import shutil
from lxml import etree

class ProvideAsuContents(CPluginScript):
    TASKNAME = 'ProvideAsuContents'

    # -- UniProt (reached through the object_method endpoint) ---------------

    def uniprotCandidates(self, text, organism=None, limit=10):
        """UniProt entries for a protein as named; see CSequence.uniprotCandidates."""
        from ccp4i2.core.CCP4ModelData import CSequence
        return CSequence.uniprotCandidates(text, organism, limit)

    def fetchUniProt(self, accession, residue_range=None, index=None, nCopies=None):
        """Fill ASU_CONTENT item ``index`` (a new item when None) from a
        UniProt entry, cut to ``residue_range`` for the construct, with
        ``nCopies`` if given; and save. Only a pending job."""
        from ccp4i2.lib.utils.jobs.editing import editable_job, save
        from ccp4i2.lib.utils.sequences import uniprot
        job, why = editable_job(self)
        if job is None:
            return {"success": False, "error": why}
        items = self.container.inputData.ASU_CONTENT
        if index is not None and not 0 <= int(index) < len(items):
            return {"success": False, "error": f"there is no ASU_CONTENT[{index}]"}
        try:  # fetched first, so a failure changes nothing
            entry = uniprot.fetch(accession, residue_range)
        except uniprot.UniProtError as err:
            return {"success": False, "error": str(err)}
        if index is None:
            items.append(items.makeItem())
            index = len(items) - 1
        index = int(index)
        added = items[index].fillFromEntry(entry)
        if nCopies is not None:
            items[index].nCopies.set(int(nCopies))
        save(self, job)
        return {"success": True, "index": index, "added": added}

    def validity(self):
      """The contents come from the list or, if it is empty, from a file.

      The interface copies a loaded file into the list; a job set up without
      the interface may give the file alone, so an empty list is not an
      error when there is a file. With neither, the list's own error stands.
      """
      from ccp4i2.core import CCP4ErrorHandling

      error = super(ProvideAsuContents, self).validity()
      if not self.container.inputData.ASUCONTENTIN.isSet():
          return error
      filtered = CCP4ErrorHandling.CErrorReport()
      for err in error.getErrors():
          if err.get('code') == 101 and str(err.get('name', '')).endswith('inputData.ASU_CONTENT'):
              continue
          filtered.append(err.get('class', ''), err.get('code', 0), err.get('details', ''),
                          err.get('name', ''), err.get('severity', 0))
      return filtered

    def startProcess(self):
      inp = self.container.inputData
      if len(inp.ASU_CONTENT) == 0 and inp.ASUCONTENTIN.isSet():
          # The interface copies a loaded file's contents into the list to be
          # edited; a job set up without it (i2run, the API, an agent) records
          # the file's contents as they are, rather than nothing.
          inp.ASUCONTENTIN.loadFile()
          inp.ASU_CONTENT.set(inp.ASUCONTENTIN.fileContent.seqList)
      asuFileObject = self.container.outputData.ASUCONTENTFILE
      asuFileObject.fileContent.seqList.set(self.container.inputData.ASU_CONTENT)
      asuFileObject.saveFile()

      xmlroot = etree.Element('ASUCONTENTMATTHEWS')
      totWeight = 0.0
      if len(self.container.inputData.ASU_CONTENT) > 0:
          entries = etree.SubElement(xmlroot,"entries")
          polymerMode = ""
          for seqObj in self.container.inputData.ASU_CONTENT:
              if seqObj.nCopies > 0:
                  if seqObj.polymerType == "PROTEIN":
                      if polymerMode == "D":
                          polymerMode = "C"
                      elif polymerMode == "":
                          polymerMode = "P"
                  if seqObj.polymerType in ["DNA","RNA"]:
                      if polymerMode == "P":
                          polymerMode = "C"
                      elif polymerMode == "":
                          polymerMode = "D"
              totWeight = totWeight + seqObj.molecularWeight(seqObj.polymerType)
              entry = etree.SubElement(entries,"entry")
              nCopies = etree.SubElement(entry,"copies")
              name = etree.SubElement(entry,"name")
              weight = etree.SubElement(entry,"weight")
              sequence = etree.SubElement(entry,"sequence")
              nCopies.text = str(seqObj.nCopies)
              etree.SubElement(entry, "polymerType").text = str(seqObj.polymerType)
              name.text = str(seqObj.name)
              weight.text = "{0:.1f}".format(float(seqObj.molecularWeight(seqObj.polymerType)))
              sequence.text = str(seqObj.sequence)
          totalWeightTag = etree.SubElement(xmlroot,"totalWeight")
          totalWeightTag.text = str(totWeight)
          # P protein, D nucleic acid, C both: it sets the density, and so
          # the range of solvent content a crystal can have
          if polymerMode:
              etree.SubElement(xmlroot, "polymerMode").text = polymerMode

      if self.container.inputData.HKLIN.isSet() and len(self.container.inputData.ASU_CONTENT) > 0:
          if totWeight > 1e-6:
              rv = self.container.inputData.HKLIN.fileContent.matthewsCoeff(molWt=totWeight,polymerMode=polymerMode)
              vol = rv.get('cell_volume','Unkown')
              volumeTag = etree.SubElement(xmlroot,"cellVolume")
              volumeTag.text = str(vol)
              matthewsComposition = etree.SubElement(xmlroot,"matthewsCompositions")
              for result in rv.get('results',[]):
                  comp = etree.SubElement(matthewsComposition,"composition")
                  nMolecules = etree.SubElement(comp,"nMolecules")
                  solventPercentage = etree.SubElement(comp,"solventPercentage")
                  matthewsCoeff = etree.SubElement(comp,"matthewsCoeff")
                  matthewsProbability = etree.SubElement(comp,"matthewsProbability")
                  nMolecules.text = str(result.get('nmol_in_asu'))
                  solventPercentage.text = "{0:.2f}".format(float(result.get('percent_solvent')))
                  matthewsCoeff.text = "{0:.2f}".format(float(result.get('matth_coef')))
                  matthewsProbability.text = "{0:.2f}".format(float(result.get('prob_matth')))

      newXml = etree.tostring(xmlroot,pretty_print=True)
      with open (self.makeFileName('PROGRAMXML')+'_tmp','w') as programXmlFile:
          programXmlFile.write(newXml.decode("utf-8"))
      shutil.move(self.makeFileName('PROGRAMXML')+'_tmp', self.makeFileName('PROGRAMXML'))

      return CPluginScript.SUCCEEDED

    def processOutputFiles(self):
      parts = []
      for seqObj in self.container.inputData.ASU_CONTENT:
          copies = int(seqObj.nCopies) if seqObj.nCopies else 1
          name = str(seqObj.name) if seqObj.name else 'unknown'
          ptype = str(seqObj.polymerType) if seqObj.polymerType else 'PROTEIN'
          seq = str(seqObj.sequence).replace(' ', '').replace('\n', '') if seqObj.sequence else ''
          resCount = len(seq)
          if copies > 1:
              parts.append('{0}x {1} ({2}, {3} res)'.format(copies, name, ptype.lower(), resCount))
          else:
              parts.append('{0} ({1}, {2} res)'.format(name, ptype.lower(), resCount))
      if parts:
          self.container.outputData.ASUCONTENTFILE.annotation = ', '.join(parts)
      from ccp4i2.core.CCP4ErrorHandling import CErrorReport
      return CErrorReport()
