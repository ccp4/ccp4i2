
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

    # -- a file chosen fills the table (any route: the page, the API, i2run) --

    SOURCES = ("ASUCONTENTIN", "XYZIN", "SEQIN")

    def afterParameterSet(self, obj):
        """A sequence, model or AU-contents file was set: add its sequences to
        ASU_CONTENT, as the page used to in the browser. A sequence already in
        the table gains copies rather than a second row; rows with no sequence
        give way to the file's. Returns what was done, for the caller to say.

        The page did this itself after setting the file, so a job set up any
        other way (an agent, i2run) had the file and an empty table, and
        agents gave up there (gamma, CDK1/cyclin B/CKS2, 2026-10-13)."""
        import re
        from ccp4i2.lib.utils.files.digest import digest_file_object
        from ccp4i2.lib.utils.sequences.asu_entries import entries_from_digest
        if obj.objectName() not in self.SOURCES or obj.parent() is not self.container.inputData:
            return None
        if not obj.isSet():
            return None  # cleared, not chosen
        label = (str(obj.annotation) if obj.annotation.isSet() else "") or str(obj.baseName) or "file"
        digest = digest_file_object(obj)
        if not digest or digest.get("status") == "Failed":
            return f"Could not read {label}" + (f": {digest.get('reason')}" if digest else "")
        entries, chains = entries_from_digest(digest, label)
        entries = [e for e in entries if e.get("polymerType") in ("PROTEIN", "DNA", "RNA")]
        if not entries:
            return f"No polymer sequences found in {label}"

        items = self.container.inputData.ASU_CONTENT
        key = lambda polymer, sequence: f"{str(polymer).upper()}:{re.sub(r'[^A-Za-z]', '', str(sequence))}"
        for i in reversed(range(len(items))):  # rows with no sequence give way
            if not str(items[i].sequence).strip():
                items.pop(i)
        present = {key(item.polymerType, item.sequence): item for item in items}
        added = folded = 0
        for entry in entries:
            hit = present.get(key(entry["polymerType"], entry["sequence"]))
            copies = int(entry.get("nCopies") or 1)
            if hit is not None:
                hit.nCopies.set(int(hit.nCopies) + copies)
                folded += 1
                continue
            item = items.makeItem()
            items.append(item)
            item.sequence.set(re.sub(r"\s", "", str(entry["sequence"])))
            item.polymerType.set(entry["polymerType"])
            item.name.set(re.sub(r"[^A-Za-z0-9_]", "_", str(entry.get("name") or "sequence")))
            if entry.get("description"):
                item.description.set(str(entry["description"]))
            item.nCopies.set(copies)
            present[key(entry["polymerType"], entry["sequence"])] = item
            added += 1
        parts = [f"{chains} chain{'' if chains == 1 else 's'} from {label}"]
        if added:
            parts.append(f"{added} new row{'' if added == 1 else 's'}")
        if folded:
            parts.append(f"{folded} added as copies of rows already present")
        return ": ".join(parts)

    def fillFromFile(self, name):
        """Fill the table from ``name`` (ASUCONTENTIN, XYZIN or SEQIN) as it is
        set now, and save: for a file that reached the job without being set
        through the API (carried over from an earlier job), which the page
        fills once when it opens. Only a pending job."""
        from ccp4i2.lib.utils.jobs.editing import editable_job, save
        job, why = editable_job(self)
        if job is None:
            return {"success": False, "error": why}
        if name not in self.SOURCES:
            return {"success": False, "error": f"{name} is not one of {', '.join(self.SOURCES)}"}
        note = self.afterParameterSet(getattr(self.container.inputData, name))
        save(self, job)
        return {"success": True, "note": note}

    def validity(self):
      """The contents come from the list or, if it is empty, from a file.

      The interface copies a loaded file into the list; a job set up without
      the interface may give the file alone, so an empty list is not an
      error when there is a file. With neither, the list's own error stands.
      """
      from ccp4i2.core import CCP4ErrorHandling

      error = super(ProvideAsuContents, self).validity()
      # A protein entry this short is rarely the crystallised construct: an
      # agent typed "MNNI" for both BetaBlip components (2026-10-13). A
      # warning, not an error: peptides are real.
      for i, item in enumerate(self.container.inputData.ASU_CONTENT):
          seq = "".join(str(item.sequence).split()) if item.sequence.isSet() else ""
          polymer = str(item.polymerType) if item.polymerType.isSet() else "PROTEIN"
          if seq and polymer == "PROTEIN" and len(seq) < 20:
              error.append(klass=self.TASKNAME, code=210,
                           details=f'Entry {i + 1} has only {len(seq)} residues: is this the whole '
                                   'sequence? Upload the sequence file into its source, or fetch it.',
                           name=f'{self.TASKNAME}.container.inputData.ASU_CONTENT[{i}].sequence',
                           severity=CCP4ErrorHandling.SEVERITY_WARNING)
      # The same chain entered twice (two routes - the page's SEQIN and an
      # entry's source - each add it): overpacked, and Crank2's own Matthews
      # then found no fit (gamma SAD, Haiku, 2026-10-13)
      seen = {}
      for i, item in enumerate(self.container.inputData.ASU_CONTENT):
          seq = "".join(str(item.sequence).split()).upper() if item.sequence.isSet() else ""
          if not seq:
              continue
          if seq in seen:
              error.append(klass=self.TASKNAME, code=211,
                           details=f'Entries {seen[seq] + 1} and {i + 1} have the same sequence: '
                                   'if it is one molecule, keep one entry and set its copies.',
                           name=f'{self.TASKNAME}.container.inputData.ASU_CONTENT[{i}].sequence',
                           severity=CCP4ErrorHandling.SEVERITY_WARNING)
          else:
              seen[seq] = i
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
