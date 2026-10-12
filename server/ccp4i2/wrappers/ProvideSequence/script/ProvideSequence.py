import os
import tempfile
from io import StringIO

from lxml import etree

from ccp4i2.core import CCP4ModelData, CCP4Utils
from ccp4i2.core.CCP4PluginScript import CPluginScript


class ProvideSequence(CPluginScript):

    TASKNAME = 'ProvideSequence'

    def previewSequences(self):
        """Parse SEQUENCETEXT and return a preview of the interpreted sequences.

        Called via the object_method endpoint so the frontend can show
        a preview table before the job is run.
        """
        from ccp4i2.wrappers.ProvideAlignment.script.ProvideAlignment import importAlignment
        import tempfile

        text = str(self.container.controlParameters.SEQUENCETEXT)
        if not text.strip():
            return {"sequences": [], "format": None, "commentary": "No sequence text provided."}

        tempFile = tempfile.NamedTemporaryFile(suffix='.txt', delete=False)
        tempFile.file.write(text.encode('utf-8'))
        tempFile.close()

        try:
            alignment, fmt, commentary = importAlignment(tempFile.name)
        finally:
            os.remove(tempFile.name)

        if alignment is None:
            return {
                "sequences": [],
                "format": None,
                "commentary": commentary.getvalue(),
            }

        sequences = []
        for seq in alignment:
            sequences.append({
                "id": seq.id,
                "name": seq.name,
                "description": seq.description,
                "sequence": str(seq.seq),
                "length": len(seq.seq),
            })

        return {
            "sequences": sequences,
            "format": fmt,
            "commentary": commentary.getvalue(),
        }

    # -- UniProt (reached through the object_method endpoint) ---------------

    def uniprotCandidates(self, text, organism=None, limit=10):
        """UniProt entries for a protein as named ("Human CDK2", "cyclin D
        from human", an accession): how the text was read and the candidates,
        best first; none is chosen."""
        return CCP4ModelData.CSequence.uniprotCandidates(text, organism, limit)

    def fetchUniProt(self, accession, residue_range=None, append=True):
        """Add a UniProt entry's sequence to SEQUENCETEXT as a FASTA record
        (its header names the accession, organism and any residue range, so
        the output sequence carries where it came from), cut to
        ``residue_range`` ("175-432") for the crystallised construct; and
        save. Only a pending job; ``append`` False replaces the text."""
        from ccp4i2.lib.utils.jobs.editing import editable_job, save
        from ccp4i2.lib.utils.sequences import uniprot
        job, why = editable_job(self)
        if job is None:
            return {"success": False, "error": why}
        try:
            entry = uniprot.fetch(accession, residue_range)
        except uniprot.UniProtError as err:
            return {"success": False, "error": str(err)}
        text = str(self.container.controlParameters.SEQUENCETEXT) if append else ""
        text = (text.rstrip() + "\n" if text.strip() else "") + entry["fasta"]
        self.container.controlParameters.SEQUENCETEXT.set(text)
        save(self, job)
        return {"success": True, "added": {k: entry[k] for k in (
            "accession", "entry_name", "protein_name", "gene", "organism", "reviewed", "range")}
            | {"length": len(entry["sequence"])}}

    # -- the text from a file, as the interface fills it -----------------

    def _files_given(self):
        inp = self.container.inputData
        return [f for f in (inp.SEQIN, inp.XYZIN) if f.isSet()]

    def _text_from_files(self):
        """What the interface puts in SEQUENCETEXT when a file is chosen:
        a sequence file's content, or a model's chain sequences as FASTA.
        The interface does it in the browser; a job set up through the API
        (an agent, i2run) never had it done, was told SEQUENCETEXT was
        required, and agents typed in "beta.seq content" or "\\n" to pass
        (2026-10-12)."""
        inp = self.container.inputData
        if inp.SEQIN.isSet() and inp.SEQIN.exists():
            with open(inp.SEQIN.getFullPath(), encoding="utf-8", errors="replace") as handle:
                return handle.read()
        if inp.XYZIN.isSet() and inp.XYZIN.exists():
            from ccp4i2.lib.utils.files.digest import digest_file_object
            sequences = (digest_file_object(inp.XYZIN) or {}).get("sequences") or {}
            return "\n".join(f">{chain}\n{seq}" for chain, seq in sequences.items())
        return ""

    def validity(self):
        # A file given stands for the text: the text is read from it at run
        # time when none was entered
        if self._files_given() and not str(self.container.controlParameters.SEQUENCETEXT).strip():
            self.container.controlParameters.SEQUENCETEXT.setQualifier('allowUndefined', True)
        return super().validity()

    def startProcess(self):
        from ccp4i2.wrappers.ProvideAlignment.script.ProvideAlignment import importAlignment

        if not str(self.container.controlParameters.SEQUENCETEXT).strip():
            text = self._text_from_files()
            if text.strip():
                self.container.controlParameters.SEQUENCETEXT.set(text)

        root = etree.Element('ProvideSequence')
        
        # Create a temporary file to store the sequence(s) that will be used
        tempFile = tempfile.NamedTemporaryFile(suffix='.txt',delete=False)
        tempFile.file.write(self.container.controlParameters.SEQUENCETEXT.__str__().encode('utf-8'))
        tempFile.close()
        
        #Attempt to interpret that as an alignment and/or stack of sequences
        alignment, format, commentary = importAlignment(tempFile.name)
        
        commentaryNode = etree.SubElement(root,"Commentary")
        commentaryNode.text = commentary.getvalue()
        
        if alignment is None:
            with open(self.makeFileName('PROGRAMXML'),'w') as programXML:
                CCP4Utils.writeXML(programXML,etree.tostring(root, pretty_print=True))
            self.reportStatus(CPluginScript.UNSATISFACTORY)
            return
        
        formatNode = etree.SubElement(root,'Format')
        formatNode.text = format
        
        from Bio import SeqIO
        for iSeq, seq in enumerate(alignment):
            outputList = self.container.outputData.SEQUENCEFILE_LIST
            outputList.append(outputList.makeItem())
            outputFile = outputList[-1]
            outputFile.setFullPath(os.path.normpath(os.path.join(self.getWorkDirectory(),'SEQUENCE'+str(iSeq)+'.fasta')))
            outputString = StringIO()
            with open(outputFile.__str__(),'w') as outputFileHandle:
                SeqIO.write([seq],outputFileHandle,'fasta')
            # Biopython's description begins with the id, so "id-description"
            # read "4HG7:A|PDBID|CHAIN|SEQUENCE-4HG7:A|PDBID|CHAIN|SEQUENCE".
            # Name it as ClustalW does (a PDB header becomes "4HG7_A"), then
            # whatever the description adds.
            from ccp4i2.wrappers.clustalw.script.clustalw import sequence_name
            rest = seq.description[len(seq.id):].strip() if seq.description.startswith(seq.id) \
                else seq.description.strip()
            outputFile.annotation = sequence_name(seq.id) + (' ' + rest if rest else '')
        
            sequenceElement = etree.SubElement(root,'Sequence')
            outputString = StringIO()
            SeqIO.write([seq],outputString,'fasta')
            sequenceElement.text = outputString.getvalue()
            for property in ['id','name','description','seq']:
                newElement = etree.SubElement(sequenceElement,property)
                newElement.text = str(getattr(seq,property,'Undefined'))

        with open(self.makeFileName('PROGRAMXML'),'w') as programXML:
            CCP4Utils.writeXML(programXML,etree.tostring(root, pretty_print=True))
        return CPluginScript.SUCCEEDED

    def processOutputFiles(self, *args, **kwargs):
        if len(self.container.outputData.SEQUENCEFILE_LIST) == 0:
            return CPluginScript.SUCCEEDED
        try:
            for iFile, sequenceFile in enumerate(self.container.outputData.SEQUENCEFILE_LIST):
                sequenceFile.loadFile()
                self.container.outputData.CASUCONTENTOUT.fileContent.seqList.append(CCP4ModelData.CAsuContentSeq())
                entry = self.container.outputData.CASUCONTENTOUT.fileContent.seqList[-1]
                entry.nCopies.set(1)
                entry.sequence.set(sequenceFile.fileContent.sequence)
                entry.name.set(sequenceFile.fileContent.identifier)
                entry.description.set(sequenceFile.fileContent.description)
                entry.autoSetPolymerType()
            self.container.outputData.CASUCONTENTOUT.saveFile()
            names = [str(f.annotation).split()[0] for f in self.container.outputData.SEQUENCEFILE_LIST
                     if str(f.annotation).strip()]
            self.container.outputData.CASUCONTENTOUT.annotation.set(
                'AU contents: ' + ', '.join(names) if names else 'AU contents')
        except Exception as err:
            print("Failed to create CASUCONTENTOUT with error", err)
        return CPluginScript.SUCCEEDED
