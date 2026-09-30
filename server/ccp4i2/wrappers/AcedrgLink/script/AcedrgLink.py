import shutil

from ccp4i2.core.CCP4PluginScript import CPluginScript


def _find_linked_dimer(work_directory, link_id):
    """The regularised linked pair AceDRG builds, as (coordinates, dictionary).

    Link mode does internally what JLigand drives libcheck and refmac to do:
    it builds the two monomers joined, and regularises the result with
    servalcat. That leaves a single merged residue whose coordinates are the
    only picture of the link a user can actually look at.

    The pair is found by suffix rather than by name. AceDRG has renamed it
    once already -- UNL_for_link -> LIG_for_link -- and the wrapper went on
    naming the old file, so both outputs silently pointed at nothing. The two
    do not share a directory: coordinates land in the work directory, the
    dictionary in <LINK_ID>_TMP.
    """
    searched = [work_directory, work_directory / f"{link_id}_TMP"]
    found = {}
    for suffix in (".pdb", ".cif"):
        for directory in searched:
            if not directory.is_dir():
                continue
            matches = sorted(directory.glob(f"*_for_link{suffix}"))
            if matches:
                found[suffix] = matches[0]
                break
    return found.get(".pdb"), found.get(".cif")


def acedrg_stop_reason(work_directory, link_id):
    """Why AceDRG stopped, in its own words, or None if it left none.

    AceDRG writes the reason to <LINK_ID>_errorInfo.txt ("Comp Lys can not be
    found in ...", "atom C in monomer GLU has a total valence of 3, which is
    not allowed!"). Without this the job reported only the log lines before
    it, which never include the one that matters.
    """
    path = work_directory / f"{link_id}_errorInfo.txt"
    try:
        text = path.read_text(errors="replace").strip()
    except OSError:
        return None
    return " ".join(text.split()) or None


class AcedrgLink(CPluginScript):
    TASKNAME = "AcedrgLink"
    TASKCOMMAND = "acedrg"
    ERROR_CODES = {
        201: {'description': 'AceDRG could not make the link'},
    }

    def postProcessCheck(self, processId=None):
        status, exit_status, exit_code = super().postProcessCheck(processId)
        if status != CPluginScript.SUCCEEDED:
            reason = acedrg_stop_reason(self.workDirectory, str(self.container.inputData.LINK_ID))
            if reason:
                print("AceDRG stopped: " + reason)
                self.appendErrorReport(201, reason)
        return status, exit_status, exit_code

    def makeCommandAndScript(self):
        inp = self.container.inputData
        par = self.container.controlParameters
        self.appendCommandLine(["-L", inp.INSTRUCTION_FILE])
        self.appendCommandLine(["-o", self.workDirectory / str(inp.LINK_ID)])
        if par.EXTRA_ACEDRG_KEYWORDS.isSet():
            for line in str(par.EXTRA_ACEDRG_KEYWORDS).splitlines():
                line = line.strip()
                if len(line) > 0 and line[0] != "#":
                    self.appendCommandLine(line)

    def processOutputFiles(self):
        print("AceDRG in link mode - processing output")
        inp = self.container.inputData
        out = self.container.outputData
        out.CIF_OUT.fullPath = self.workDirectory / f"{inp.LINK_ID}_link.cif"
        out.CIF_OUT.annotation = f"Link dictionary: {inp.ANNOTATION or inp.LINK_ID}"
        # The linked pair, so the job's Moorhen view can show the link in 3D.
        # Left unset when AceDRG wrote none: an output pointing at a file that
        # is not there is gleaned as nothing and reported as nothing.
        dimer_pdb, dimer_cif = _find_linked_dimer(self.workDirectory, str(inp.LINK_ID))
        if dimer_pdb is not None:
            out.UNL_PDB.fullPath = self._adopt(dimer_pdb)
            out.UNL_PDB.annotation = f"Linked pair: {inp.ANNOTATION or inp.LINK_ID}"
        if dimer_cif is not None:
            out.UNL_CIF.fullPath = self._adopt(dimer_cif)
            out.UNL_CIF.annotation = f"Linked pair dictionary: {inp.ANNOTATION or inp.LINK_ID}"

    def _adopt(self, path):
        """Bring an output out of AceDRG's scratch directory into the job's own.

        A CDataFile records a base name against the job directory, so it
        cannot express a file one level down in <LINK_ID>_TMP: the path comes
        back pointing at the job directory, where nothing of that name exists,
        and the gleaner publishes nothing. Copying is also the safer place for
        it -- the scratch directory is AceDRG's to delete.
        """
        if path.parent == self.workDirectory:
            return path
        destination = self.workDirectory / path.name
        shutil.copyfile(path, destination)
        return destination
