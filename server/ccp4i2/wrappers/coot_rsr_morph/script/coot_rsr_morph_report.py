from ccp4i2.report import Report


class coot_rsr_morph_report(Report):
    TASKNAME = 'coot_rsr_morph'
    RUNNING = False
    USEPROGRAMXML = True

    def __init__(self, xmlnode=None, jobInfo={}, jobStatus=None, **kw):
        super().__init__(xmlnode=xmlnode, jobInfo=jobInfo, **kw)
        # How far the model moved, measured from the input and output models.
        # (The report said "Full reporting is not yet available in this task".)
        moved = self.xmlnode.find(".//Shifts") if self.xmlnode is not None else None
        if moved is None or not moved.get("atoms"):
            self.addText(text="Coot real-space morphing finished.")
            return
        self.addDiv().addText(text=(
            "Coot real-space morphing moved %s atoms: RMS shift %s Å, mean %s Å, "
            "largest %s Å." % (moved.get("atoms"), moved.get("rms"),
                                   moved.get("mean"), moved.get("max"))))
        residues = moved.findall("Residue")
        if residues:
            self.addDiv().addText(text="The residues that moved most:")
            table = self.addTable()
            table.addData(title="Residue", data=[r.get("name") for r in residues])
            table.addData(title="Largest atom shift (Å)", data=[r.get("max") for r in residues])
