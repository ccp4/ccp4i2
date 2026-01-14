from ccp4i2.report import Report


class coot_find_ligand_report(Report):
    TASKNAME = "coot_find_ligand"
    RUNNING = False

    def __init__(self, xmlnode=None, jobInfo={}, jobStatus=None, **kw):
        super().__init__(xmlnode=xmlnode, jobInfo=jobInfo, **kw)

        added = xmlnode.findall(".//AddedLigand")
        self.addText(text=f"Number of ligands added: {len(added)}")

        table = self.addTable()
        indices = []
        volumes = []
        scores = []
        for i, ligand in enumerate(added):
            indices.append(i + 1)
            volume = float(ligand.get("clusterVolume", "nan"))
            score = float(ligand.get("fittingScore", "nan"))
            volumes.append(f"{volume:.1f}")
            scores.append(f"{score:.3f}")
        table.addData(title="#", data=indices)
        table.addData(title="Cluster Volume", data=volumes)
        table.addData(title="Fitting Score", data=scores)
