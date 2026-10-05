from ccp4i2.report import Report

LABELS = (("added", "Added as atoms"), ("converted", "Methionines made MSE (Se on SD)"),
          ("on_model", "Already model atoms (not doubled)"),
          ("not_placed", "Not placed (S or Se with no residue to belong to, or frame not found)"),
          ("clashes", "Left out: too close to another model atom"),
          ("removed", "Residues removed: built into a heavy atom's density"))


class add_substructure_report(Report):
    TASKNAME = "add_substructure"
    RUNNING = False

    def __init__(self, xmlnode=None, jobInfo={}, jobStatus=None, **kw):
        Report.__init__(self, xmlnode=xmlnode, jobInfo=jobInfo, jobStatus=jobStatus, **kw)
        node = xmlnode.find(".//CompleteModel") if xmlnode is not None else None
        if node is None:
            return
        fold = self.addFold(label="Sites", brief="Sites", initiallyOpen=True)
        origin = node.find("Origin")
        if origin is not None:
            text = origin.get("note", "")
            if origin.get("transform"):
                text += f" ({origin.get('transform')})"
            fold.addText(text=f"Origin and hand: {text}")
            fold.addDiv(style="clear:both;")
        if node.get("all_mse", "0") != "0":
            fold.addText(text=f"Every methionine made selenomethionine: {node.get('all_mse')}")
            fold.addDiv(style="clear:both;")
        for kind, label in LABELS:
            items = [e.text for e in node.findall(kind)]
            if items:
                table = fold.addTable(transpose=False)
                table.addData(title=f"{label} ({len(items)})", data=items)
                fold.addDiv(style="clear:both;")
