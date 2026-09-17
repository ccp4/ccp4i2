from ccp4i2.report import Report


class moorhen_report(Report):
    TASKNAME = "moorhen"
    RUNNING = False

    def __init__(self, xmlnode=None, jobInfo={}, jobStatus=None, **kw):
        Report.__init__(self, xmlnode=xmlnode, jobInfo=jobInfo, **kw)
        self.addText(text="Moorhen session finished")
        models = xmlnode.findall(".//moorhen/number_output_files")
        if models:
            self.append("<br/>")
            self.addText(text="Models saved: " + (models[0].text or "0"))
        dicts = xmlnode.findall(".//moorhen/number_output_dicts")
        if dicts:
            self.append("<br/>")
            self.addText(text="Dictionaries saved: " + (dicts[0].text or "0"))
