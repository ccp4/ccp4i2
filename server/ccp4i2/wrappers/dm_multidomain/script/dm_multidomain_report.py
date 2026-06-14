from ccp4i2.report import Report


class dm_multidomain_report(Report):
    """Report for multi-domain NCS averaging (dm).

    Renders the program.xml substrate written by dm_multidomain.processOutputFiles:
      - per-domain NCS averaging correlations (initial -> final, with verdict)
      - per-cycle convergence note
      - the native dm loggraph tables (completeness, mean FOM & dphi vs
        resolution, Free-R vs cycle)
    """
    TASKNAME = 'dm_multidomain'
    RUNNING = False

    def __init__(self, xmlnode=None, jobInfo={}, jobStatus=None, **kw):
        Report.__init__(self, xmlnode=xmlnode, jobInfo=jobInfo, **kw)
        if jobStatus is None or jobStatus.lower() == 'nooutput':
            return
        self.defaultReport()

    def defaultReport(self, parent=None):
        if parent is None:
            parent = self
        parent.addResults()

        # ---- per-domain NCS averaging correlations -------------------------
        try:
            domains = self.xmlnode.findall('NCSCorrelations/Domain')
            if domains:
                parent.append("<h3>NCS averaging correlation by domain</h3>")
                rows = ["<tr><th>Domain</th><th>Initial</th>"
                        "<th>Final</th><th>Status</th></tr>"]
                warned = False
                for d in domains:
                    status = d.get('status', 'OK')
                    warn = 'WARNING' in status.upper()
                    warned = warned or warn
                    style = " style='color:#b26a00;'" if warn else ""
                    rows.append(
                        f"<tr{style}><td>{d.get('number')}</td>"
                        f"<td>{d.get('initial')}</td><td>{d.get('final')}</td>"
                        f"<td>{status}</td></tr>")
                parent.append(
                    "<table border='1' cellpadding='4' "
                    "style='border-collapse:collapse;'>" + "".join(rows)
                    + "</table>")
                if warned:
                    parent.append(
                        "<p><b>Note:</b> a low initial NCS correlation usually "
                        "means the operators or mask for that domain need "
                        "checking (wrong copy assignment, mask overlap, or a "
                        "domain that does not obey the NCS).</p>")
                else:
                    parent.append("<p>All averaged domains showed improving NCS "
                                  "correlation through density modification.</p>")
        except Exception as e:
            parent.append(f"<p>Problem reporting NCS correlations: {e}</p>")

        # ---- per-cycle progress plots -------------------------------------
        try:
            cycles = self.xmlnode.findall('PerCycle/Cycle')
            if cycles:
                # columns present (Number first), in child order
                cols = [c.tag for c in cycles[0]]
                idx = {tag: i + 1 for i, tag in enumerate(cols)}
                graph = parent.addFlotGraph(
                    title="Progress by cycle", xmlnode=self.xmlnode,
                    select=".//PerCycle/Cycle",
                    style="width:500px;height:320px;margin:0 auto;border:0px;")
                for tag in cols:
                    graph.addData(title=tag, select=tag)

                # Plot 1: map quality (FOM + per-domain NCS correlation)
                quality = [t for t in cols
                           if t == 'FOM' or t.startswith('Corr_')]
                if quality:
                    p = graph.addPlotObject()
                    p.append('title', 'FOM and per-domain NCS correlation')
                    p.append('plottype', 'xy')
                    p.append('xintegral', 'true')
                    p.append('xlabel', 'Cycle')
                    p.append('yrange', min='0.0', max='1.0')
                    for tag in quality:
                        p.append('plotline', xcol=idx['Number'], ycol=idx[tag])

                # Plot 2: convergence (perturbation gamma)
                if 'Gamma' in idx:
                    p = graph.addPlotObject()
                    p.append('title', 'Perturbation gamma (convergence)')
                    p.append('plottype', 'xy')
                    p.append('xintegral', 'true')
                    p.append('xlabel', 'Cycle')
                    p.append('plotline', xcol=idx['Number'], ycol=idx['Gamma'])
                parent.addDiv(style="clear:both;")
        except Exception as e:
            parent.append(f"<p>Problem reporting per-cycle graphs: {e}</p>")

        # ---- native dm loggraph tables ------------------------------------
        try:
            fold = parent.addFold(label='Graphs from the dm log')
            group = fold.addFlotGraphGroup(
                style="width:450px;height:300px;margin:0 auto;border:0px;")
            for table in self.xmlnode.findall(
                    ".//CCP4ApplicationOutput/CCP4Table"):
                graph = group.addFlotGraph(xmlnode=table,
                                           title=table.get("title"))
                graph.addPimpleData(xmlnode=table)
            parent.addDiv(style="clear:both;")
        except Exception as e:
            parent.append(f"<p>Problem reporting log graphs: {e}</p>")
