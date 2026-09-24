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

    # Mid-tone and distinguishable on a light report page; deliberately the
    # same order as the task interface's palette so a body keeps its colour
    # from setting the job up to reading it back.
    BODY_COLOURS = ['#1976d2', '#9c27b0', '#2e7d32', '#ed6c02', '#0288d1',
                    '#d32f2f']

    def defaultReport(self, parent=None):
        if parent is None:
            parent = self
        parent.addResults()

        self._add_body_map(parent)

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

    # ------------------------------------------------------------------ #
    def _add_body_map(self, parent):
        """Draw the rigid bodies on the reference copy's residues.

        A report that only gives correlations per domain number leaves the
        reader to remember which domain was which. Here the partition is the
        picture: a track per entity, a coloured block per body, a hatched
        block where two bodies claimed the same residues, and the bare track
        showing residues no body averaged.
        """
        try:
            tracks = self.xmlnode.findall('BodyMap/Track')
            bodies = self.xmlnode.findall('BodyMap/Body')
            if not tracks or not bodies:
                return

            width, label_w, pad = 640, 110, 30
            track_h, gap = 18, 16
            plot_w = width - label_w - pad
            height = len(tracks) * (track_h + gap)

            parts = [
                "<h3>What was averaged</h3>",
                f"<svg viewBox='0 0 {width} {height + 6}' "
                f"style='width:100%;max-width:{width}px;height:auto;' "
                f"role='img' aria-label='Rigid bodies on the reference copy'>",
                "<defs><pattern id='dmclash' width='6' height='6' "
                "patternUnits='userSpaceOnUse' patternTransform='rotate(45)'>"
                "<rect width='6' height='6' fill='#d32f2f' fill-opacity='0.18'/>"
                "<line x1='0' y1='0' x2='0' y2='6' stroke='#d32f2f' "
                "stroke-width='2'/></pattern></defs>",
            ]

            for row, track in enumerate(tracks):
                role = track.get('role')
                lo, hi = int(track.get('first')), int(track.get('last'))
                span = max(1, hi - lo)
                y = row * (track_h + gap)

                def x(residue):
                    clamped = min(max(residue, lo), hi)
                    return label_w + (clamped - lo) / span * plot_w

                name = (f"chain {track.get('chain')}" if role == '_'
                        else f"{role} ({track.get('chain')})")
                parts.append(
                    f"<text x='0' y='{y + track_h - 4}' font-size='12' "
                    f"fill='currentColor'>{name}</text>")
                parts.append(
                    f"<rect x='{label_w}' y='{y}' width='{plot_w}' "
                    f"height='{track_h}' rx='3' fill='currentColor' "
                    f"fill-opacity='0.06' stroke='currentColor' "
                    f"stroke-opacity='0.25'/>")

                spans = []
                for i, body in enumerate(bodies):
                    for segment in body.findall('Segment'):
                        if segment.get('role') != role:
                            continue
                        spans.append((int(segment.get('lo')),
                                      int(segment.get('hi')), i))
                for s_lo, s_hi, i in spans:
                    colour = self.BODY_COLOURS[i % len(self.BODY_COLOURS)]
                    opacity = ('0.35' if bodies[i].get('mode') == 'exclude'
                               else '0.85')
                    parts.append(
                        f"<rect x='{x(s_lo):.1f}' y='{y}' "
                        f"width='{max(2, x(s_hi) - x(s_lo)):.1f}' "
                        f"height='{track_h}' rx='3' fill='{colour}' "
                        f"fill-opacity='{opacity}'>"
                        f"<title>body {i + 1}: {s_lo}-{s_hi}</title></rect>")
                for a in range(len(spans)):
                    for b in range(a + 1, len(spans)):
                        if spans[a][2] == spans[b][2]:
                            continue
                        start = max(spans[a][0], spans[b][0])
                        end = min(spans[a][1], spans[b][1])
                        if start > end:
                            continue
                        parts.append(
                            f"<rect x='{x(start):.1f}' y='{y}' "
                            f"width='{max(2, x(end) - x(start)):.1f}' "
                            f"height='{track_h}' rx='3' fill='url(#dmclash)' "
                            f"stroke='#d32f2f'/>")

                parts.append(
                    f"<text x='{label_w}' y='{y + track_h + 11}' "
                    f"font-size='10' fill='currentColor' "
                    f"fill-opacity='0.7'>{lo}</text>")
                parts.append(
                    f"<text x='{label_w + plot_w}' y='{y + track_h + 11}' "
                    f"font-size='10' text-anchor='end' fill='currentColor' "
                    f"fill-opacity='0.7'>{hi}</text>")

            parts.append("</svg>")

            rows = ["<tr><th>Body</th><th>Residues</th><th>Mode</th>"
                    "<th>Superposition RMSD on each copy</th></tr>"]
            for i, body in enumerate(bodies):
                colour = self.BODY_COLOURS[i % len(self.BODY_COLOURS)]
                segments = ", ".join(
                    (f"{s.get('lo')}-{s.get('hi')}" if s.get('role') == '_'
                     else f"{s.get('role')}:{s.get('lo')}-{s.get('hi')}")
                    for s in body.findall('Segment'))
                fits = ", ".join(f"{f.get('copy')} {f.get('rmsd')} &#8491;"
                                 for f in body.findall('Fit')) or "&#8212;"
                swatch = (f"<span style='display:inline-block;width:10px;"
                          f"height:10px;border-radius:2px;background:{colour};"
                          f"margin-right:6px;'></span>")
                rows.append(
                    f"<tr><td>{swatch}{i + 1}</td><td>{segments}</td>"
                    f"<td>{body.get('mode')}</td><td>{fits}</td></tr>")
            parts.append(
                "<table border='1' cellpadding='4' "
                "style='border-collapse:collapse;'>" + "".join(rows)
                + "</table>")
            parent.append("".join(parts))
        except Exception as e:
            parent.append(f"<p>Problem drawing the rigid-body map: {e}</p>")
