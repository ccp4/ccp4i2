import sys

from ....report.CCP4ReportParser import Report


class validate_protein_report(Report):
    TASKNAME = 'validate_protein'

    def __init__(self, xmlnode=None, jobInfo={ }, **kw):
        Report.__init__(self, xmlnode=xmlnode, jobInfo=jobInfo, **kw)

        self.append('<title>Multimetric Validation Report</title>')

        results = self.addResults()

        if len(self.xmlnode.findall('.//Iris'))>0:
            self.add_iris_panel()

        if len(self.xmlnode.findall('.//Molprobity'))>0:
            self.add_molprobity()
        
        if len(self.xmlnode.findall('.//B_factors'))>0:
            self.add_b_factors()
        
        if len(self.xmlnode.findall('.//Ramachandran'))>0:
            self.add_ramachandran()


    def add_iris_panel(self, parent=None):
        if parent is None:
            parent = self
        section_div = parent.addDiv()
        iris_panel_svg = self.xmlnode.findall('.//Iris/Panel_svg')[0].text
        section_div.append(iris_panel_svg)
        section_div.addDiv(style='clear:both;')
        parent.append('<br/><br/>')


    def add_molprobity(self, parent=None):
        if parent is None:
            parent = self
        section_div = parent.addDiv()
        section_div.append('<p style="font-weight:bold; font-size: 14px;"><b>MolProbity Analyses</b></p>')

        fold = section_div.addFold(label='Summary', initiallyOpen=True)
        table_summary = fold.addTable(title='MolProbity Summary', select='.//Molprobity/Summary', transpose=True)
        table_summary.addData(title='Ramachandran outliers', select='Ramachandran_outliers')
        table_summary.addData(title='Ramachandran favoured', select='Ramachandran_favoured')
        table_summary.addData(title='Rotamer outliers', select='Rotamer_outliers')
        table_summary.addData(title='CBeta deviations', select='CBeta_deviations')
        table_summary.addData(title='Clashscore', select='Clashscore')
        table_summary.addData(title='RMS bonds', select='RMS_bonds')
        table_summary.addData(title='RMS angles', select='RMS_angles')
        table_summary.addData(title='Molprobity score', select='Molprobity_score')
        section_div.append('<br/>')

        # Show tables conditionally
        n_rmchn = len(self.xmlnode.findall('.//Molprobity/Ramachandran_outliers/Outlier'))
        n_omega = len(self.xmlnode.findall('.//Molprobity/Nonplanar_omegas/Outlier'))
        n_rotam = len(self.xmlnode.findall('.//Molprobity/Rotamer_outliers/Outlier'))
        n_cbeta = len(self.xmlnode.findall('.//Molprobity/CBeta_outliers/Outlier'))
        n_clash = len(self.xmlnode.findall('.//Molprobity/Clashes/Outlier'))
        n_flips = len(self.xmlnode.findall('.//Molprobity/Side_chain_flips/Outlier'))

        fold = section_div.addFold(label='Detailed breakdown of outliers', initiallyOpen=False)

        rmchn_div = fold.addDiv()
        if n_rmchn > 0:
            rmchn_div.append('<p><b>Ramachandran outliers: </b></p>')
            table_outliers = rmchn_div.addTable(title='List of outliers', transpose=False)
            table_outliers.addData(title='Chain', data=[ x.attrib["chain"] for x in self.xmlnode.findall('.//Molprobity/Ramachandran_outliers/Outlier') if "chain" in x.attrib])
            table_outliers.addData(title='Seqnum', data=[ x.attrib["seqnum"] for x in self.xmlnode.findall('.//Molprobity/Ramachandran_outliers/Outlier') if "seqnum" in x.attrib])
            table_outliers.addData(title='Name', data=[ x.attrib["name"] for x in self.xmlnode.findall('.//Molprobity/Ramachandran_outliers/Outlier') if "name" in x.attrib])
            table_outliers.addData(title='Score', data=[ x.attrib["score"] for x in self.xmlnode.findall('.//Molprobity/Ramachandran_outliers/Outlier') if "score" in x.attrib])

        omega_div = fold.addDiv()
        if n_omega > 0:
            omega_div.append('<p><b>Non-planar omega (peptide bond) angles: </b></p>')
            table_outliers = omega_div.addTable(title='List of outliers', transpose=False)
            table_outliers.addData(title='Chain', data=[ x.attrib["chain"] for x in self.xmlnode.findall('.//Molprobity/Nonplanar_omegas/Outlier') if "chain" in x.attrib])
            table_outliers.addData(title='Seqnum', data=[ x.attrib["seqnum"] for x in self.xmlnode.findall('.//Molprobity/Nonplanar_omegas/Outlier') if "seqnum" in x.attrib])
            table_outliers.addData(title='Name', data=[ x.attrib["name"] for x in self.xmlnode.findall('.//Molprobity/Nonplanar_omegas/Outlier') if "name" in x.attrib])

        rotam_div = fold.addDiv()
        if n_rotam > 0:
            rotam_div.append('<p><b>Rotamer outliers: </b></p>')
            table_outliers = rotam_div.addTable(title='List of outliers', transpose=False)
            table_outliers.addData(title='Chain', data=[ x.attrib["chain"] for x in self.xmlnode.findall('.//Molprobity/Rotamer_outliers/Outlier') if "chain" in x.attrib])
            table_outliers.addData(title='Seqnum', data=[ x.attrib["seqnum"] for x in self.xmlnode.findall('.//Molprobity/Rotamer_outliers/Outlier') if "seqnum" in x.attrib])
            table_outliers.addData(title='Name', data=[ x.attrib["name"] for x in self.xmlnode.findall('.//Molprobity/Rotamer_outliers/Outlier') if "name" in x.attrib])
            table_outliers.addData(title='Score', data=[ x.attrib["score"] for x in self.xmlnode.findall('.//Molprobity/Rotamer_outliers/Outlier') if "score" in x.attrib])

        cbeta_div = fold.addDiv()
        if n_cbeta > 0:
            cbeta_div.append('<p><b>C<sup><i>&#946;</i></sup> outliers: </b></p>')
            table_outliers = cbeta_div.addTable(title='List of outliers', transpose=False)
            table_outliers.addData(title='Chain', data=[ x.attrib["chain"] for x in self.xmlnode.findall('.//Molprobity/CBeta_outliers/Outlier') if "chain" in x.attrib])
            table_outliers.addData(title='Seqnum', data=[ x.attrib["seqnum"] for x in self.xmlnode.findall('.//Molprobity/CBeta_outliers/Outlier') if "seqnum" in x.attrib])
            table_outliers.addData(title='Name', data=[ x.attrib["name"] for x in self.xmlnode.findall('.//Molprobity/CBeta_outliers/Outlier') if "name" in x.attrib])
            table_outliers.addData(title='Deviation', data=[ x.attrib["score"] for x in self.xmlnode.findall('.//Molprobity/CBeta_outliers/Outlier') if "score" in x.attrib])

        flips_div = fold.addDiv()
        if n_flips > 0:
            flips_div.append('<p><b>Suggested side-chain flips: </b></p>')
            table_outliers = flips_div.addTable(title='List of flips', transpose=False)
            table_outliers.addData(title='Chain', data=[ x.attrib["chain"] for x in self.xmlnode.findall('.//Molprobity/Side_chain_flips/Outlier') if "chain" in x.attrib])
            table_outliers.addData(title='Seqnum', data=[ x.attrib["seqnum"] for x in self.xmlnode.findall('.//Molprobity/Side_chain_flips/Outlier') if "seqnum" in x.attrib])
            table_outliers.addData(title='Name', data=[ x.attrib["name"] for x in self.xmlnode.findall('.//Molprobity/Side_chain_flips/Outlier') if "name" in x.attrib])

        clash_div = fold.addDiv()
        if n_clash > 0:
            clash_div.append('<p><b>Atomic clashes: </b></p>')
            table_outliers = clash_div.addTable(title='List of clashes', transpose=False)
            table_outliers.addData(title='First atom', data=[ x.attrib["first_atom"] for x in self.xmlnode.findall('.//Molprobity/Clashes/Outlier') if "first_atom" in x.attrib])
            table_outliers.addData(title='Second atom', data=[ x.attrib["second_atom"] for x in self.xmlnode.findall('.//Molprobity/Clashes/Outlier') if "second_atom" in x.attrib])
            table_outliers.addData(title='Overlap', data=[ x.attrib["overlap"] for x in self.xmlnode.findall('.//Molprobity/Clashes/Outlier') if "overlap" in x.attrib])

        section_div.addDiv(style='clear:both;')
        parent.append('<br/><br/>')


    def add_b_factors(self, parent=None):
        if parent is None:
            parent = self
        b_factor_averages = { }
        for list_name in ('all', 'amino_acids', 'main_chains', 'side_chains', 'non_amino_acids', 'waters', 'ligands', 'ions'):
            b_factor_averages[list_name] = { }
            chain_ids = [ x.attrib["chain"] for x in self.xmlnode.findall('.//B_factors/' + list_name) if "chain" in x.attrib]
            means = [ x.attrib["mean"] for x in self.xmlnode.findall('.//B_factors/' + list_name) if "mean" in x.attrib]
            stds = [ x.attrib["std"] for x in self.xmlnode.findall('.//B_factors/' + list_name) if "std" in x.attrib]
            ns = [ x.attrib["n"] for x in self.xmlnode.findall('.//B_factors/' + list_name) if "n" in x.attrib]
            for chain_id, mean, std, n in zip(chain_ids, means, stds, ns):
                chain_id = int(chain_id) if chain_id != 'All' else 'All'
                mean = round(float(mean), 2) if mean not in ('None', 'nan') else None
                std = round(float(std), 2) if std not in ('None', 'nan') else None
                n = int(n) if n not in ('None', 'nan') else None
                b_factor_averages[list_name][chain_id] = (mean, std, n)

        section_div = parent.addDiv()
        section_div.append('<style> th { white-space: nowrap; } </style>')
        section_div.append('<p style="font-weight:bold; font-size: 14px;"><b>B-factor Analyses</b></p>')
        #table_div = section_div.addDiv(style='float:left; margin-left: 30px; width:25%;')
        fold = section_div.addFold(label='Whole model', initiallyOpen=True)
        table = fold.addTable(title='B-factor analysis: whole model', transpose=False)
        table.addData(title='', data=('Mean', 'Stdev.', 'N'))
        table.addData(title='All monomers', data=b_factor_averages['all']['All'])
        table.addData(title='All amino acids', data=b_factor_averages['amino_acids']['All'])
        table.addData(title='Main chains', data=b_factor_averages['main_chains']['All'])
        table.addData(title='Side chains', data=b_factor_averages['side_chains']['All'])
        table.addData(title='All non-amino acids', data=b_factor_averages['non_amino_acids']['All'])
        table.addData(title='Waters', data=b_factor_averages['waters']['All'])
        table.addData(title='Ligands', data=b_factor_averages['ligands']['All'])
        table.addData(title='Ions', data=b_factor_averages['ions']['All'])
        section_div.append('<br/>')

        chain_count = int(self.xmlnode.findall('.//Model_info/Chain_count')[0].text)
        for chain_id in range(chain_count):
            if chain_id in b_factor_averages['all']:
                fold = section_div.addFold(label='Chain ' + str(chain_id+1), initiallyOpen=False)
                table = fold.addTable(title='B-factor analysis: chain ' + str(chain_id+1), transpose=False)
                table.addData(title='', data=('Mean', 'Stdev.', 'N'))
                table.addData(title='All monomers', data=b_factor_averages['all'][chain_id])
                table.addData(title='All amino acids', data=b_factor_averages['amino_acids'][chain_id])
                table.addData(title='Main chains', data=b_factor_averages['main_chains'][chain_id])
                table.addData(title='Side chains', data=b_factor_averages['side_chains'][chain_id])
                table.addData(title='All non-amino acids', data=b_factor_averages['non_amino_acids'][chain_id])
                table.addData(title='Waters', data=b_factor_averages['waters'][chain_id])
                table.addData(title='Ligands', data=b_factor_averages['ligands'][chain_id])
                table.addData(title='Ions', data=b_factor_averages['ions'][chain_id])
                section_div.append('<br/>')

        section_div.addDiv(style='clear:both;')
        parent.append('<br/><br/>')


    def add_ramachandran(self, parent=None):
        if parent is None:
            parent = self
        section_div = parent.addDiv()
        section_div.append('<p style="font-weight:bold; font-size:14px;"><b>Ramachandran Analyses</b></p>')

        background_PRO = 'http://127.0.0.1:43434/report_files/0.1.0/rama_pro.png'
        background_GLY = 'http://127.0.0.1:43434/report_files/0.1.0/rama_gly.png'
        background_RST = 'http://127.0.0.1:43434/report_files/0.1.0/rama_rst.png'

        graph_div = parent.addDiv(style='float:left; margin-left:50px;')
        rama_graph = graph_div.addFlotGraph(title='Ramachandran Plot', select='.//Ramachandran', style='height:500px; width:500px; border:0px; float:left; padding:10px; padding-left:15px;')

        for resName, resCheck in [
            ("PRO", lambda x: x.attrib.get("type") == "PRO"),
            ("GLY", lambda x: x.attrib.get("type") == "GLY"),
            ("RST", lambda x: x.attrib.get("type") not in {"PRO", "GLY"}),
        ]:
            for section in ["Favoured", "Allowed", "Outliers"]:
                for angle in ["Phi", "Psi"]:
                    title = f"{resName}_{section}_{angle}"
                    xmlPath = f".//Ramachandran/{section}/Residue"
                    data = [
                        float(x.findall(angle)[0].text.strip())
                        for x in self.xmlnode.findall(xmlPath)
                        if resCheck(x)
                    ]
                    rama_graph.addData(title=title, data=data)

        for title, colStart, background in [
            ("Non-Pro/Gly", 13, background_RST),
            ("Proline", 1, background_PRO),
            ("Glycine", 7, background_GLY),
        ]:
            p = rama_graph.addPlotObject()
            description = (
                "This graph shows the dihedral Phi and Psi angles for all residues, "
                "coloured according to Ramachandran's criterion. "
                "Source: Richardsons' Top 500 structures."
            )
            p.append('description', description)
            p.append('background').text = background
            p.append('title', f'Ramachandran plot [{title}]')
            p.append('plottype', 'xy')
            p.append('showlegend', 'false')
            p.append('xintegral', 'true')
            p.append('yintegral', 'true')
            p.append('xlabel', 'Phi')
            p.append('ylabel', 'Psi')
            p.append('xrange', min=-180.0, max=180.0)
            p.append('yrange', min=-180.0, max=180.0)
            l = p.append('plotline', xcol=colStart+0, ycol=colStart+1)
            l.append('colour', 'green')
            l.append('linestyle', '.')
            l.append('symbolsize', '1')
            l.append('symbol', '.')
            l = p.append('plotline', xcol=colStart+2, ycol=colStart+3)
            l.append('colour', 'orange')
            l.append('linestyle', '.')
            l.append('symbolsize', '3')
            l.append('symbol', 'o')
            l = p.append('plotline', xcol=colStart+4, ycol=colStart+5)
            l.append('colour', 'red')
            l.append('linestyle', '.')
            l.append('symbolsize', '10')
            l.append('symbol', 'x')

        n_residues = int(self.xmlnode.findall('.//Ramachandran/Totals/Residues')[0].text)
        n_favoured = int(self.xmlnode.findall('.//Ramachandran/Totals/Favoured')[0].text)
        n_allowed  = int(self.xmlnode.findall('.//Ramachandran/Totals/Allowed')[0].text)
        n_outliers = int(self.xmlnode.findall('.//Ramachandran/Totals/Outliers')[0].text)
        n_na = int(self.xmlnode.findall('.//Ramachandran/Totals/NA')[0].text)

        tab_div = parent.addDiv(style='float:left; margin-left:50px;')

        if int(n_outliers) > 0 :
            tab_div.append('<p ><b>Outliers: </b></p>')
            table_outliers = tab_div.addTable(title='List of outliers', transpose=False, downloadable=True)
            table_outliers.addData(title='Chain', data=[ x.attrib["chain"] for x in self.xmlnode.findall('.//Ramachandran/Outliers/Residue') if "chain" in x.attrib])
            table_outliers.addData(title='Name',  data=[ x.attrib["type"] for x in self.xmlnode.findall('.//Ramachandran/Outliers/Residue') if "type" in x.attrib])
            table_outliers.addData(title='Residue',    data=[ x.attrib["seqnum"] for x in self.xmlnode.findall('.//Ramachandran/Outliers/Residue') if "seqnum" in x.attrib])

        text ='<b>%s</b> residues have been analysed.<br/><br/>' % n_residues
        text +='In <b><font color=\'green\'>favoured</font></b> regions: <b>%s (%0.2f%%)</b><br/>' % (n_favoured, 100 * n_favoured / n_residues)
        text +='In <b><font color=\'orange\'>allowed</font></b> regions: <b>%s (%0.2f%%)</b><br/>' % (n_allowed, 100 * n_allowed / n_residues)
        text +='In <b><font color=\'red\'>high-energy</font></b> backbone conformations (outliers): <b>%s (%0.2f%%)</b>' % (n_outliers, 100 * n_outliers / n_residues)
        tab_div.append(text)

        section_div.addDiv(style='clear:both;')
        parent.append('<br/><br/><br/><br/>')


if __name__ == '__main__':
    validate_protein_report(xmlFile=sys.argv[1], jobId=sys.argv[2])
