"""
    dnatco_report.py: CCP4 GUI Project
    Copyright (C) 2025 MRC-LMB
    Author: Martin Maly
    
    This library is free software: you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public License
    version 3, modified in accordance with the provisions of the
    license to address the requirements of UK law.
    
    You should have received a copy of the modified GNU Lesser General
    Public License along with this library.  If not, copies may be
    downloaded from http://www.ccp4.ac.uk/ccp4license.php
    
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.
"""

from report.CCP4ReportParser import *
from core import CCP4Modules
from pathlib import Path
import gemmi


class dnatco_report(Report):
    TASKNAME= 'dnatco'
    RUNNING = True

    def __init__(self, xmlnode=None, jobInfo={}, jobStatus=None, **kw):
        Report.__init__( self, xmlnode=xmlnode, jobInfo=jobInfo, jobStatus=jobStatus, **kw)

        if jobStatus is None or jobStatus.lower() == 'nooutput': return

        projectid = self.jobInfo.get("projectid", None)
        jobNumber = self.jobInfo.get("jobnumber", None)
        jobId = self.jobInfo.get("jobid", None)
        jobDirectory = CCP4Modules.PROJECTSMANAGER().jobDirectory(jobId=jobId)
        self.jobLog = str(Path(jobDirectory) / "log.txt")
        if jobStatus is not None and jobStatus.lower() == "running":
            self.runningReport(parent=self)
        else:
            self.defaultReport(parent=self)


    def runningReport(self, parent=None):
        if parent is None:
            parent = self
        if Path(self.jobLog).is_file():
            jobLogFold = parent.addFold(label="DNATCO log", initiallyOpen=True)
            jobLogFold.addPre("DNATCO is running...")


    def defaultReport(self, parent=None, ciffilePaths=[]):

        if parent is None:
            parent = self
        self.addDiv(style="clear:both;")  # gives space for the title

        if ciffilePaths:
            if len(ciffilePaths) == 2:
                compare_two = True
            else:
                compare_two = False
        else:
            # only one mmCIF file from this job
            compare_two = False
            if self.jobInfo:
                if 'filenames' in self.jobInfo and 'CIFOUT' in self.jobInfo['filenames']:
                    ciffilePaths = [self.jobInfo['filenames']['CIFOUT']]
            
        cif_data1 = self.read_data_from_cif(ciffilePaths[0], parent)
        if compare_two:
            cif_data2 = self.read_data_from_cif(ciffilePaths[1], parent)

        # draw report
        overallFold = parent.addFold(label="Overall structure quality", initiallyOpen=True)
        noteDiv = overallFold.addDiv(style='font-size:110%;')
        noteDiv.append(
            "Assignment of conformer class (NtC) to dinucleotide steps and the quality of fit between each"
            " dinucleotide and the NtC reference structure, measured by RMSD and confal score."
            )
        table = overallFold.addTable(title="Overall structure quality", transpose=True)
        if compare_two:
            table.addData(title="", data=["Model 1", "Model 2"])
            table.addData(title="No. NtC assigned", data=[cif_data1['overall_num_classified'], cif_data2['overall_num_classified']])
            table.addData(title="No. NtC close", data=[cif_data1['overall_num_unclassified_rmsd_close'], cif_data2['overall_num_unclassified_rmsd_close']])
            table.addData(title="No. NtC unassigned", data=[cif_data1['overall_num_unclassified'], cif_data2['overall_num_unclassified']])
            table.addData(title="No. steps with RMSD < 0.5 &#197;", data=[cif_data1["overall_num_rmsd_below_0p5"], cif_data2["overall_num_rmsd_below_0p5"]])
            table.addData(title="No. steps with RMSD 0.5-1 &#197;", data=[cif_data1["overall_num_rmsd_between_0p5_and_1"], cif_data2["overall_num_rmsd_between_0p5_and_1"]])
            table.addData(title="No. steps with RMSD > 1 &#197;", data=[cif_data1["overall_num_rmsd_higher_1"], cif_data2["overall_num_rmsd_higher_1"]])
            table.addData(title="Confal score", data=[cif_data1['overall_confal_score'], cif_data2['overall_confal_score']])
            table.addData(title="Confal score percentile", data=[cif_data1['overall_confal_percentile'], cif_data2['overall_confal_percentile']])
        else:
            table.addData(title="No. NtC assigned", data=[cif_data1['overall_num_classified']])
            table.addData(title="No. NtC close", data=[cif_data1['overall_num_unclassified_rmsd_close']])
            table.addData(title="No. NtC unassigned", data=[cif_data1['overall_num_unclassified']])
            table.addData(title="No. steps with RMSD < 0.5 &#197;", data=[cif_data1["overall_num_rmsd_below_0p5"]])
            table.addData(title="No. steps with RMSD 0.5-1 &#197;", data=[cif_data1["overall_num_rmsd_between_0p5_and_1"]])
            table.addData(title="No. steps with RMSD > 1 &#197;", data=[cif_data1["overall_num_rmsd_higher_1"]])
            table.addData(title="Confal score", data=[cif_data1['overall_confal_score']])
            table.addData(title="Confal score percentile", data=[cif_data1['overall_confal_percentile']])

        if compare_two:
            model_label = " (model 1)"
        else:
            model_label = ""
        outliersFold = parent.addFold(label="Dinucleotides outliers", initiallyOpen=True)
        noteDiv = outliersFold.addDiv(style='font-size:110%;')
        noteDiv.append(
            "List all unassigned dinucleotide steps. Dinucleotide conformer (NtC),"
            " resp. conformational alphabet of nucleic acids (CANA) classes in"
            " the table below represent the closest NtC, resp. CANA class that would be assigned to the given"
            " dinucleotide if all assignment criteria were met."
            )
        if (
            (compare_two and len(cif_data2['idx_outliers']) > 0)
            or (not compare_two and len(cif_data1['idx_outliers']) > 0)
            ):
            table = outliersFold.addTable(title="Dinucleotides outliers")
            table.addData(title="Step ID", data=cif_data1['step_id_outliers'])
            table.addData(title="Chain", data=cif_data1['chain_display_outliers'])
            table.addData(title="Step", data=cif_data1['steps_outliers'])
            table.addData(title=f"Closest CANA{model_label}", data=cif_data1['assigned_CANA_outliers'])
            if compare_two:
                table.addData(title="Closest CANA (model 2)", data=cif_data2['assigned_CANA_outliers'])
            table.addData(title=f"Closest NtC1{model_label}", data=cif_data1['assigned_NtC_outliers'])
            if compare_two:
                table.addData(title="Closest NtC (model 2)", data=cif_data2['assigned_NtC_outliers'])
            table.addData(title=f"RMSD to closest NtC representative (&#197;){model_label}", data=cif_data1['rmsd_NtC_outliers'])
            if compare_two:
                table.addData(title="RMSD to closest NtC representative (&#197;) (model 1)", data=cif_data2['rmsd_NtC_outliers'])
            noteDiv = outliersFold.addDiv(style='font-size:110%;')
            noteDiv.append(
            "A more detailed analysis could be performed at the"
            " <a href='https://dnatco.datmos.org'>DNATCO web server</a>."
            )
        else:
            noteDiv = outliersFold.addDiv(style='font-size:110%;')
            noteDiv.append("No dinucleotide outliers found.")

        improvablesFold = parent.addFold(label="Improvable dinucleotide outliers", initiallyOpen=True)
        noteDiv = improvablesFold.addDiv(style='font-size:110%;')
        noteDiv.append(
            "List of unassigned dinucleotide steps that are considered"
            " sufficiently close to a representative from the Golden Set. Closeness criterion is defined as RMSD"
            " value less or equal to 0.5 Å."
            )
        if (
            (compare_two and len(cif_data2['idx_improvables']) > 0)
            or (not compare_two and len(cif_data1['idx_improvables']) > 0)
            ):
            table = improvablesFold.addTable(title="Improvable dinucleotide outliers")
            table.addData(title="Step ID", data=cif_data1['step_id_improvables'])
            table.addData(title="Chain", data=cif_data1['chain_display_improvables'])
            table.addData(title="Step", data=cif_data1['steps_improvables'])
            table.addData(title=f"Closest CANA{model_label}", data=cif_data1['assigned_CANA_improvables'])
            if compare_two:
                table.addData(title="Closest CANA (model 2)", data=cif_data2['assigned_CANA_improvables'])
            table.addData(title=f"Closest NtC{model_label}", data=cif_data1['assigned_NtC_improvables'])
            if compare_two:
                table.addData(title="Closest NtC (model 2)", data=cif_data2['assigned_NtC_improvables'])
            table.addData(title=f"RMSD to closest NtC representative (&#197;){model_label}", data=cif_data1['rmsd_NtC_improvables'])
            if compare_two:
                table.addData(title="RMSD to closest NtC representative (&#197;) (model 1)", data=cif_data2['rmsd_NtC_improvables'])
            noteDiv = improvablesFold.addDiv(style='font-size:110%;')
            noteDiv.append(
            "A more detailed analysis could be performed at the"
            " <a href='https://dnatco.datmos.org'>DNATCO web server</a>."
            )
        else:
            noteDiv = improvablesFold.addDiv(style='font-size:110%;')
            noteDiv.append("No improvable dinucleotide outliers found.")

        allDinucleotidesFold = parent.addFold(label="All dinucleotides", initiallyOpen=True)

        graphPerStep = allDinucleotidesFold.addFlotGraph(
            title="RMSD per step",
            style="height:330px; width:585px; border:0px; padding:10px; padding-left:15px; margin-right:15px;")
        graphPerStep.addData(title="step", data=cif_data1['step_id'])
        graphPerStep.addData(title="RMSD_model1(A)", data=cif_data1['rmsd'])
        if compare_two:
            graphPerStep.addData(title="RMSD_model2(A)", data=cif_data2['rmsd'])
        plotPerStep = graphPerStep.addPlotObject()
        plotPerStep.append('title', 'RMSD per step')
        plotPerStep.append('plottype', 'xy')
        plotPerStep.append('xlabel', 'step id')
        plotPerStep.append('ylabel', 'RMSD (A)')
        plotPerStep.append('legendposition', x=1, y=0)  # right bottom corner
        plotLine = plotPerStep.append('plotline', xcol=1, ycol=2)
        plotLine.append('symbolsize', '1')
        plotLine.append('linestyle', '.')
        plotLine.append('colour', 'gray')
        if compare_two:
            plotLine2 = plotPerStep.append('plotline', xcol=1, ycol=3)
            plotLine2.append('symbolsize', '1')
            plotLine2.append('linestyle', '.')
            plotLine2.append('colour', 'blue')

        table = allDinucleotidesFold.addTable(title="All dinucleotides")
        table.addData(title="Step ID", data=cif_data1['step_id'])
        table.addData(title="Chain", data=cif_data1['chain_display'])
        table.addData(title="Step", data=cif_data1['steps'])
        table.addData(title=f"Assigned CANA{model_label}", data=cif_data1['assigned_CANA'])
        if compare_two:
            table.addData(title="Assigned CANA (model 2)", data=cif_data2['assigned_CANA'])
        table.addData(title=f"Assigned NtC{model_label}", data=cif_data1['assigned_NtC'])
        if compare_two:
            table.addData(title="Assigned NtC (model 2)", data=cif_data2['assigned_NtC'])
        table.addData(title=f"Confal score{model_label}", data=cif_data1['confal_score'])
        if compare_two:
            table.addData(title="Confal score (model 2)", data=cif_data2['rmsd_NtC_assigned'])
        table.addData(title=f"RMSD to closest NtC representative (&#197;){model_label}", data=cif_data1['rmsd_NtC_assigned'])
        if compare_two:
            table.addData(title="RMSD to closest NtC representative (&#197;) (model 2)", data=cif_data2['rmsd_NtC_assigned'])

        self.addDiv(style="clear:both;")


    def read_data_from_cif(self, ciffilePath, parent=None):
        if not Path(ciffilePath).is_file():
            noteDiv = parent.addDiv(style='font-size:110%;color:red;')
            noteDiv.append(
                f"DNATCO did not generate the extended CIF file (expected path: {ciffilePath})."
                " No report could be generated.<br>"
                " Please check the log files for more information."
            )
            return None
        ciffile = gemmi.cif.read_file(ciffilePath)
        cif_data = {}

        cif_data['overall_confal_score'] = ciffile[0].find_value('_ndb_struct_ntc_overall.confal_score')
        cif_data['overall_confal_percentile'] = ciffile[0].find_value('_ndb_struct_ntc_overall.confal_percentile')
        cif_data['overall_num_classified'] = ciffile[0].find_value('_ndb_struct_ntc_overall.num_classified')
        cif_data['overall_num_unclassified'] = ciffile[0].find_value('_ndb_struct_ntc_overall.num_unclassified')
        cif_data['overall_num_unclassified_rmsd_close'] = ciffile[0].find_value(
            '_ndb_struct_ntc_overall.num_unclassified_rmsd_close'
        )

        table_step_id = ciffile[0].find(
                '_ndb_struct_ntc_step.',
                ['id',
                 'name',
                 'PDB_model_number',
                 'label_entity_id_1',
                 'label_asym_id_1',
                 'label_seq_id_1',
                 'label_comp_id_1',
                 'label_alt_id_1',
                 'label_entity_id_2',
                 'label_asym_id_2',
                 'label_seq_id_2',
                 'label_comp_id_2',
                 'label_alt_id_2',
                 'auth_asym_id_1',
                 'auth_seq_id_1',
                 'auth_asym_id_2',
                 'auth_seq_id_2',
                 'PDB_ins_code_1',
                 'PDB_ins_code_2']
            )
        cif_data['chain_label'] = table_step_id.find_column('label_asym_id_1')
        cif_data['chain_auth'] = table_step_id.find_column('auth_asym_id_1')
        cif_data['chain_display'] = [
            f"{label}"
            if label == auth
            else f"{label} (auth: {auth})" for label, auth in zip(cif_data['chain_label'], cif_data['chain_auth'])
            ]
        cif_data['comp_id_1'] = table_step_id.find_column('label_comp_id_1')
        cif_data['auth_seq_id_1'] = table_step_id.find_column('auth_seq_id_1')
        cif_data['comp_id_2'] = table_step_id.find_column('label_comp_id_2')
        cif_data['auth_seq_id_2'] = table_step_id.find_column('auth_seq_id_2')
        cif_data['steps'] = [f"{comp_id_1}{seq_id_1} {comp_id_2}{seq_id_2}" for comp_id_1, seq_id_1, comp_id_2, seq_id_2 in zip(
                  cif_data['comp_id_1'], cif_data['auth_seq_id_1'], cif_data['comp_id_2'], cif_data['auth_seq_id_2'])]

        table_step = ciffile[0].find(
            '_ndb_struct_ntc_step_summary.',
            ['step_id',
             'assigned_CANA',
             'assigned_NtC',
             'confal_score',
             'euclidean_distance_NtC_ideal',
             'cartesian_rmsd_closest_NtC_representative',
             'closest_CANA',
             'closest_NtC',
             'closest_step_golden']
            )
        cif_data['step_id'] = table_step.find_column('step_id')
        cif_data['assigned_CANA'] = table_step.find_column('assigned_CANA')
        cif_data['assigned_NtC'] = table_step.find_column('assigned_NtC')
        cif_data['confal_score'] = table_step.find_column('confal_score')
        # cif_data['euclidean_distance_NtC_ideal'] = table_step.find_column('euclidean_distance_NtC_ideal')
        cif_data['rmsd'] = table_step.find_column('cartesian_rmsd_closest_NtC_representative')
        cif_data['closest_NtC'] = table_step.find_column('closest_NtC')
        cif_data['closest_CANA'] = table_step.find_column('closest_CANA')
        cif_data['rmsd_NtC_assigned'] = [
            f"{float(rmsd):.2f} ({closest_NtC})"
            if closest_NtC != '.'
            else f"{float(rmsd):.2f}" for rmsd, closest_NtC in zip(cif_data['rmsd'], cif_data['closest_NtC'])
        ]
        cif_data['idx_rmsd_below_0p5'] = [
            i for i, rmsd in enumerate(cif_data['rmsd'])
            if rmsd not in ('.', '?') and float(rmsd) <= 0.5
        ]
        cif_data['idx_rmsd_between_0p5_and_1'] = [
            i for i, rmsd in enumerate(cif_data['rmsd'])
            if rmsd not in ('.', '?') and 0.5 < float(rmsd) <= 1.0
        ]
        cif_data['idx_rmsd_higher_1'] = [
            i for i, rmsd in enumerate(cif_data['rmsd'])
            if rmsd not in ('.', '?') and float(rmsd) > 1.0
        ]
        cif_data["overall_num_rmsd_below_0p5"] = len(cif_data['idx_rmsd_below_0p5'])
        cif_data["overall_num_rmsd_between_0p5_and_1"] = len(cif_data['idx_rmsd_between_0p5_and_1'])
        cif_data["overall_num_rmsd_higher_1"]= len(cif_data['idx_rmsd_higher_1'])

        # Filter data based on CANA/NtC assignment and RMSD values
        cif_data['idx_improvables'] = []
        cif_data['idx_outliers'] = []
        for i, (ntc, rmsd) in enumerate(zip(cif_data['assigned_NtC'], cif_data['rmsd'])):
            if ntc == 'NANT':
                if rmsd not in ('.', '?') and rmsd.isnumeric() and float(rmsd) <= 0.5:
                    cif_data['idx_improvables'].append(i)
                else:
                    cif_data['idx_outliers'].append(i)
        cif_data['step_id_improvables'] = [cif_data['step_id'][i] for i in cif_data['idx_improvables']]
        cif_data['chain_display_improvables'] = [cif_data['chain_display'][i] for i in cif_data['idx_improvables']]
        cif_data['steps_improvables'] = [cif_data['steps'][i] for i in cif_data['idx_improvables']]
        cif_data['assigned_CANA_improvables'] = [cif_data['closest_CANA'][i] for i in cif_data['idx_improvables']]
        cif_data['assigned_NtC_improvables'] = [cif_data['closest_NtC'][i] for i in cif_data['idx_improvables']]
        cif_data['rmsd_NtC_improvables'] = [cif_data['rmsd'][i] for i in cif_data['idx_improvables']]

        cif_data['step_id_outliers'] = [cif_data['step_id'][i] for i in cif_data['idx_outliers']]
        cif_data['chain_display_outliers'] = [cif_data['chain_display'][i] for i in cif_data['idx_outliers']]
        cif_data['steps_outliers'] = [cif_data['steps'][i] for i in cif_data['idx_outliers']]
        cif_data['assigned_CANA_outliers'] = [ cif_data['closest_CANA'][i] for i in cif_data['idx_outliers']]
        cif_data['assigned_NtC_outliers'] = [cif_data['closest_NtC'][i] for i in cif_data['idx_outliers']]
        cif_data['rmsd_NtC_outliers'] = [cif_data['rmsd'][i] for i in cif_data['idx_outliers']]

        cif_data['assigned_CANA'] = [cana if cana != 'NAN' else '-' for cana in cif_data['assigned_CANA']]
        cif_data['assigned_NtC'] = [ntc if ntc != 'NANT' else '-' for ntc in cif_data['assigned_NtC']]

        return cif_data
