"""
     validate_protein.py: CCP4 GUI 2 Project
     Copyright (C) 2022 William Rochira

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

## @package validate_protein
# This package computes various metrics used in validation and deposition,
# for example average B-factors by chain, ligands, etc. and Ramachandran statistics.

import os
import sys
import time
import shutil
from math import pi

# Ignore NumPy warnings
import warnings
warnings.filterwarnings('ignore')

import numpy as np
#from lxml import etree
from xml.etree import ElementTree as ET

try:
    from core import CCP4Utils
    from core.CCP4PluginScript import CPluginScript
except ImportError:
    if 'CCP4' not in os.environ:
        sys.exit('Error: CCP4 environment variable must be set')
    sys.path.append(os.path.join(os.environ['CCP4'], 'share', 'ccp4i2'))
    try:
        from core import CCP4Utils
        from core.CCP4PluginScript import CPluginScript
    except ImportError:
        sys.exit('Error: Failed to import CCP4 core modules')

from iris_validation.graphics import Panel
from iris_validation.metrics import metrics_model_series_from_files


class validate_protein(CPluginScript):
    TASKNAME = 'validate_protein'
    WHATNEXT = [ 'coot_rebuild' ]
    MAINTAINER = 'jon.agirre@york.ac.uk'

    def process(self):
        from core import CCP4XtalData
        log_string = ''
        self.xml_root = ET.Element('Validate_geometry_CCP4i2')


        F_SIGF_1,errReport = self.makeHklin([['F_SIGF_1',CCP4XtalData.CObsDataFile.CONTENT_FLAG_FMEAN]],hklin="hklin_1")
        F_SIGF_2,errReport = self.makeHklin([['F_SIGF_2',CCP4XtalData.CObsDataFile.CONTENT_FLAG_FMEAN]],hklin="hklin_2")

        #F_SIGF_1,errReport = self.container.inputData.F_SIGF_1.convert(targetContent=CCP4XtalData.CObsDataFile.CONTENT_FLAG_FMEAN)
        #F_SIGF_2,errReport = self.container.inputData.F_SIGF_2.convert(targetContent=CCP4XtalData.CObsDataFile.CONTENT_FLAG_FMEAN)

        self.latest_model_path = str(self.container.inputData.XYZIN_1)
        self.latest_reflections_path = str(F_SIGF_1)
        self.previous_model_path = str(self.container.inputData.XYZIN_2)
        self.previous_reflections_path = str(F_SIGF_2)

        log_string += _print_and_return('\n\n######### Calculating metrics using Iris #########\n')
        try:
            fn_log, fn_xml = self.calculate_iris_metrics()
            log_string += fn_log
            self.xml_root.append(fn_xml)
        except Exception as e:
            log_string += _print_and_return('UNHANDLED ERROR:\n%s' % e)

        if self.container.controlParameters.DO_IRIS:
            log_string += _print_and_return('\n\n######### Generating Iris report SVG #########\n')
            try:
                fn_log, fn_xml = self.generate_iris_report()
                log_string += fn_log
                self.xml_root.append(fn_xml)
            except Exception as e:
                log_string += _print_and_return('UNHANDLED ERROR:\n%s' % e)

        if self.container.controlParameters.DO_MOLPROBITY:
            log_string += _print_and_return('\n\n######### Compiling MolProbity data #########\n')
            try:
                fn_log, fn_xml = self.compile_molprobity_data()
                log_string += fn_log
                self.xml_root.append(fn_xml)
            except Exception as e:
                log_string += _print_and_return('UNHANDLED ERROR:\n%s' % e)

        if self.container.controlParameters.DO_BFACT:
            log_string += _print_and_return('\n\n######### Summarising B-factor data #########\n')
            try:
                fn_log, fn_xml = self.b_averages()
                log_string += fn_log
                self.xml_root.append(fn_xml)
            except Exception as e:
                log_string += _print_and_return('UNHANDLED ERROR:\n%s' % e)

        if self.container.controlParameters.DO_RAMA:
            log_string += _print_and_return('\n\n######### Generating Ramachandran plots #########\n')
            try:
                fn_log, fn_xml = self.ramachandran_maps()
                log_string += fn_log
                self.xml_root.append(fn_xml)
            except Exception as e:
                log_string += _print_and_return('UNHANDLED ERROR:\n%s' % e)

        with open(self.makeFileName('PROGRAMXML'), 'w') as xml_file:
            ET.indent(self.xml_root)
            xml_file.write(ET.tostring(self.xml_root).decode('utf8'))

        self.reportStatus(CPluginScript.SUCCEEDED)
        return CPluginScript.SUCCEEDED


    def calculate_iris_metrics(self):
        log_string = ''
        xml_root = ET.Element('Model_info')
        self.model_series = metrics_model_series_from_files(model_paths=(self.previous_model_path, self.latest_model_path),
                                                            reflections_paths=(self.previous_reflections_path, self.latest_reflections_path),
                                                            sequence_paths=(None,),
                                                            distpred_paths=(None,),
                                                            model_json_paths=(None,),
                                                            run_covariance=False,
                                                            run_molprobity=self.container.controlParameters.DO_MOLPROBITY,
                                                            multiprocessing=None)
        self.latest_model = self.model_series.metrics_models[-1]
        ET.SubElement(xml_root, 'Chain_count').text = str(self.latest_model.chain_count)
        return log_string, xml_root


    def generate_iris_report(self):
        log_string = ''
        xml_root = ET.Element('Iris')
        model_series_data = self.model_series.get_raw_data()
        panel = Panel(model_series_data)
        panel.dwg.attribs['style'] += ' margin-top: -20px;'
        panel_string = panel.dwg.tostring()
        ET.SubElement(xml_root, 'Panel_svg').text = panel_string
        return log_string, xml_root


    def compile_molprobity_data(self):
        log_string = ''
        xml_root = ET.Element('Molprobity')
        molprobity_data = self.latest_model.molprobity_data['model_wide']

        xml_summary = ET.SubElement(xml_root, 'Summary')
        ET.SubElement(xml_summary, 'Ramachandran_outliers').text = str(round(molprobity_data['summary']['ramachandran_outliers'], 2)) + '%'
        ET.SubElement(xml_summary, 'Ramachandran_favoured').text = str(round(molprobity_data['summary']['ramachandran_favoured'], 2)) + '%' 
        ET.SubElement(xml_summary, 'Rotamer_outliers').text = str(round(molprobity_data['summary']['rotamer_outliers'], 2)) + '%'
        ET.SubElement(xml_summary, 'CBeta_deviations').text = str(molprobity_data['summary']['cbeta_deviations'])
        ET.SubElement(xml_summary, 'Clashscore').text = str(round(molprobity_data['summary']['clashscore'], 2))
        ET.SubElement(xml_summary, 'RMS_bonds').text = str(round(molprobity_data['summary']['rms_bonds'], 4))
        ET.SubElement(xml_summary, 'RMS_angles').text = str(round(molprobity_data['summary']['rms_angles'], 2))
        ET.SubElement(xml_summary, 'Molprobity_score').text = str(round(molprobity_data['summary']['molprobity_score'], 2))

        for category, element_name in (('ramachandran', 'Ramachandran_outliers'),
                                       ('omega', 'Nonplanar_omegas'),
                                       ('rotamer', 'Rotamer_outliers'),
                                       ('c-beta', 'CBeta_outliers'),
                                       ('nqh_flips', 'Side_chain_flips')):
            category_root = ET.SubElement(xml_root, element_name)
            for row in molprobity_data['details'][category]:
                chain, seqnum, name, score = row
                if score is None:
                    score = 'N/A'
                else:
                    score = str(round(score, 3))
                item = ET.SubElement(category_root, 'Outlier', chain=chain, seqnum=seqnum, name=name, score=score)

        clash_root = ET.SubElement(xml_root, 'Clashes')
        for row in molprobity_data['details']['clash']:
            item = ET.SubElement(clash_root, 'Outlier', first_atom=row[0], second_atom=row[1], overlap=str(row[2]))

        return log_string, xml_root


    def b_averages(self):
        log_string = ''
        xml_root = ET.Element('B_factors')

        b_factor_list_names = ('all', 'amino_acids', 'main_chains', 'side_chains', 'non_amino_acids', 'waters', 'ligands', 'ions')
        for chain_id, chain in enumerate(self.latest_model.chains):
            for list_name, b_factor_list in zip(b_factor_list_names, chain.b_factor_lists()):
                mean, std = np.mean(b_factor_list), np.std(b_factor_list)
                n = len(b_factor_list)
                ET.SubElement(xml_root, list_name, chain=str(chain_id), mean=str(mean), std=str(std), n=str(n))
        for list_name, b_factor_list in zip(b_factor_list_names, self.latest_model.b_factor_lists()):
            mean, std = np.mean(b_factor_list), np.std(b_factor_list)
            n = len(b_factor_list)
            ET.SubElement(xml_root, list_name, chain='All', mean=str(mean), std=str(std), n=str(n))
        return log_string, xml_root


    def ramachandran_maps(self):
        log_string = ''
        xml_root = ET.Element('Ramachandran')

        img_dir_from = os.path.join(os.path.normpath(CCP4Utils.getCCP4I2Dir()),
                                    os.path.normpath('wrappers/validate_protein/script/img/'))
        img_dir_to = os.path.join(os.path.normpath(self.workDirectory), 'img')
        shutil.copytree(img_dir_from, img_dir_to)

        favoured_root = ET.SubElement(xml_root, 'Favoured')
        allowed_root = ET.SubElement(xml_root, 'Allowed')
        outliers_root = ET.SubElement(xml_root, 'Outliers')
        n_residues, n_favoured, n_allowed, n_outliers, n_na = 0, 0, 0, 0, 0

        for chain in self.latest_model:
            for residue in chain:
                if not residue.is_aa:
                    continue
                n_residues += 1
                if None in (residue.phi, residue.psi):
                    n_na += 1
                    continue
                if residue.ramachandran_favoured:
                    residue_root = favoured_root
                    n_favoured += 1
                elif residue.ramachandran_allowed:
                    residue_root = allowed_root
                    n_allowed += 1
                else:
                    residue_root = outliers_root
                    n_outliers += 1
                residue_rama = ET.SubElement(residue_root, 'Residue', chain=str(chain.chain_id), seqnum=str(residue.sequence_number), type=residue.code)
                ET.SubElement(residue_rama, 'Phi').text = str(residue.phi * 180 / pi)
                ET.SubElement(residue_rama, 'Psi').text = str(residue.psi * 180 / pi)

        totals = ET.SubElement(xml_root, 'Totals')
        ET.SubElement(totals, 'Residues').text = str(n_residues)
        ET.SubElement(totals, 'Favoured').text = str(n_favoured)
        ET.SubElement(totals, 'Allowed').text = str(n_allowed)
        ET.SubElement(totals, 'Outliers').text = str(n_outliers)
        ET.SubElement(totals, 'NA').text = str(n_na)

        return log_string, xml_root


def _print_and_return(log_string):
    print(log_string)
    return log_string
