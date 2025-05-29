"""
     validate_protein.py: CCP4i2 validation task
     Copyright (C) 2022-2024 William Rochira & Jon Agirre

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
import traceback

# Ignore NumPy warnings
import warnings
warnings.filterwarnings('ignore')

import numpy as np
from lxml import etree

try:
    from core import CCP4Utils, CCP4XtalData
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
import iris_validation

class validate_protein(CPluginScript):
    TASKNAME = 'validate_protein'
    WHATNEXT = [ 'coot_rebuild' ]
    MAINTAINER = 'jon.agirre@york.ac.uk'

    def process(self):
        log_string = ''
        self.xml_root = etree.Element('Validate_geometry_CCP4i2')

        print ("Using Iris installation: " + iris_validation.__file__.replace(f"/__init__.py", ""))

        self.latest_model_path = str(self.container.inputData.XYZIN_1)
        if self.container.inputData.F_SIGF_1.isSet() :
            self.latest_reflections_path, _ = self.makeHklin([['F_SIGF_1',
                                                             CCP4XtalData.CObsDataFile.CONTENT_FLAG_FMEAN]], 
                                                             hklin='F_SIGF_1')
        else: 
            self.latest_reflections_path = None
        
        if self.container.controlParameters.TWO_DATASETS :
            self.previous_model_path = str(self.container.inputData.XYZIN_2)
            if self.container.inputData.F_SIGF_2.isSet() :
                self.previous_reflections_path, _ = self.makeHklin([['F_SIGF_2',
                                                                CCP4XtalData.CObsDataFile.CONTENT_FLAG_FMEAN]], 
                                                                hklin='F_SIGF_2')
            else:
                self.previous_reflections_path = None
            log_string += _print_and_return('\n\nCalculating metrics for two datasets...\n\n')
        else :
            log_string += _print_and_return('\n\nCalculating metrics for one dataset...\n\n')
        log_string += _print_and_return('\n\n######### Calculating metrics using Iris #########\n')

        try:
            fn_log, fn_xml = self.calculate_iris_metrics()
            log_string += fn_log
            self.xml_root.append(fn_xml)
        except Exception as e:
            traceback.print_exc()
            log_string += _print_and_return('UNHANDLED ERROR:\n%s' % e)

        if self.container.controlParameters.DO_IRIS:
            log_string += _print_and_return('\n\n######### Generating Iris report SVG #########\n')
            try:
                fn_log, fn_xml = self.generate_iris_report()
                log_string += fn_log
                self.xml_root.append(fn_xml)
            except Exception as e:
                traceback.print_exc()
                log_string += _print_and_return('UNHANDLED ERROR:\n%s' % e)

        if self.container.controlParameters.DO_MOLPROBITY:
            log_string += _print_and_return('\n\n######### Compiling MolProbity data #########\n')
            try:
                fn_log, fn_xml = self.compile_molprobity_data()
                log_string += fn_log
                self.xml_root.append(fn_xml)
            except Exception as e:
                traceback.print_exc()
                log_string += _print_and_return('UNHANDLED ERROR:\n%s' % e)

        if self.container.controlParameters.DO_BFACT:
            log_string += _print_and_return('\n\n######### Summarising B-factor data #########\n')
            try:
                fn_log, fn_xml = self.b_averages()
                log_string += fn_log
                self.xml_root.append(fn_xml)
            except Exception as e:
                traceback.print_exc()
                log_string += _print_and_return('UNHANDLED ERROR:\n%s' % e)

        if self.container.controlParameters.DO_RAMA:
            log_string += _print_and_return('\n\n######### Generating Ramachandran plots #########\n')
            try:
                fn_log, fn_xml = self.ramachandran_maps()
                log_string += fn_log
                self.xml_root.append(fn_xml)
            except Exception as e:
                traceback.print_exc()
                log_string += _print_and_return('UNHANDLED ERROR:\n%s' % e)

        with open(self.makeFileName('PROGRAMXML'), 'w') as xml_file:
            xml_file.write(etree.tostring(self.xml_root, pretty_print=True).decode('utf8'))

        self.reportStatus(CPluginScript.SUCCEEDED)
        return CPluginScript.SUCCEEDED


    def calculate_iris_metrics(self):
        log_string = ''
        xml_root = etree.Element('Model_info')
        if self.container.controlParameters.TWO_DATASETS :
            print("PATHS\n",self.latest_model_path, self.previous_model_path,self.latest_reflections_path, self.previous_reflections_path)
            self.model_series = metrics_model_series_from_files(model_paths=(self.latest_model_path, self.previous_model_path),
                                                                reflections_paths=(self.latest_reflections_path, self.previous_reflections_path),
                                                                sequence_paths=(None,),
                                                                distpred_paths=(None,),
                                                                model_json_paths=(None,),
                                                                run_covariance=False,
                                                                calculate_rama_z=self.container.controlParameters.DO_TORTOIZE,
                                                                run_molprobity=self.container.controlParameters.DO_MOLPROBITY,
                                                                multiprocessing=False)
            self.latest_model = self.model_series.metrics_models[-1]
            etree.SubElement(xml_root, 'Chain_count').text = str(self.latest_model.chain_count)
            return log_string, xml_root
        else :
            print("PATHS\n",self.latest_model_path, self.latest_reflections_path )
            self.model_series = metrics_model_series_from_files(model_paths=(self.latest_model_path,),
                                                                reflections_paths=(self.latest_reflections_path,),
                                                                sequence_paths=(None,),
                                                                distpred_paths=(None,),
                                                                model_json_paths=(None,),
                                                                run_covariance=False,
                                                                calculate_rama_z=self.container.controlParameters.DO_TORTOIZE,
                                                                run_molprobity=self.container.controlParameters.DO_MOLPROBITY,
                                                                multiprocessing=False)
            self.latest_model = self.model_series.metrics_models[-1]
            etree.SubElement(xml_root, 'Chain_count').text = str(self.latest_model.chain_count)
            return log_string, xml_root

    def generate_iris_report(self):
        log_string = ''
        xml_root = etree.Element('Iris')
        model_series_data = self.model_series.get_raw_data()
        panel = Panel(model_series_data,
                      custom_labels={'Latest': self.container.inputData.NAME_1,
                                     'Previous': self.container.inputData.NAME_2
                                     }
                     )

        panel.dwg.attribs['style'] += ' margin-top: -20px;'
        panel_string = panel.dwg.tostring()
        etree.SubElement(xml_root, 'Panel_svg').text = panel_string
        return log_string, xml_root


    def compile_molprobity_data(self):
        log_string = ''
        xml_root = etree.Element('Molprobity')
        molprobity_data = self.latest_model.molprobity_data['model_wide']

        xml_summary = etree.SubElement(xml_root, 'Summary')
        etree.SubElement(xml_summary, 'Ramachandran_outliers').text = str(round(molprobity_data['summary']['ramachandran_outliers'], 2)) + '%'
        etree.SubElement(xml_summary, 'Ramachandran_favoured').text = str(round(molprobity_data['summary']['ramachandran_favoured'], 2)) + '%' 
        etree.SubElement(xml_summary, 'Rotamer_outliers').text = str(round(molprobity_data['summary']['rotamer_outliers'], 2)) + '%'
        etree.SubElement(xml_summary, 'CBeta_deviations').text = str(molprobity_data['summary']['cbeta_deviations'])
        etree.SubElement(xml_summary, 'Clashscore').text = str(round(molprobity_data['summary']['clashscore'], 2))
        etree.SubElement(xml_summary, 'RMS_bonds').text = str(round(molprobity_data['summary']['rms_bonds'], 4))
        etree.SubElement(xml_summary, 'RMS_angles').text = str(round(molprobity_data['summary']['rms_angles'], 2))
        etree.SubElement(xml_summary, 'Molprobity_score').text = str(round(molprobity_data['summary']['molprobity_score'], 2))

        for category, element_name in (('ramachandran', 'Ramachandran_outliers'),
                                       ('omega', 'Nonplanar_omegas'),
                                       ('rotamer', 'Rotamer_outliers'),
                                       ('c-beta', 'CBeta_outliers'),
                                       ('nqh_flips', 'Side_chain_flips')):
            category_root = etree.SubElement(xml_root, element_name)
            for row in molprobity_data['details'][category]:
                chain, seqnum, name, score = row
                if score is None:
                    score = 'N/A'
                else:
                    score = str(round(score, 3))
                item = etree.SubElement(category_root, 'Outlier', chain=chain, seqnum=seqnum, name=name, score=score)

        clash_root = etree.SubElement(xml_root, 'Clashes')
        for row in molprobity_data['details']['clash']:
            item = etree.SubElement(clash_root, 'Outlier', first_atom=row[0], second_atom=row[1], overlap=str(row[2]))

        return log_string, xml_root


    def b_averages(self):
        log_string = ''
        xml_root = etree.Element('B_factors')

        b_factor_list_names = ('all', 'amino_acids', 'main_chains', 'side_chains', 'non_amino_acids', 'waters', 'ligands', 'ions')
        for chain_id, chain in enumerate(self.latest_model.chains):
            for list_name, b_factor_list in zip(b_factor_list_names, chain.b_factor_lists()):
                mean, std = np.mean(b_factor_list), np.std(b_factor_list)
                n = len(b_factor_list)
                etree.SubElement(xml_root, list_name, chain=str(chain_id), mean=str(mean), std=str(std), n=str(n))
        for list_name, b_factor_list in zip(b_factor_list_names, self.latest_model.b_factor_lists()):
            mean, std = np.mean(b_factor_list), np.std(b_factor_list)
            n = len(b_factor_list)
            etree.SubElement(xml_root, list_name, chain='All', mean=str(mean), std=str(std), n=str(n))
        return log_string, xml_root


    def ramachandran_maps(self):
        log_string = ''
        xml_root = etree.Element('Ramachandran')

        img_dir_from = os.path.join(os.path.normpath(CCP4Utils.getCCP4I2Dir()),
                                    os.path.normpath('wrappers/validate_protein/script/img/'))
        img_dir_to = os.path.join(os.path.normpath(self.workDirectory), 'img')
        shutil.copytree(img_dir_from, img_dir_to)

        favoured_root = etree.SubElement(xml_root, 'Favoured')
        allowed_root = etree.SubElement(xml_root, 'Allowed')
        outliers_root = etree.SubElement(xml_root, 'Outliers')
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
                residue_rama = etree.SubElement(residue_root, 'Residue', chain=str(chain.chain_id), seqnum=str(residue.sequence_number), type=residue.code)
                etree.SubElement(residue_rama, 'Phi').text = str(residue.phi * 180 / pi)
                etree.SubElement(residue_rama, 'Psi').text = str(residue.psi * 180 / pi)

        totals = etree.SubElement(xml_root, 'Totals')
        etree.SubElement(totals, 'Residues').text = str(n_residues)
        etree.SubElement(totals, 'Favoured').text = str(n_favoured)
        etree.SubElement(totals, 'Allowed').text = str(n_allowed)
        etree.SubElement(totals, 'Outliers').text = str(n_outliers)
        etree.SubElement(totals, 'NA').text = str(n_na)

        return log_string, xml_root


def _print_and_return(log_string):
    print(log_string)
    return log_string
