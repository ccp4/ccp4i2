"""
    dnatco.py: CCP4 GUI Project
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

DNATCO (https://dnatco.datmos.org) assigns a dinucleotide conformer class
(NtC) to every dinucleotide step of a nucleic acid model, validates bond
lengths and angles against the NAVAL reference (ProSco probability scores),
and can write NtC-based external restraints for Refmac/Servalcat.

One run produces, with ``--prefix dnatco`` in the job directory:

    dnatco_extended.cif                       CIFOUT     model + NtC categories
    dnatco_angles_lengths_by_residue.json     JSONOUT    NAVAL validation
    dnatco_restraints_refmac.txt              RESTRAINTS (GENERATE_RESTRAINTS)
"""

import os
from pathlib import Path
import xml.etree.ElementTree as ET

from ccp4i2.config.program_discovery import resolve_program
from ccp4i2.core import CCP4ErrorHandling
from ccp4i2.core.CCP4PluginScript import CPluginScript
from . import dnatco_data


def dnatco_launcher_name():
    """The launcher CCP4 ships in its bin/: a shell script, or a batch file on Windows."""
    return "dnatco.bat" if os.name == "nt" else "dnatco.sh"


def dnatco_js_path():
    """$CCP4/dnatco/bin/dnatco.js when CCP4 is set up and ships DNATCO, else None."""
    ccp4 = os.environ.get("CCP4")
    if not ccp4:
        return None
    path = Path(ccp4) / "dnatco" / "bin" / "dnatco.js"
    return path if path.is_file() else None


def find_dnatco_command():
    """What to launch DNATCO with.

    Prefer the CCP4 launcher script, resolved through program discovery so a
    user can relocate it from Preferences -> Program locations. When no
    launcher resolves (a bundle without one, or Windows) but the bundled
    dnatco.js is present, run it with CCP4's ``node`` directly -- that is
    all the launcher does. Falls through to the launcher name so a missing
    install is reported as "dnatco.sh not found" rather than "node not
    found".
    """
    launcher = dnatco_launcher_name()
    if resolve_program(launcher):
        return launcher
    if dnatco_js_path() is not None and resolve_program("node"):
        return "node"
    return launcher


class dnatco(CPluginScript):

    TASKNAME = 'dnatco'
    TASKVERSION = 0.1
    TASKCOMMAND = find_dnatco_command()
    MAINTAINER = 'martin.maly@mrc-lmb.cam.ac.uk'
    WHATNEXT = ['servalcat_pipe']

    # --prefix: DNATCO names every output file <prefix>_<kind>.<ext>
    OUTPUT_PREFIX = 'dnatco'

    ERROR_CODES = {
        201: {'description': 'No output restraints file from DNATCO'},
        202: {'description': 'No output extended mmCIF file from DNATCO'},
        204: {'description': 'DNATCO failed to process the structure (see the log file)'},
        205: {'description': 'No NAVAL validation JSON file from DNATCO'},
        206: {'description': 'The input model contains no nucleic acid residues; DNATCO has nothing to validate'},
        207: {'description': 'Extended mmCIF from DNATCO carries no NtC assignment categories',
              'severity': CCP4ErrorHandling.SEVERITY_WARNING},
    }

    # DNATCO reports a structure it cannot process to its log and may still
    # exit cleanly; the base postProcessCheck() fails the job on this line.
    LOG_FAILURES = (
        (r'Failed to dnatcoify structure', 204, None),
    )

    def runTimeValidity(self):
        error = super().runTimeValidity()
        if error.maxSeverity() >= CCP4ErrorHandling.SEVERITY_ERROR:
            return error
        xyzin = self.container.inputData.XYZIN
        if xyzin.isSet() and not dnatco_data.model_has_nucleic_acid(str(xyzin.fullPath)):
            error.append(
                klass=self.TASKNAME, code=206,
                details=self.ERROR_CODES[206]['description'],
                name=f'{self.TASKNAME}.container.inputData.XYZIN',
                severity=CCP4ErrorHandling.SEVERITY_ERROR,
            )
        return error

    def makeCommandAndScript(self):
        if self.TASKCOMMAND == 'node':
            # No launcher script: drive the bundled dnatco.js ourselves.
            self.appendCommandLine([str(dnatco_js_path())])
        self.appendCommandLine(['--coords', str(self.container.inputData.XYZIN.fullPath)])
        self.appendCommandLine(['--outputDir', str(self.workDirectory)])
        self.appendCommandLine(['--prefix', self.OUTPUT_PREFIX])
        self.appendCommandLine(['--extendedCIF'])
        self.appendCommandLine(['--anglesLengthsByResidueJson'])
        params = self.container.controlParameters
        if bool(params.GENERATE_RESTRAINTS):
            self.appendCommandLine(['--refmacRestraints'])
            self.appendCommandLine(['--restraintsRmsd', str(float(params.MAX_RMSD))])
            self.appendCommandLine(['--restraintsSigmaFactor', str(float(params.RESTRAINTS_SIGMA))])
        return CPluginScript.SUCCEEDED

    def _output_path(self, kind):
        return Path(self.workDirectory) / f'{self.OUTPUT_PREFIX}_{kind}'

    def processOutputFiles(self):
        outputs = self.container.outputData
        params = self.container.controlParameters

        cif_path = self._output_path('extended.cif')
        if not cif_path.is_file():
            self.appendErrorReport(202, str(cif_path))
            return CPluginScript.FAILED
        outputs.CIFOUT.setFullPath(str(cif_path))
        outputs.CIFOUT.annotation.set('Model with DNATCO NtC assignments (mmCIF format)')

        json_path = self._output_path('angles_lengths_by_residue.json')
        if not json_path.is_file():
            self.appendErrorReport(205, str(json_path))
            return CPluginScript.FAILED
        outputs.JSONOUT.setFullPath(str(json_path))
        outputs.JSONOUT.annotation.set('NAVAL bond length and angle validation (JSON)')

        if bool(params.GENERATE_RESTRAINTS):
            restraints_path = self._output_path('restraints_refmac.txt')
            if not restraints_path.is_file():
                self.appendErrorReport(201, str(restraints_path))
                return CPluginScript.FAILED
            outputs.RESTRAINTS.setFullPath(str(restraints_path))
            outputs.RESTRAINTS.annotation.set('DNATCO NtC restraints for Refmac/Servalcat')

        self.writeProgramXml(cif_path)
        return CPluginScript.SUCCEEDED

    def writeProgramXml(self, cif_path):
        """program.xml: the overall NtC statistics, so a pipeline (or a later
        report) can summarise the run without re-reading the mmCIF."""
        root = ET.Element('DNATCO')
        ntc = dnatco_data.read_ntc_data(cif_path)
        if ntc is None:
            self.appendErrorReport(207, str(cif_path))
        else:
            overall = ET.SubElement(root, 'overall')
            for key, value in ntc['overall'].items():
                if value is not None:
                    ET.SubElement(overall, key).text = str(value)
            below, between, above = dnatco_data.rmsd_bins(ntc['steps'])
            ET.SubElement(overall, 'num_rmsd_below_0p5').text = str(below)
            ET.SubElement(overall, 'num_rmsd_between_0p5_and_1').text = str(between)
            ET.SubElement(overall, 'num_rmsd_above_1').text = str(above)
            ET.SubElement(overall, 'num_outliers').text = str(
                sum(1 for step in ntc['steps'] if step['is_outlier']))
            ET.SubElement(overall, 'num_improvable').text = str(
                sum(1 for step in ntc['steps'] if step['is_improvable']))
        ET.indent(root, space='  ')
        ET.ElementTree(root).write(self.makeFileName('PROGRAMXML'), encoding='utf-8', xml_declaration=True)
