"""
    dnatco_pipe.py: CCP4 GUI Project
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

Runs DNATCO on one nucleic acid model, or on two so the report can set them
side by side (typically the same structure before and after refinement).
NtC restraints, when requested, come from the second model when there is
one -- that is the model about to be refined further.
"""

import os
import shutil
import traceback
import xml.etree.ElementTree as ET

from ccp4i2.core import CCP4ErrorHandling
from ccp4i2.core.CCP4PluginScript import CPluginScript
from ccp4i2.wrappers.dnatco.script import dnatco_data
from ccp4i2.wrappers.dnatco.script.dnatco import dnatco as dnatco_wrapper


class dnatco_pipe(CPluginScript):

    TASKNAME = 'dnatco_pipe'
    TASKVERSION = 0.1
    MAINTAINER = 'martin.maly@mrc-lmb.cam.ac.uk'
    WHATNEXT = ['servalcat_pipe']
    # The dnatco wrapper is what actually runs. Declaring its program here gets
    # a missing DNATCO reported before the job starts, and lists it on
    # Preferences -> Program locations.
    AUXILIARY_PROGRAMS = (dnatco_wrapper.TASKCOMMAND,)

    ERROR_CODES = {
        201: {'description': 'No output restraints file from DNATCO'},
        202: {'description': 'No output extended mmCIF file from DNATCO'},
        203: {'description': 'DNATCO failed on the first structure model'},
        204: {'description': 'DNATCO failed on the second structure model'},
        205: {'description': 'No NAVAL validation JSON file from DNATCO'},
        206: {'description': 'A structure model contains no nucleic acid residues'},
        207: {'description': 'A second structure model was requested but none was given'},
        208: {'description': 'Failed to harvest output from a DNATCO sub-job'},
        209: {'description': 'Failed to create the dnatco sub-job'},
    }

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.dnatco1 = None
        self.dnatco2 = None

    # ------------------------------------------------------------------ inputs

    def compareTwoModels(self):
        return bool(self.container.controlParameters.TOGGLE_XYZIN2)

    def validity(self):
        error = super().validity()
        if self.compareTwoModels() and not self.container.inputData.XYZIN2.isSet():
            error.append(
                klass=self.TASKNAME, code=207,
                details=self.ERROR_CODES[207]['description'],
                name=f'{self.TASKNAME}.container.inputData.XYZIN2',
                severity=CCP4ErrorHandling.SEVERITY_ERROR,
            )
        return error

    def runTimeValidity(self):
        error = super().runTimeValidity()
        if error.maxSeverity() >= CCP4ErrorHandling.SEVERITY_ERROR:
            return error
        models = [('XYZIN1', self.container.inputData.XYZIN1)]
        if self.compareTwoModels():
            models.append(('XYZIN2', self.container.inputData.XYZIN2))
        for name, model in models:
            if model.isSet() and not dnatco_data.model_has_nucleic_acid(str(model.fullPath)):
                error.append(
                    klass=self.TASKNAME, code=206,
                    details=f'{self.ERROR_CODES[206]["description"]}: {model.fullPath}',
                    name=f'{self.TASKNAME}.container.inputData.{name}',
                    severity=CCP4ErrorHandling.SEVERITY_ERROR,
                )
        return error

    # --------------------------------------------------------------- execution

    def startProcess(self):
        self.dnatco1 = self._runDnatco(self.container.inputData.XYZIN1, code=203)
        if self.dnatco1 is None:
            return CPluginScript.FAILED
        if self.compareTwoModels():
            self.dnatco2 = self._runDnatco(self.container.inputData.XYZIN2, code=204)
            if self.dnatco2 is None:
                return CPluginScript.FAILED
        return CPluginScript.SUCCEEDED

    def _runDnatco(self, model, code):
        """Run the dnatco wrapper on one model synchronously.

        Returns the finished sub-plugin, or None after filing the failure.
        """
        try:
            plugin = self.makePluginObject('dnatco')
        except Exception as exc:
            self.appendErrorReport(209, f'{exc}\n{traceback.format_exc()}')
            return None
        try:
            plugin.container.inputData.XYZIN = model
            params = self.container.controlParameters
            plugin.container.controlParameters.GENERATE_RESTRAINTS.set(bool(params.GENERATE_RESTRAINTS))
            plugin.container.controlParameters.MAX_RMSD.set(float(params.MAX_RMSD))
            plugin.container.controlParameters.RESTRAINTS_SIGMA.set(float(params.RESTRAINTS_SIGMA))
            plugin.doAsync = False
            status = plugin.process()
        except Exception as exc:
            self.appendErrorReport(code, f'{exc}\n{traceback.format_exc()}')
            return None
        if status != CPluginScript.SUCCEEDED:
            self.appendErrorReport(code, str(model.fullPath))
            return None
        return plugin

    # ----------------------------------------------------------------- outputs

    def processOutputFiles(self):
        outputs = self.container.outputData
        harvested = []
        for label, plugin, cif_out, json_out in [
            ('model 1', self.dnatco1, outputs.CIFOUT1, outputs.JSONOUT1),
            ('model 2', self.dnatco2, outputs.CIFOUT2, outputs.JSONOUT2),
        ]:
            if plugin is None:
                continue
            sub_outputs = plugin.container.outputData
            if not self._harvest(sub_outputs.CIFOUT, cif_out,
                                 f'{label.capitalize()} with DNATCO NtC assignments (mmCIF format)'):
                self.appendErrorReport(202, str(sub_outputs.CIFOUT.fullPath))
                return CPluginScript.FAILED
            if not self._harvest(sub_outputs.JSONOUT, json_out,
                                 f'NAVAL bond length and angle validation of {label} (JSON)'):
                self.appendErrorReport(205, str(sub_outputs.JSONOUT.fullPath))
                return CPluginScript.FAILED
            harvested.append((label, plugin))

        if bool(self.container.controlParameters.GENERATE_RESTRAINTS):
            # Restraints from the model that will be refined next: the second
            # one when comparing, else the only one.
            label, plugin = harvested[-1]
            annotation = 'DNATCO NtC restraints for Refmac/Servalcat'
            if self.compareTwoModels():
                annotation += f' (from {label})'
            if not self._harvest(plugin.container.outputData.RESTRAINTS, outputs.RESTRAINTS, annotation):
                self.appendErrorReport(201, str(plugin.container.outputData.RESTRAINTS.fullPath))
                return CPluginScript.FAILED

        self._writeProgramXml(harvested)
        return CPluginScript.SUCCEEDED

    def _harvest(self, source, destination, annotation):
        """Copy a sub-job output file into this job under its own output name."""
        source_path = str(source.fullPath) if source.isSet() else ''
        if not source_path or not os.path.isfile(source_path):
            return False
        try:
            shutil.copyfile(source_path, str(destination.fullPath))
        except OSError as exc:
            self.appendErrorReport(208, f'{source_path} -> {destination.fullPath}: {exc}')
            return False
        destination.annotation.set(annotation)
        return True

    def _writeProgramXml(self, harvested):
        """program.xml: each sub-job's DNATCO summary, tagged with its model."""
        root = ET.Element('DNATCO_PIPE')
        for index, (label, plugin) in enumerate(harvested, start=1):
            try:
                summary = ET.parse(plugin.makeFileName('PROGRAMXML')).getroot()
            except (OSError, ET.ParseError):
                summary = ET.Element('DNATCO')
            summary.set('model', str(index))
            root.append(summary)
        ET.indent(root, space='  ')
        ET.ElementTree(root).write(self.makeFileName('PROGRAMXML'), encoding='utf-8', xml_declaration=True)
