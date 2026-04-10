"""
Martin Maly, MRC-LMB
"""

import os
import json
import gemmi
import xml.etree.ElementTree as ET
from math import degrees
from core.CCP4PluginScript import CPluginScript
from core.CCP4ErrorHandling import SEVERITY_WARNING
from core import CCP4Utils
from wrappers.servalcat.script.json2xml import json2xml
from . import json2restraints


class metalCoord(CPluginScript):

    TASKMODULE = 'wrappers'        # Where this plugin will appear on gui
    TASKNAME = 'metalCoord'        # Task name - should be same as class name
    TASKVERSION = 0.3              # Version of this plugin
    TASKCOMMAND = 'metalCoord'     # The command to run the executable
    MAINTAINER = 'martin.maly@mrc-lmb.cam.ac.uk'

    ERROR_CODES = {
        201: {'description': 'No output JSON file from metalCoord'},
        202: {'description': 'Log file does not report successful job completion', 'severity': SEVERITY_WARNING},
    }

    @staticmethod
    def _atom_site_to_seqid(site):
        icode = site.get('icode', '').strip()
        if not icode or icode == '.':
            return str(site['sequence'])
        return f"{site['sequence']}{icode}"

    def _get_atom(self, model, site, atom_key='name'):
        try:
            residue = model[site['chain']][self._atom_site_to_seqid(site)][site['residue']]
        except Exception:
            return None
        atom_name = site.get(atom_key, '')
        if not atom_name:
            return None

        altloc = site.get('altloc', '').strip()
        if altloc:
            atom = residue.find_atom(atom_name, altloc)
            if atom:
                return atom
        return residue.find_atom(atom_name, '*')

    @staticmethod
    def _image_position(cell, ref_pos, atom_pos, symmetry):
        try:
            symm_idx = int(symmetry)
        except Exception:
            symm_idx = 0
        if symm_idx <= 0:
            return atom_pos
        try:
            return cell.find_nearest_pbc_position(ref_pos, atom_pos, symm_idx)
        except Exception:
            return atom_pos

    @staticmethod
    def _repeat_for_multiple_ideal_values(model_value, ideal_values):
        if isinstance(ideal_values, (list, tuple)) and len(ideal_values) > 1:
            return [model_value for _ in ideal_values]
        return model_value


    def _add_model_geometry_to_json(self, output_json_stats):

        structure = gemmi.read_structure( str(self.container.inputData.XYZIN.fullPath))
        model = structure[0]

        for site in output_json_stats:
            if not isinstance(site, dict):
                continue
            metal_atom = self._get_atom(model, site, atom_key='metal')
            if not metal_atom:
                continue

            for ligand_class in site.get('ligands', []):
                # Distances reported in "base" and "pdb"
                for entry in ligand_class.get('base', []) + ligand_class.get('pdb', []):
                    ligand_site = entry.get('ligand', {})
                    ligand_atom = self._get_atom(model, ligand_site, atom_key='name')
                    if not ligand_atom:
                        entry['distance_model'] = None
                    else:
                        ligand_pos = self._image_position(
                            cell=structure.cell,
                            ref_pos=metal_atom.pos,
                            atom_pos=ligand_atom.pos,
                            symmetry=ligand_site.get('symmetry', 0))
                        entry['distance_model'] = round(metal_atom.pos.dist(ligand_pos), 2)
                    entry['distance_model'] = self._repeat_for_multiple_ideal_values(
                        model_value=entry.get('distance_model'),
                        ideal_values=entry.get('distance'))

                # Angles listed for ligand pairs
                for angle_entry in ligand_class.get('angles', []):
                    ligand1_site = angle_entry.get('ligand1', {})
                    ligand2_site = angle_entry.get('ligand2', {})
                    ligand1_atom = self._get_atom(model, ligand1_site, atom_key='name')
                    ligand2_atom = self._get_atom(model, ligand2_site, atom_key='name')
                    if not ligand1_atom or not ligand2_atom:
                        angle_entry['angle_model'] = None
                    else:
                        ligand1_pos = self._image_position(
                            cell=structure.cell,
                            ref_pos=metal_atom.pos,
                            atom_pos=ligand1_atom.pos,
                            symmetry=ligand1_site.get('symmetry', 0))
                        ligand2_pos = self._image_position(
                            cell=structure.cell,
                            ref_pos=metal_atom.pos,
                            atom_pos=ligand2_atom.pos,
                            symmetry=ligand2_site.get('symmetry', 0))
                        angle_model = degrees(gemmi.calculate_angle(ligand1_pos, metal_atom.pos, ligand2_pos))
                        angle_entry['angle_model'] = None if angle_model is None else round(angle_model, 2)
                    angle_entry['angle_model'] = self._repeat_for_multiple_ideal_values(
                        model_value=angle_entry.get('angle_model'),
                        ideal_values=angle_entry.get('angle'))


    def makeCommandAndScript(self):
        self.appendCommandLine('--no-progress')
        self.appendCommandLine('stats')
        self.appendCommandLine(['-p', self.container.inputData.XYZIN.fullPath])
        if self.container.controlParameters.MAXIMUM_COORDINATION_NUMBER.isSet():
            if str(self.container.controlParameters.MAXIMUM_COORDINATION_NUMBER) != "auto":
                self.appendCommandLine(['-c', self.container.controlParameters.MAXIMUM_COORDINATION_NUMBER])
                option = f"COORD{int(str(self.container.controlParameters.MAXIMUM_COORDINATION_NUMBER)):02d}"
                if hasattr(self.container.coordination, option):
                    if getattr(self.container.coordination, option) != "auto":
                        self.appendCommandLine(['--cl', getattr(self.container.coordination, option)])
        if self.container.controlParameters.MINIMUM_SAMPLE_SIZE.isSet():
            self.appendCommandLine(['-m', self.container.controlParameters.MINIMUM_SAMPLE_SIZE])
        if self.container.controlParameters.DISTANCE_THRESHOLD.isSet():
            self.appendCommandLine(['-d', self.container.controlParameters.DISTANCE_THRESHOLD])
        if self.container.controlParameters.PROCRUSTES_DISTANCE_THRESHOLD.isSet():
            self.appendCommandLine(['-t', self.container.controlParameters.PROCRUSTES_DISTANCE_THRESHOLD])
        if self.container.controlParameters.IDEAL_ANGLES:
            self.appendCommandLine('--ideal-angles')
        if self.container.controlParameters.SIMPLE:
            self.appendCommandLine('--simple')
        if self.container.controlParameters.USE_PDB:
            self.appendCommandLine('--use-pdb')
        self.appendCommandLine(['-l', self.container.inputData.LIGAND_CODE])
        outputJsonFilename = f"{self.container.inputData.LIGAND_CODE}.json"
        self.outputJsonPath = os.path.join(self.getWorkDirectory(), outputJsonFilename)
        self.appendCommandLine(['-o', outputJsonFilename])

    def processOutputFiles(self):
        # sanity check that metalCoord finished successfully - from .json.status.json
        outputJsonStatusFilename = f"{self.container.inputData.LIGAND_CODE}.json.status.json"
        outputJsonStatusPath = os.path.join(self.getWorkDirectory(), outputJsonStatusFilename)
        with open(outputJsonStatusPath, encoding="utf-8") as file:
            data = json.load(file)
        if data['status'] != "Success":
            self.appendErrorReport(202)

        # sanity check that metalCoord finished successfully - from log
        for line in self.logFileText().splitlines():
            if "Report written" in line:
                break
        else:
            self.appendErrorReport(202)

        # Load JSON file content
        if os.path.isfile(self.outputJsonPath):
            self.container.outputData.JSON.setFullPath(self.outputJsonPath)
            self.container.outputData.JSON.annotation = f'Full analysis for monomer {self.container.inputData.LIGAND_CODE}'
            with open(self.outputJsonPath, encoding="utf-8") as outputJsonFile:
                outputJsonText = outputJsonFile.read()
        else:
            self.appendErrorReport(201, str(self.container.outputData.JSON))
            return CPluginScript.FAILED
        outputJsonStats = json.loads(outputJsonText)

        # Add distances and angles from input model
        try:
            if isinstance(outputJsonStats, list) and len(outputJsonStats) > 0:
                self._add_model_geometry_to_json(outputJsonStats)
                with open(self.outputJsonPath, "w", encoding="utf-8") as outputJsonFile:
                    json.dump(outputJsonStats, outputJsonFile, indent=4)
        except Exception as e:
            print(f"Warning: Failed to add model geometry to JSON stats: {e}")

        outputJsonText = json.dumps(outputJsonStats)
        outputRestraintsPrefix = f"{self.container.inputData.LIGAND_CODE}_restraints"
        outputRestraintsFilename = f"{outputRestraintsPrefix}.txt"
        outputRestraintsMmcifFilename = f"{outputRestraintsPrefix}.mmcif"
        outputRestraintsPathPrefix = os.path.join(self.getWorkDirectory(), outputRestraintsPrefix)
        outputRestraintsPath = os.path.join(self.getWorkDirectory(), outputRestraintsFilename)
        outputRestraintsMmcifPath = os.path.join(self.getWorkDirectory(), outputRestraintsMmcifFilename)
        if self.container.controlParameters.SAVE_PDBMMCIF:
            stPath = str(self.container.inputData.XYZIN.fullPath)
        else:
            stPath = None

        # Convert JSON to external restraint keywords
        json2restraints.main(
            jsonPaths=[self.outputJsonPath],
            stPath=stPath,
            outputPrefix=outputRestraintsPathPrefix,
            jsonEquivalentsPath=None,
            keep_links=bool(self.container.controlParameters.KEEP_LINKS))
        if os.path.isfile(outputRestraintsPath):
            self.container.outputData.RESTRAINTS.setFullPath(outputRestraintsPath)
            self.container.outputData.RESTRAINTS.annotation = f'Restraints for {self.container.inputData.LIGAND_CODE}'
        if self.container.controlParameters.SAVE_PDBMMCIF:
            if os.path.isfile(outputRestraintsMmcifPath) and stPath:
                self.container.outputData.XYZOUT.setFullPath(outputRestraintsMmcifPath)
                self.container.outputData.XYZOUT.annotation = 'Structure model with links from MetalCoord (mmCIF format)'

        # Convert JSON to program.xml for i2 report
        xmlText = json2xml(list(outputJsonStats), tag_name_subroot="site")
        xmlroot = ET.fromstringlist(["<METALCOORD>", xmlText, "</METALCOORD>"])
        ET.indent(xmlroot, space="\t", level=0)
        with open(self.makeFileName('PROGRAMXML'), 'w', encoding="utf-8") as programXML:
            CCP4Utils.writeXML(programXML, ET.tostring(xmlroot))

        return CPluginScript.SUCCEEDED
