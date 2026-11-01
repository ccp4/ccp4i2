"""Submit to OneDep: create a wwPDB deposition from deposition-ready files.

The task takes the model mmCIF and reflections mmCIF the deposition-preparation
task wrote, creates a OneDep deposition for the depositor (or replaces the files
on an existing one), waits for OneDep to process the upload, and records the
deposition id. Everything descriptive (title, authors, sources, citation) is
read by OneDep from the model file, so there is nothing to type here beyond
who is depositing. The depositor finishes the remaining pages in the OneDep
web interface.

The API key is a credential (config/credentials.py), never a task parameter.
"""

import os
import shutil
import sys

from lxml import etree

from ccp4i2.config import credentials
from ccp4i2.core import CCP4ErrorHandling, CCP4Utils
from ccp4i2.core.CCP4PluginScript import CPluginScript

# Object path the missing-credential error is reported against. The CREDENTIAL_
# prefix is a contract with the frontend's run dialog, which offers a
# "Set token..." button beside the message.
CREDENTIAL_ERROR_PATH = 'CREDENTIAL_ONEDEP'


def _text(param):
    value = param.get() if hasattr(param, 'get') else param
    return '' if value is None else str(value).strip()


def _list(param):
    """The non-empty strings in a CList of CString, in order."""
    try:
        items = list(param)
    except TypeError:
        return []
    return [s for s in (_text(item) for item in items) if s]


class onedep_submit(CPluginScript):

    TASKNAME = 'onedep_submit'
    MAINTAINER = 'martin.noble@newcastle.ac.uk'
    ERROR_CODES = {
        200: {'description': 'OneDep needs an API key before it can submit'},
        201: {'description': 'OneDep rejected the API key'},
        202: {'description': 'Depositor details are incomplete or malformed'},
        203: {'description': 'The onedep_lib package is not installed'},
        204: {'description': 'OneDep rejected the API key during the run'},
        205: {'description': 'Lost contact with OneDep'},
        206: {'description': 'OneDep rejected the uploaded files'},
        207: {'description': 'Unexpected error while talking to OneDep'},
    }

    def __init__(self, *args, **kws):
        super().__init__(*args, **kws)
        self._ephemeral = None

    # ------------------------------------------------------------------
    # Credential
    # ------------------------------------------------------------------

    def _key(self):
        found = credentials.get_credential('onedep')
        return found['api_key'] if found else None

    def validity(self):
        """Cheap, polled: is a key configured, are the depositor details sane?"""
        error = super().validity()
        from . import onedep_client

        if self._key() is None:
            error.append(
                klass=self.TASKNAME, code=200,
                details=(
                    'Submitting to OneDep needs an API key. Generate one in the '
                    'wwPDB deposition interface (your account page, "Generate '
                    'Key") and set it here; it stays on this machine.'
                ),
                name=f'{self.TASKNAME}.container.inputData.{CREDENTIAL_ERROR_PATH}',
                severity=CCP4ErrorHandling.SEVERITY_ERROR,
            )
        cpar = self.container.controlParameters
        problems = []
        if not onedep_client.valid_email(_text(cpar.EMAIL)):
            problems.append('the depositor e-mail is not an e-mail address')
        orcids = _list(cpar.ORCIDS)
        if not orcids:
            problems.append('at least one ORCID iD is needed')
        for orcid in orcids:
            if not onedep_client.valid_orcid(orcid):
                problems.append(f'{orcid} is not an ORCID iD (0000-0000-0000-0000)')
        dep_id = _text(cpar.REPLACE_DEPOSITION_ID)
        if dep_id and not (dep_id.startswith('D_') and dep_id[2:].isdigit()):
            problems.append(f'{dep_id} is not a OneDep deposition id (D_ followed by digits)')
        if problems:
            error.append(
                klass=self.TASKNAME, code=202,
                details='; '.join(problems),
                name=f'{self.TASKNAME}.container.controlParameters.EMAIL',
                severity=CCP4ErrorHandling.SEVERITY_ERROR,
            )
        return error

    def runTimeValidity(self):
        """Once, at submission: is the library here, does OneDep accept the key?"""
        error = super().runTimeValidity()
        if error.maxSeverity() >= CCP4ErrorHandling.SEVERITY_ERROR:
            return error
        from . import onedep_client

        if not onedep_client.library_available():
            error.append(klass=self.TASKNAME, code=203, details=onedep_client.INSTALL_HINT,
                         name=f'{self.TASKNAME}.container.inputData.{CREDENTIAL_ERROR_PATH}',
                         severity=CCP4ErrorHandling.SEVERITY_ERROR)
            return error
        try:
            onedep_client.country_from_name(_text(self.container.controlParameters.COUNTRY))
        except ValueError as err:
            error.append(klass=self.TASKNAME, code=202, details=str(err),
                         name=f'{self.TASKNAME}.container.controlParameters.COUNTRY',
                         severity=CCP4ErrorHandling.SEVERITY_ERROR)
            return error
        ok, message = credentials.validate_credential('onedep')
        if not ok:
            error.append(klass=self.TASKNAME, code=201, details=message,
                         name=f'{self.TASKNAME}.container.inputData.{CREDENTIAL_ERROR_PATH}',
                         severity=CCP4ErrorHandling.SEVERITY_ERROR)
        return error

    # ------------------------------------------------------------------
    # Run
    # ------------------------------------------------------------------

    def _state_dir(self):
        """Where onedep_lib keeps its config and session for this run.

        Per-user and outside the project on a desktop; a throwaway directory
        on a hosted server, removed in processOutputFiles, so a rotated token
        never persists where other users' jobs run.
        """
        from . import onedep_client

        if os.environ.get('CCP4I2_LOCAL_SESSION_TOKEN'):
            return onedep_client.private_dir()
        self._ephemeral = onedep_client.ephemeral_dir()
        return self._ephemeral

    def _cleanup(self):
        if self._ephemeral is not None:
            shutil.rmtree(str(self._ephemeral), ignore_errors=True)
            self._ephemeral = None

    def _fail(self, code, message):
        print('ERROR:', message)
        sys.stdout.flush()
        self.appendErrorReport(code, message, severity=CCP4ErrorHandling.SEVERITY_ERROR)
        self._cleanup()
        return CPluginScript.FAILED

    def _write_xml(self, result=None, error=None):
        root = etree.Element('onedep_submit')
        if result:
            for tag in ('dep_id', 'site_url', 'session_id', 'status', 'details'):
                if result.get(tag):
                    etree.SubElement(root, tag.upper()).text = str(result[tag])
            files = etree.SubElement(root, 'FILES')
            for f in result.get('files', []):
                node = etree.SubElement(files, 'FILE')
                for key in ('name', 'type', 'remote_id'):
                    if f.get(key) is not None:
                        etree.SubElement(node, key.upper()).text = str(f[key])
        if error:
            etree.SubElement(root, 'ERROR').text = str(error)
        with open(self.makeFileName('PROGRAMXML'), 'w') as handle:
            CCP4Utils.writeXML(handle, etree.tostring(root, pretty_print=True))

    def startProcess(self):
        from . import onedep_client

        key = self._key()
        if key is None:
            return self._fail(200, 'No OneDep API key is configured. Set one in '
                                   'Preferences -> Credentials, or via ONEDEP_API_KEY.')
        try:
            onedep_client.lib()
        except onedep_client.OneDepMissing as err:
            return self._fail(203, str(err))

        inp = self.container.inputData
        cpar = self.container.controlParameters
        model = str(inp.MODEL_CIF.fullPath)
        sf = str(inp.SF_CIF.fullPath)
        email = _text(cpar.EMAIL)
        orcids = _list(cpar.ORCIDS)
        country = _text(cpar.COUNTRY)
        experiment = _text(cpar.EXPERIMENT_TYPE) or 'xray'
        replace = _text(cpar.REPLACE_DEPOSITION_ID) or None

        self._write_xml({'files': [], 'dep_id': replace})
        base = self._state_dir()
        config = onedep_client.make_config(key, base)

        def log(line):
            print(line)
            sys.stdout.flush()

        try:
            result = onedep_client.submit(model, sf, email, orcids, country, config,
                                          experiment_type=experiment, replace_dep_id=replace,
                                          log=log)
        except onedep_client.OneDepAuthError as err:
            self._write_xml({'dep_id': err.dep_id, 'site_url': err.site_url}, error=str(err))
            return self._fail(204, f'{err} The key is rotated on every use: if it was '
                                   'used from another machine since it was set here, '
                                   'generate a new one in the OneDep interface.')
        except onedep_client.OneDepUnreachable as err:
            self._write_xml({'dep_id': err.dep_id, 'site_url': err.site_url}, error=str(err))
            where = f' Deposition {err.dep_id} exists at OneDep;' if err.dep_id else ''
            return self._fail(205, f'{err}{where} open it in the deposition interface '
                                   'to see how processing ended. Nothing is lost at the '
                                   'OneDep end.')
        except onedep_client.OneDepProcessingError as err:
            self._write_xml({'dep_id': err.dep_id, 'site_url': err.site_url}, error=str(err))
            return self._fail(206, f'{err} Deposition {err.dep_id} is open in the OneDep '
                                   'interface; fix the files and re-run this task with its '
                                   'id in "Replace files on existing deposition".')
        except Exception as err:
            self._write_xml(error=f'{err.__class__.__name__}: {err}')
            return self._fail(207, f'Unexpected error while talking to OneDep '
                                   f'({err.__class__.__name__}: {err}).')
        finally:
            self._remember_rotated_key(key, config)

        self._write_xml(result)
        self._result = result
        log(f'OneDep deposition {result["dep_id"]} processed: {result["details"]}')
        return CPluginScript.SUCCEEDED

    def _remember_rotated_key(self, old_key, config):
        """Store the key OneDep rotated to, if the store is writable.

        Best effort and silent on failure: on a hosted server the store is
        environment-driven and the user generates a new key next time, which
        the task description says.
        """
        from . import onedep_client

        try:
            new_key = onedep_client.rotated_key(config, old_key)
            if new_key:
                credentials.remember_rotation('onedep', {'api_key': new_key})
        except Exception:
            print('WARNING: could not store the rotated OneDep key; generate a new '
                  'key before the next submission.')

    def processOutputFiles(self):
        self._cleanup()
        return CPluginScript.SUCCEEDED
