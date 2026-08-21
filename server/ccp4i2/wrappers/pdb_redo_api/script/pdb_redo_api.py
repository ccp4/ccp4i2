import json
import os
import shutil
import sys
import zipfile

from lxml import etree

from ccp4i2.config import credentials
from ccp4i2.core import CCP4Utils
from ccp4i2.core import CCP4ErrorHandling
from ccp4i2.core.CCP4ErrorHandling import SEVERITY_WARNING
from ccp4i2.core.CCP4PluginScript import CPluginScript
from ccp4i2.core.CCP4XtalData import CObsDataFile

from . import test_api

# Object path the missing-credential error is reported against. The CREDENTIAL_
# prefix is a contract with the frontend: the run dialog's quick-action table
# keys on it to offer a "Set PDB-REDO token..." button beside the message, so the
# error that blocks the run also carries its own fix.
CREDENTIAL_ERROR_PATH = 'CREDENTIAL_PDB_REDO'


def _flag(param):
    """Render a CBoolean parameter as the 0/1 int the PDB-REDO API expects.

    Three traps, all of which this avoids:

    * ``int(param)`` raises ``TypeError`` - CBoolean defines ``__bool__`` but
      not ``__int__``, so there is nothing for ``int()`` to call.
    * ``bool(param)`` conflates two questions: CBoolean.__bool__ is
      ``isSet(allowDefault=True) and value``, so it answers "set AND true"
      rather than "what is the value". Here we only want the value.
    * ``param.value`` may arrive as the string 'True'/'False' from some paths.
      ``.get()`` returns the coerced primitive, which is what the wire format
      wants - the same idiom arp_warp_classic uses for its CBoolean flags.
    """
    return int(bool(param.get()))


class pdb_redo_api(CPluginScript):

    TASKNAME = 'pdb_redo_api'
    PERFORMANCECLASS = 'CRefinementPerformance'

    def _token(self):
        """Resolve the PDB-REDO token pair, or None if not fully configured.

        Never falls back to a partial pair: signing a request with half a token
        produced an opaque HTTP 401 several minutes into a job, which is the
        failure mode the credential store exists to remove.
        """
        return credentials.get_credential('pdb_redo')

    def validity(self):
        """Cheap, polled validation — is a token configured at all?

        No network access here; this runs on every parameter edit. The live
        check belongs in runTimeValidity().
        """
        error = super(pdb_redo_api, self).validity()

        if self._token() is None:
            error.append(
                klass=self.TASKNAME, code=200,
                details=(
                    'PDB-REDO needs an API token before it can submit a job. '
                    'Tokens are free from pdb-redo.eu/token; your token secret '
                    'stays on this machine.'
                ),
                name=f'{self.TASKNAME}.container.inputData.{CREDENTIAL_ERROR_PATH}',
                severity=CCP4ErrorHandling.SEVERITY_ERROR,
            )

        return error

    def runTimeValidity(self):
        """Pre-flight validation: does PDB-REDO actually accept this token?

        Runs once at submission, so it can afford a round trip. Tokens can be
        revoked at either end, and expiry is not readable from the API, so the
        only way to know a token still works is to ask.
        """
        error = super(pdb_redo_api, self).runTimeValidity()
        if error.maxSeverity() >= CCP4ErrorHandling.SEVERITY_ERROR:
            return error  # already failing; do not spend a round trip

        ok, message = credentials.validate_credential('pdb_redo')
        if not ok:
            error.append(
                klass=self.TASKNAME, code=201,
                details=message,
                name=f'{self.TASKNAME}.container.inputData.{CREDENTIAL_ERROR_PATH}',
                severity=CCP4ErrorHandling.SEVERITY_ERROR,
            )
        return error

    def startProcess(self):
        print("pdb_redo_api.startProcess")
        inp = self.container.inputData
        xyzin = str( inp.XYZIN.fullPath )
        print(self.hklin, xyzin)

        token = self._token()
        if token is None:
            print(
                "ERROR: no PDB-REDO API token is configured. Set one in "
                "Preferences -> Credentials, or via the PDB_REDO_TOKEN_ID and "
                "PDB_REDO_TOKEN_SECRET environment variables."
            )
            return CPluginScript.FAILED
        token_id = token['token_id']
        token_secret = token['token_secret']

        sequence=None
        restraints=None

        if inp.SEQIN.isSet():
            sequence = str(inp.SEQIN.fullPath)

        if inp.DICT.isSet():
            restraints = str(inp.DICT.fullPath)

        cpar = self.container.controlParameters
        params = {
            'paired': _flag(cpar.PAIRED),
            'noloops': _flag(cpar.NOLOOPS),
            'nopepflip': _flag(cpar.NOPEPFLIP),
            'noscbuild': _flag(cpar.NOSCBUILD),
            'nocentrifuge': _flag(cpar.NOCENTRIFUGE),
            'nosugarbuild': _flag(cpar.NOSUGARBUILD),
            'norebuild': _flag(cpar.NOREBUILD),
            'newmodel': _flag(cpar.NEWMODEL),
            'isotropic': _flag(cpar.ISOTROPIC),
            'anisotropic': _flag(cpar.ANISOTROPIC),
            'notls': _flag(cpar.NOTLS),
            'tighter': _flag(cpar.TIGHTER),
            'looser': _flag(cpar.LOOSER),
        }

        with open(self.makeFileName("PROGRAMXML"),"w") as programXMLFile:
            xmlStructure = etree.Element("pdb_redo_api")
            CCP4Utils.writeXML(programXMLFile,etree.tostring(xmlStructure))

        print("Submit PDB-REDO job")
        self.pdb_redo_job_id = test_api.submit(xyzin,self.hklin,token_id,token_secret,sequence=sequence,restraints=restraints,params=params)
        with open(self.makeFileName("PROGRAMXML"),"w") as programXMLFile:
            xmlStructure = etree.Element("pdb_redo_api")
            pdbRedoId = etree.SubElement(xmlStructure,'PDB_REDO_JOB_ID')
            pdbRedoId.text = str(self.pdb_redo_job_id)
            CCP4Utils.writeXML(programXMLFile,etree.tostring(xmlStructure))

        print("Wait for PDB-REDO job")
        try:
            test_api.monitor(self.pdb_redo_job_id, token_id, token_secret)

        except test_api.PDBRedoAuthError as err:
            # The token stopped working mid-run. Nothing to salvage (fetching
            # the results would fail the same way), and the fix is a UI action.
            return self._fail(
                203,
                "{err} Set a working token in Preferences -> Credentials, then "
                "re-run. Your PDB-REDO run may have completed regardless: see "
                "run {run_id} at {uri}/jobs".format(
                    err=err, run_id=self.pdb_redo_job_id, uri=test_api.PDBREDO_URI))

        except test_api.PDBRedoUnreachable as err:
            # We gave up waiting; this says nothing about the remote run, so
            # the important thing is to hand back the run id.
            return self._fail(
                204,
                "{err} Nothing has been lost at the PDB-REDO end: look for run "
                "{run_id} at {uri}/jobs and, once it has finished, re-run this "
                "task or download the results there.".format(
                    err=err, run_id=self.pdb_redo_job_id, uri=test_api.PDBREDO_URI))

        except test_api.PDBRedoJobStopped as err:
            # A real failure of the science at their end. The results zip
            # normally exists and holds the logs, so pull it in for diagnosis.
            salvaged = self._salvage_results(token_id, token_secret)
            detail = ("Its log files have been retrieved into this job's "
                      "directory." if salvaged else
                      "The results could not be retrieved for diagnosis.")
            return self._fail(
                205,
                "{err} {detail} PDB-REDO run {run_id} ({uri}/jobs)".format(
                    err=err, detail=detail, run_id=self.pdb_redo_job_id,
                    uri=test_api.PDBREDO_URI))

        except Exception as err:
            # Genuinely unexpected: report the class and message rather than
            # swallowing it, and still try to salvage anything downloadable.
            self._salvage_results(token_id, token_secret)
            return self._fail(
                206,
                "Unexpected error while waiting for PDB-REDO ({cls}: {err}). "
                "PDB-REDO run {run_id} may still be available at {uri}/jobs"
                .format(cls=err.__class__.__name__, err=err,
                        run_id=self.pdb_redo_job_id, uri=test_api.PDBREDO_URI))

        print("PDB-REDO job has finished, leave rest to processOutputFiles")

        return CPluginScript.SUCCEEDED

    def _fail(self, code, message):
        """Record a job failure that says what happened and what to do next.

        Both channels matter: appendErrorReport surfaces the message in the
        CCP4i2 UI, while the print keeps it in the job log next to whatever
        else was happening at the time.
        """
        print("ERROR:", message)
        sys.stdout.flush()
        self.appendErrorReport(code, message,
                               severity=CCP4ErrorHandling.SEVERITY_ERROR)
        return CPluginScript.FAILED

    def _salvage_results(self, token_id, token_secret):
        """Try to download and unpack the results zip after a failure.

        PDB-REDO writes its logs into the same zip as its results, so a run
        that failed scientifically still has the explanation in it. Best
        effort: this runs on an error path and must not raise.

        Returns:
            bool: True if a zip was retrieved and unpacked.
        """
        output_zip = os.path.join(self.getWorkDirectory(), "pdb_redo_results.zip")
        try:
            print("Fetching results after job failure, for diagnosis")
            sys.stdout.flush()
            test_api.do_fetch(self.pdb_redo_job_id, token_id, token_secret,
                              output_zip)
            print("Got", output_zip)
            sys.stdout.flush()

            with open(self.makeFileName("PROGRAMXML"), "w") as programXMLFile:
                xmlStructure = etree.Element("pdb_redo_api")
                pdbRedoId = etree.SubElement(xmlStructure, 'PDB_REDO_JOB_ID')
                pdbRedoId.text = str(self.pdb_redo_job_id)
                with zipfile.ZipFile(output_zip) as myzip:
                    for f in myzip.namelist():
                        myzip.extract(f, self.getWorkDirectory())
                        pdbRedoDir = etree.SubElement(xmlStructure,
                                                      'PDB_REDO_RESULTS_DIR')
                        pdbRedoDir.text = str(os.path.dirname(f))
                CCP4Utils.writeXML(programXMLFile, etree.tostring(xmlStructure))
            return True
        except Exception as err:
            print("Could not retrieve PDB-REDO results for diagnosis ({cls}: "
                  "{err})".format(cls=err.__class__.__name__, err=err))
            sys.stdout.flush()
            return False

    def processInputFiles(self):
        miniMtzs = [
            ["F_SIGF", CObsDataFile.CONTENT_FLAG_FMEAN],
            ["FREERFLAG", None],
        ]
        self.hklin, self.columns, error = self.makeHklin0(miniMtzs)
        if error.maxSeverity() > SEVERITY_WARNING:
            return CPluginScript.FAILED
        return CPluginScript.SUCCEEDED

    def processOutputFiles(self):
        token = self._token()
        if token is None:
            print("ERROR: PDB-REDO API token is no longer available")
            return CPluginScript.FAILED
        token_id = token['token_id']
        token_secret = token['token_secret']

        print("Extracting from zip"); sys.stdout.flush()
        output_zip = os.path.join(self.getWorkDirectory(),"pdb_redo_results.zip")
        test_api.do_fetch(self.pdb_redo_job_id,token_id,token_secret,output_zip)
        print("Got",output_zip); sys.stdout.flush()
        
        redoDir = None
        finalRefmacLog = None
        pdbRedoLog = None

        with zipfile.ZipFile(output_zip) as myzip:
            outputColumns = ['FWT,PHWT','DELFWT,PHDELWT']
            infolist = myzip.infolist()
            for finfo in infolist:
                if finfo.filename.endswith("_besttls.mtz"):
                    outputFiles = ['FPHIOUT_BESTTLS','DIFFPHIOUT_BESTTLS']
                    print("Extracting",finfo.filename); sys.stdout.flush()
                    myzip.extract(finfo,self.getWorkDirectory())
                    print("Extracted",finfo.filename); sys.stdout.flush()
                    hkloutFile= os.path.join(self.getWorkDirectory(), finfo.filename)
                    print("Set hkloutFile",hkloutFile); sys.stdout.flush()
                    print("Splitting..."); sys.stdout.flush()
                    self.splitHklout(outputFiles,outputColumns,infile=hkloutFile)
                    print("Split",finfo.filename); sys.stdout.flush()
                if finfo.filename.endswith("_final.mtz"):
                    outputFiles = ['FPHIOUT_FINAL','DIFFPHIOUT_FINAL']
                    print("Extracting",finfo.filename); sys.stdout.flush()
                    myzip.extract(finfo,self.getWorkDirectory())
                    print("Extracted",finfo.filename); sys.stdout.flush()
                    hkloutFile=os.path.join(self.getWorkDirectory(), finfo.filename)
                    print("Set hkloutFile",hkloutFile); sys.stdout.flush()
                    print("Splitting..."); sys.stdout.flush()
                    self.splitHklout(outputFiles,outputColumns,infile=hkloutFile)
                    print("Split",finfo.filename); sys.stdout.flush()
                if finfo.filename.endswith("_besttls.pdb"):
                    print("Extracting",finfo.filename); sys.stdout.flush()
                    myzip.extract(finfo,self.getWorkDirectory())
                    print("Extracted",finfo.filename); sys.stdout.flush()
                    outputPDB      = os.path.normpath(os.path.join(self.getWorkDirectory(),finfo.filename))
                    outputFilePath = os.path.normpath(os.path.join(self.getWorkDirectory(),os.path.basename(finfo.filename)))
                    shutil.copyfile(outputPDB, outputFilePath)
                    self.container.outputData.XYZOUT_BESTTLS.setFullPath(outputFilePath)
                if finfo.filename.endswith("_final.pdb"):
                    print("Extracting",finfo.filename); sys.stdout.flush()
                    myzip.extract(finfo,self.getWorkDirectory())
                    print("Extracted",finfo.filename); sys.stdout.flush()
                    outputPDB      = os.path.normpath(os.path.join(self.getWorkDirectory(),finfo.filename))
                    outputFilePath = os.path.normpath(os.path.join(self.getWorkDirectory(),os.path.basename(finfo.filename)))
                    shutil.copyfile(outputPDB, outputFilePath)
                    self.container.outputData.XYZOUT_FINAL.setFullPath(outputFilePath)
                if finfo.filename.endswith("_final.log"):
                    print("Extracting",finfo.filename); sys.stdout.flush()
                    myzip.extract(finfo,self.getWorkDirectory())
                    print("Extracted",finfo.filename); sys.stdout.flush()
                    finalRefmacLog = finfo.filename
                if finfo.filename.endswith(".log") and not finfo.filename.endswith("_final.log") and os.path.basename(finfo.filename) != "process.log":
                    print("Extracting log",finfo.filename); sys.stdout.flush()
                    myzip.extract(finfo,self.getWorkDirectory())
                    print("Extracted log",finfo.filename); sys.stdout.flush()
                    pdbRedoLog = finfo.filename
                if finfo.filename.endswith(".html"):
                    print("Extracting",finfo.filename); sys.stdout.flush()
                    myzip.extract(finfo,self.getWorkDirectory())
                    print("Extracted",finfo.filename); sys.stdout.flush()
                if finfo.filename.endswith(".json"):
                    print("Extracting",finfo.filename); sys.stdout.flush()
                    myzip.extract(finfo,self.getWorkDirectory())
                    print("Extracted",finfo.filename); sys.stdout.flush()
                    if os.path.basename(finfo.filename) == "data.json":
                        print("Reading",os.path.join(self.getWorkDirectory(),finfo.filename))
                        with open(os.path.join(self.getWorkDirectory(),finfo.filename)) as f:
                            print("Loading as JSON")
                            redoDir = os.path.dirname(finfo.filename)
                            j = json.load(f)
                            if "properties" in j and "RFIN" in j["properties"]:
                                print("Set performance R")
                                self.container.outputData.PERFORMANCE.RFactor.set(j["properties"]["RFIN"])
                            if "properties" in j and "RFFIN" in j["properties"]:
                                print("Set performance RFree")
                                self.container.outputData.PERFORMANCE.RFree.set(j["properties"]["RFFIN"])

        with open(self.makeFileName("PROGRAMXML"),"w") as programXMLFile:
            xmlStructure = etree.Element("pdb_redo_api")
            if "properties" in j:
                for k,v in j["properties"].items():
                    ele = etree.SubElement(xmlStructure,k)
                    ele.text = str(v)
            pdbRedoDir = etree.SubElement(xmlStructure,'PDB_REDO_RESULTS_DIR')
            pdbRedoDir.text = str(redoDir)
            pdbRedoId = etree.SubElement(xmlStructure,'PDB_REDO_JOB_ID')
            pdbRedoId.text = str(self.pdb_redo_job_id)
            if finalRefmacLog:
                finalRefmacLogEle = etree.SubElement(xmlStructure,'PDB_REDO_FINAL_REFMAC_LOG_FILE')
                finalRefmacLogEle.text = str(finalRefmacLog)
            if pdbRedoLog:
                pdbRedoLogEle = etree.SubElement(xmlStructure,'PDB_REDO_LOG_FILE')
                pdbRedoLogEle.text = str(pdbRedoLog)
            CCP4Utils.writeXML(programXMLFile,etree.tostring(xmlStructure,encoding='utf-8', pretty_print=True))

        return CPluginScript.SUCCEEDED
