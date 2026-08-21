# SPDX-License-Identifier: BSD-2-Clause
# 
# Copyright (c) 2020 NKI/AVL, Netherlands Cancer Institute
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
# 
# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
# ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

# ---------------------------------------------------------------------------
# Local modifications (CCP4i2), on top of the NKI/AVL sample client above:
#   * typed exceptions, so the wrapper can tell "PDB-REDO says your run failed"
#     apart from "we lost contact" and "your token was rejected";
#   * monitor() tolerates transient errors and backs its polling interval off
#     from 5s to 60s, rather than hammering pdb-redo.eu every 5s for the hours a
#     run can take;
#   * request timeouts throughout, so a hung connection cannot wedge a job.
# ---------------------------------------------------------------------------

import json
import requests
import argparse
import time

from . import PDBRedoAPIAuth

# Due to a bug? in the server implementation, the :port is required here...
PDBREDO_URI = 'https://pdb-redo.eu'

# Per-request timeout. Generous, because the results zip can be large.
REQUEST_TIMEOUT = 60

# Polling: start responsive for short runs, then back off. A PDB-REDO run can
# take hours; at the original flat 5s that is ~720 requests an hour per job,
# which is neither kind to their server nor readable in a job log.
POLL_INITIAL_SECONDS = 5
POLL_MAX_SECONDS = 60
POLL_BACKOFF = 1.5

# Consecutive transient failures tolerated before we give up waiting. With the
# backoff above this spans roughly ten minutes of unreachability.
MAX_CONSECUTIVE_FAILURES = 12

# How often to say "still running" when nothing has changed, so a long wait
# does not look like a hang.
HEARTBEAT_SECONDS = 600


class PDBRedoError(RuntimeError):
    """Base class for PDB-REDO API failures.

    Carries the run id where one exists, because the single most useful thing
    to tell a user whose job failed locally is which remote run to go and look
    at -- it may well have completed regardless.
    """

    def __init__(self, message, run_id=None):
        super().__init__(message)
        self.run_id = run_id


class PDBRedoAuthError(PDBRedoError):
    """PDB-REDO rejected our credentials (401/403), e.g. a revoked token."""


class PDBRedoJobStopped(PDBRedoError):
    """PDB-REDO ran the job and it stopped/failed at their end.

    A real scientific failure, not an infrastructure one: the results zip
    usually exists and contains the logs that explain why.
    """


class PDBRedoUnreachable(PDBRedoError):
    """We lost contact with pdb-redo.eu while waiting.

    Says nothing about the remote run, which is very probably still going.
    """


def _raise_for_auth(response, run_id=None):
    """Map 401/403 onto PDBRedoAuthError; other statuses are left alone."""
    if response.status_code in (401, 403):
        raise PDBRedoAuthError(
            "PDB-REDO rejected the API token (HTTP {code}). It may have been "
            "revoked, or the token ID and secret may not match.".format(
                code=response.status_code),
            run_id=run_id)

def submit(xyzin,hklin,token_id,token_secret,sequence=None,restraints=None,params=None):

    files = {
        'pdb-file': open(xyzin, 'rb'),
        'mtz-file': open(hklin, 'rb')
    }
    if (restraints != None):
        files['restraints-file'] = open(restraints, 'rb')
    
    if (sequence != None):
        files['sequence-file'] = open(sequence, 'rb')

# Optional parameters, currently there's only one:
    _params = {
        'paired': False
    }
    if params:
       _params = params

# Create a new job/run
    auth = PDBRedoAPIAuth.PDBRedoAPIAuth(token_id, token_secret)
    print(auth)
    r = requests.post(PDBREDO_URI + "/api/run", auth = auth, files = files,
                      data = {'parameters': json.dumps(_params)},
                      timeout = REQUEST_TIMEOUT)
    _raise_for_auth(r)
    r.raise_for_status()
    
    if (not r.ok):
        raise ValueError('Could not submit job to server: ' + r.text)

    run_id = r.json()['id']
    print("Job submitted with id", run_id)
    return run_id

def monitor(run_id, token_id, token_secret,
            poll_initial=POLL_INITIAL_SECONDS,
            poll_max=POLL_MAX_SECONDS,
            max_consecutive_failures=MAX_CONSECUTIVE_FAILURES,
            heartbeat=HEARTBEAT_SECONDS):
    """Wait for a PDB-REDO run to finish, tolerating transient failures.

    Polls the run's status until it reaches 'ended', backing the interval off
    from ``poll_initial`` to ``poll_max``.

    A blip -- a 502, a dropped connection, a sleeping laptop -- must not fail a
    job whose remote run is proceeding normally, so transient errors are
    retried; only ``max_consecutive_failures`` of them in a row gives up. A
    rejected token, by contrast, will not fix itself and is raised immediately.

    Raises:
        PDBRedoAuthError:   the token was rejected mid-run.
        PDBRedoJobStopped:  PDB-REDO reports the run stopped/failed.
        PDBRedoUnreachable: we could not reach pdb-redo.eu for long enough to
            give up waiting. Says nothing about the remote run, which is
            probably still going.
    """
    auth = PDBRedoAPIAuth.PDBRedoAPIAuth(token_id, token_secret)
    url = PDBREDO_URI + "/api/run/{run_id}".format(run_id=run_id)

    delay = poll_initial
    failures = 0
    waited = 0.0
    since_heartbeat = 0.0
    last_status = None

    while True:
        try:
            r = requests.get(url, auth=auth, timeout=REQUEST_TIMEOUT)
            _raise_for_auth(r, run_id=run_id)
            r.raise_for_status()
            status = r.json()['status']
        except PDBRedoAuthError:
            raise
        except (requests.RequestException, ValueError, KeyError) as err:
            failures += 1
            if failures >= max_consecutive_failures:
                raise PDBRedoUnreachable(
                    "Lost contact with pdb-redo.eu: {n} consecutive failures "
                    "over {mins:.0f} minutes (last: {cls}). The PDB-REDO run "
                    "itself may still be running.".format(
                        n=failures, mins=waited / 60.0,
                        cls=err.__class__.__name__),
                    run_id=run_id) from err
            print("Could not reach pdb-redo.eu ({cls}); retry {n}/{m} in {d:.0f}s"
                  .format(cls=err.__class__.__name__, n=failures,
                          m=max_consecutive_failures, d=delay))
            time.sleep(delay)
            waited += delay
            delay = min(delay * POLL_BACKOFF, poll_max)
            continue

        failures = 0

        if (status == 'stopped'):
            raise PDBRedoJobStopped(
                "PDB-REDO run {run_id} stopped without completing. Its log "
                "files usually explain why.".format(run_id=run_id),
                run_id=run_id)

        if (status == 'ended'):
            return status

        if status != last_status:
            print("Job status is", status)
            last_status = status
            # A transition is a good moment to be responsive again.
            delay = poll_initial
            since_heartbeat = 0.0
        elif since_heartbeat >= heartbeat:
            print("Still {status} after {mins:.0f} minutes (PDB-REDO run "
                  "{run_id})".format(status=status, mins=waited / 60.0,
                                     run_id=run_id))
            since_heartbeat = 0.0

        time.sleep(delay)
        waited += delay
        since_heartbeat += delay
        delay = min(delay * POLL_BACKOFF, poll_max)



def do_fetch(run_id,token_id,token_secret,output):
    auth = PDBRedoAPIAuth.PDBRedoAPIAuth(token_id, token_secret)
    r = requests.get(PDBREDO_URI + "/api/run/{run_id}/output/zipped".format(run_id = run_id),
                     auth = auth, timeout = REQUEST_TIMEOUT)
    _raise_for_auth(r, run_id = run_id)
    r.raise_for_status()

    with open(output, 'wb') as f:
        f.write(r.content)

if __name__ == "__main__":

# collect arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--token-id', help='The token ID', required=True)
    parser.add_argument('--token-secret', help='The token secret', required=True)
    parser.add_argument('--xyzin', help='The coordinates file', required=True)
    parser.add_argument('--hklin', help='The diffraction data file', required=True)
    parser.add_argument('--paired', help='Do a paired refinement', action='store_true')

    args = parser.parse_args()

# The token id and secret for a session at PDB-REDO    
    token_id = args.token_id
    token_secret = args.token_secret

# The files to submit
    xyzin = args.xyzin
    hklin = args.hklin
    paired = args.paired
    params = {
        'paired': paired
    }

    submit(xyzin,hklin,token_id,token_secret,params=paired)
