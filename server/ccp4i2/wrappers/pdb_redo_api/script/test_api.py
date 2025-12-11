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

import json
import requests
import argparse
import time

from . import PDBRedoAPIAuth

# Due to a bug? in the server implementation, the :port is required here...
PDBREDO_URI = 'https://pdb-redo.eu'

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
    r = requests.post(PDBREDO_URI + "/api/run".format(token_id = token_id), auth = auth, files = files, data = {'parameters': json.dumps(_params)})
    r.raise_for_status()
    
    if (not r.ok):
        raise ValueError('Could not submit job to server: ' + r.text)

    run_id = r.json()['id']
    print("Job submitted with id", run_id)
    return run_id

def monitor(run_id,token_id,token_secret):
# Loop until job is done
    auth = PDBRedoAPIAuth.PDBRedoAPIAuth(token_id, token_secret)
    while(True):
        r = requests.get(PDBREDO_URI + "/api/run/{run_id}".format(token_id = token_id, run_id = run_id), auth = auth)
        status = r.json()['status']
    
        if (status == 'stopped'):
            raise ValueError('The job somehow failed after submitting')

        if (status == 'ended'):
            break

        print("Job status is", status)
        time.sleep(5)

def do_fetch(run_id,token_id,token_secret,output):
    auth = PDBRedoAPIAuth.PDBRedoAPIAuth(token_id, token_secret)
    r = requests.get(PDBREDO_URI + "/api/run/{run_id}/output/zipped".format(token_id = token_id, run_id = run_id), auth = auth)
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
