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

from requests.auth import AuthBase
import hashlib
import hmac
import base64
import re
from urllib.parse import urlparse
from datetime import datetime, timezone

class PDBRedoAPIAuth(AuthBase):
    def __init__(self, token_id, token_secret):
        self.token_id = token_id
        self.token_secret = token_secret

    def hash_s(self, s):
        if (not isinstance(s, bytes)):
            s = bytes(s, 'utf-8')
        m = hashlib.sha256()
        m.update(s)
        h = m.digest()
        return bytes.decode(base64.b64encode(h), 'utf-8')

    def signed(self, secret, message):
        if (not isinstance(message, bytes)):
            message = bytes(message, 'utf-8')
        if (not isinstance(secret, bytes)):
            secret = bytes(secret, 'utf-8')
        return hmac.new(secret, message, digestmod=hashlib.sha256).digest()

    def signed_s(self, secret, message):
        signature = self.signed(secret, message)
        return bytes.decode(base64.b64encode(signature), 'utf-8')

    def __call__(self, r):
        body = r.body
        if (body == None):
            body = ''

        contentDigest = self.hash_s(body)

        uri = urlparse(r.url)
        path = uri.path
        host = uri.hostname
        
        if (re.match('(([0-9a-fA-F]{1,4}:){7,7}[0-9a-fA-F]{1,4}|([0-9a-fA-F]{1,4}:){1,7}:|([0-9a-fA-F]{1,4}:){1,6}:[0-9a-fA-F]{1,4}|([0-9a-fA-F]{1,4}:){1,5}(:[0-9a-fA-F]{1,4}){1,2}|([0-9a-fA-F]{1,4}:){1,4}(:[0-9a-fA-F]{1,4}){1,3}|([0-9a-fA-F]{1,4}:){1,3}(:[0-9a-fA-F]{1,4}){1,4}|([0-9a-fA-F]{1,4}:){1,2}(:[0-9a-fA-F]{1,4}){1,5}|[0-9a-fA-F]{1,4}:((:[0-9a-fA-F]{1,4}){1,6})|:((:[0-9a-fA-F]{1,4}){1,7}|:)|fe80:(:[0-9a-fA-F]{0,4}){0,4}%[0-9a-zA-Z]{1,}|::(ffff(:0{1,4}){0,1}:){0,1}((25[0-5]|(2[0-4]|1{0,1}[0-9]){0,1}[0-9])\.){3,3}(25[0-5]|(2[0-4]|1{0,1}[0-9]){0,1}[0-9])|([0-9a-fA-F]{1,4}:){1,4}:((25[0-5]|(2[0-4]|1{0,1}[0-9]){0,1}[0-9])\.){3,3}(25[0-5]|(2[0-4]|1{0,1}[0-9]){0,1}[0-9]))', host)):
            host = "[" + host + "]"

        if (uri.port != 80 and uri.port != 443 and uri.port != None):
            host = "{host}:{port}".format(host=host, port = uri.port)
            
        query = uri.query

        canonicalRequest = "\n".join([r.method, path, query, host, contentDigest])
        canonicalRequestHash = self.hash_s(canonicalRequest)

        now = datetime.utcnow()
        timestamp = now.strftime("%Y-%m-%dT%H:%M:%SZ")
        date = now.strftime("%Y%m%d")

        credential = "{token_id}/{date}/pdb-redo-api".format(token_id = self.token_id, date = date)
        stringToSign = "\n".join(["PDB-REDO-api", timestamp, credential, canonicalRequestHash])

        key = self.signed("PDB-REDO" + self.token_secret, date)
        signature = self.signed_s(key, stringToSign)

        r.headers.update({
            'X-PDB-REDO-Date': timestamp,
            'Authorization': "PDB-REDO-api Credential={credential},SignedHeaders=host;x-pdb-redo-content,Signature={signature}"
                .format(credential = credential, signature = signature)
        })

        return r
