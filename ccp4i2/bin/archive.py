import argparse
import getpass
import http.cookiejar
import io
import itertools
import json
import mimetypes
import os
import re
import sys
import tempfile
import unicodedata
import urllib.parse
import urllib.request

from lxml import etree
import mimetools

from ..utils.startup import startProjectsManager


class MultiPartForm(object):
    """Accumulate the data to be used when posting a form."""

    def __init__(self):
        super(MultiPartForm, self).__init__()
        self.form_fields = []
        self.files = []
        self.boundary = mimetools.choose_boundary()
        return
    
    def get_content_type(self):
        return 'multipart/form-data; boundary=%s' % self.boundary

    def add_field(self, name, value):
        """Add a simple field to the form data."""
        self.form_fields.append((name, value))
        return

    def add_file(self, fieldname, filename, fileHandle, mimetype=None):
        """Add a file to be uploaded."""
        body = fileHandle.read()
        if mimetype is None:
            mimetype = mimetypes.guess_type(filename)[0] or 'application/octet-stream'
        self.files.append((fieldname, filename, mimetype, body))
        return

    def __str__(self):
        """Return a string representing the form data, including attached files."""
        # Build a list of lists, each containing "lines" of the
        # request.  Each part is separated by a boundary string.
        # Once the list is built, return a string where each
        # line is separated by '\r\n'.  
        parts = []
        part_boundary = '--' + self.boundary
        
        # Add the form fields
        parts.extend(
            [ part_boundary,
              'Content-Disposition: form-data; name="%s"' % name,
              '',
              value,
            ]
            for name, value in self.form_fields
            )
        
        # Add the files to upload
        parts.extend(
            [ part_boundary,
              'Content-Disposition: file; name="%s"; filename="%s"' % \
                 (field_name, filename),
              'Content-Type: %s' % content_type,
              '',
              body,
            ]
            for field_name, filename, content_type, body in self.files
            )
        
        # Flatten the list and add closing boundary marker,
        # then return CR+LF separated data
        flattened = list(itertools.chain(*parts))
        flattened.append('--' + self.boundary + '--')
        flattened.append('')
        return '\r\n'.join(flattened)


class DjangoMultiPartForm(MultiPartForm):
    def __init__(self, baseURL):
        super (DjangoMultiPartForm,self).__init__()
        self.baseURL = baseURL
        
    def get_request(self):
        # Build the request
        request = urllib.request.Request(self.baseURL)
        request.add_header('User-agent', 'PyMOTW (http://www.doughellmann.com/PyMOTW/)')
        body = str(self)
        request.add_header('Content-type', self.get_content_type())
        request.add_header('Content-length', len(body))
        request.add_data(body)
        return request

class DjangoSession (object):
    def __init__(self, *args, **kws):
        self.baseURL = args[0]
        self.username = args[1]
        self.password=args[2]
        self.baseURL = args[0]
        self.cookieLookup = {}
        handlers = self.createHandlers()
        self.opener = urllib.request.build_opener(*handlers)
        self.openSession()

    def createHandlers(self):
        handlers = [urllib.request.HTTPHandler(), urllib.request.HTTPSHandler()]
        handlers += [self.basicAuthHandler()]
        handlers += [self.cookieHandler()]
        return handlers

    def basicAuthHandler(self):
        password_mgr = urllib.request.HTTPPasswordMgrWithDefaultRealm()
        top_level_url = ""
        password_mgr.add_password(None, top_level_url, "username","PASSWORD")
        handler = urllib2.HTTPhandler = urllib.request.HTTPBasicAuthHandler(password_mgr)
        return handler

    def cookieHandler(self):
        self.cookies = http.cookiejar.LWPCookieJar()
        handler = urllib.request.HTTPCookieProcessor(self.cookies)
        return handler
    
    def openRequest(self, request):
        print(self.cookieLookup)
        response = self.opener.open(request)
        for cookie in self.cookies:
            self.cookieLookup[cookie.name] = cookie.value
        return response

    def openURL(self, url):
        return self.openRequest(url)

    def postURLWithValues(self, url, values):
        data = urllib.parse.urlencode(values)
        req = urllib.request.Request(url, data)
        response = self.openURL(req)
        return response

    def getURLWithValues(self, url, values):
        print('url',url)
        url_values = urllib.parse.urlencode(values)
        full_url = url + '?' + url_values
        print('full_url',full_url)
        req = urllib.request.Request(full_url)
        response = self.openURL(req)
        return response

    def openSession(self):
        response = self.openURL(self.baseURL)
        responseHTML=response.read()
        formURL = response.geturl()
        #Here check that this looks like a properly crafted login form
        parser = etree.HTMLParser()
        responseAsEtree = etree.parse(io.StringIO(responseHTML), parser)
        if len(responseAsEtree.xpath("//input[@name='username']")) == 0:
            raise Exception("No username field in login page retrieved from ", response.geturl())
        if len(responseAsEtree.xpath("//input[@name='password']")) == 0:
            raise Exception("No password field in login page retrieved from ", response.geturl())
        else:
            print('Retrieved login form with url', formURL)
        #Now post to this url, using the provided csrftoken and our stored credentials
        values = {'username':self.username,'password':self.password,'csrfmiddlewaretoken':self.cookieLookup['csrftoken'],next:self.baseURL}
        secondResponse = self.postURLWithValues(formURL, values)
        secondResponse.read()

    def multiPartForm(self, baseURL):
        form = DjangoMultiPartForm(baseURL)
        form.add_field('csrfmiddlewaretoken',self.cookieLookup['csrftoken'])
        return form


class CCP4i2DjangoSession(DjangoSession):
    def __init__(self, *args, **kws):
        super().__init__(*args, **kws)
        self.pm = self.myStartProjectsManager()

    def myStartProjectsManager(self):
        return startProjectsManager()

    def projectIdForName(self, projectName=None, strict=False):
        try:
            if strict: projectIdList = [projectTuple[0] for projectTuple in self.pm.db().listProjects() if projectName == projectTuple[1]]
            else: projectIdList = [projectTuple[0] for projectTuple in self.pm.db().listProjects() if projectName in projectTuple[1]]
            return projectIdList[0]
        except:
            return None

    def slugify(self, value):
        """
            Normalizes string, converts to lowercase, removes non-alpha characters,
            and converts spaces to hyphens.
            """
        value = unicodedata.normalize('NFKD', value.decode('unicode-escape')).encode('ascii', 'ignore')
        value = str(re.sub(r'[^\w\s-]', '', value).strip().lower())
        result = re.sub(r'[-\s]+', '-', value)
        return result

    def pushProject(self,projectName):
        projectId = self.projectIdForName(projectName=projectName, strict=False)
        
        projectList = self.pm.db().listProjects()
        projectNameList = [projectTuple[1] for projectTuple in projectList if projectId in projectTuple[0]]
        projectName = self.slugify(projectNameList[-1])
        
        tmpArchive = tempfile.NamedTemporaryFile(suffix='.ccp4_project.zip',delete=False)
        tmpArchive.close()
        fullPath=tmpArchive.name
        self.pm.compressProject(projectId,after=None,excludeI2files=False,fileName=fullPath,blocking=True)
        
        form = self.multiPartForm(self.baseURL+'/importProject')

        with open(fullPath,"r") as fileHandle:
            form.add_file('file', os.path.split(fullPath)[1], fileHandle)
        
        request = form.get_request()
        response = self.openRequest(request)
        print(response.read())

    def pushZip(self,zipName):
        
        form = self.multiPartForm(self.baseURL+'/importProject')

        with open(zipName,"r") as fileHandle:
            form.add_file('file', os.path.split(zipName)[1], fileHandle)
        
        request = form.get_request()
        response = self.openRequest(request)
        print(response.read())

    def fetchProject(self, projectName):
        response = self.getURLWithValues(self.baseURL+"/ProjectZipForProjName/"+projectName,{})
        tmpArchive = tempfile.NamedTemporaryFile(suffix='.ccp4_project.zip',delete=False)
        print(tmpArchive.name)
        CHUNK = 16 * 1024
        while True:
            chunk = response.read(CHUNK)
            if not chunk:
                break
            tmpArchive.write(chunk)
        tmpArchive.close()
        pm = self.myStartProjectsManager()
        from ..core import CCP4NonGuiProjectUtils
        CCP4NonGuiProjectUtils.CCP4NonGuiProjectUtils(tmpArchive.name)
        pm.db().commit()
        pm.db().close()
        return tmpArchive.name

    def queryDb(self, command, values=None):
        print('command',command,'values',values)
        response = self.getURLWithValues(self.baseURL+'?'+command, values)
        return json.loads(response.read())


def main():
    #First check for old syntax: a list contaiing password etc
    dashedArgs = [arg for arg in sys.argv[1:] if arg.startswith('-')]
    if len(sys.argv) == 5 and len(dashedArgs) == 4:
        parameterNamespace = Namespace(server=sys.argv[1], username=sys.argv[2], password=sys.argv[3],depositProject=[sys.argv[4]],fetchProject=[],depositZip=[])
    else:
        parser = argparse.ArgumentParser(description='Interact with CCP4i2 Archive.')
        parser.add_argument('-s','--server',help='URL of server e.g. http://myserver.local:8080/ManageCCP4i2Archive')
        parser.add_argument('-u','--username','--user',help='username of CCP4i2Archive e.g. djangouser. Defaults to login of current session')
        parser.add_argument('-p','--password','--pass',help='Password of CCP4i2Archive e.g. djangopassword')
        parser.add_argument('-f','--fetchProject',help='Fetch project from archive',action='append')
        parser.add_argument('-d','--depositProject',help='Deposit project in archive',action='append')
        parser.add_argument('-q','--queryDb',help='Perform database request',action='append')
        parser.add_argument('-z','--depositZip',help='Deposit zipped project in archive',action='append')
        parameterNamespace = parser.parse_args()

    print(parameterNamespace)
    if parameterNamespace.username is None: parameterNamespace.username = getpass.getuser()
    if parameterNamespace.password is None: parameterNamespace.password = getpass.getpass()
    ccp4i2DjangoSession = CCP4i2DjangoSession(parameterNamespace.server, parameterNamespace.username, parameterNamespace.password)

    if parameterNamespace.depositProject is not None:
        for depositProject in parameterNamespace.depositProject:
            ccp4i2DjangoSession.pushProject(depositProject)
    if parameterNamespace.fetchProject is not None:
        for fetchProject in parameterNamespace.fetchProject:
            ccp4i2DjangoSession.fetchProject(fetchProject)
    if parameterNamespace.depositZip is not None:
        for depositZip in parameterNamespace.depositZip:
            ccp4i2DjangoSession.pushZip(depositZip)
    if parameterNamespace.queryDb is not None:
        for queryDb in parameterNamespace.queryDb:
            tokens = queryDb.split('?')
            print('tokens',tokens)
            command = tokens[0]
            values = {}
            for valuePair in tokens[1:]:
                values[valuePair.split('=')[0]] = valuePair.split('=')[1] 
            print(ccp4i2DjangoSession.queryDb(command,values))
