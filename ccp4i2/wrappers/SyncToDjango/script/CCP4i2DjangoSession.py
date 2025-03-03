from __future__ import print_function

import itertools
import email.generator
import mimetypes
from io import StringIO
import urllib.request, urllib.parse
from html.parser import HTMLParser

class MultiPartForm(object):
    """Accumulate the data to be used when posting a form."""

    def __init__(self):
        super(MultiPartForm, self).__init__()
        self.form_fields = []
        self.files = []
        self.boundary = email.generator._make_boundary()
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
        request = urllib.request.urlopen(self.baseURL)
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
        handler = urllib.request.HTTPhandler = urllib.request.HTTPBasicAuthHandler(password_mgr)
        return handler

    def cookieHandler(self):
        import cookielib
        self.cookies = cookielib.LWPCookieJar()
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
        req = urllib.request.urlopen(url, data)
        response = self.openURL(req)
        return response

    def getURLWithValues(self, url, values):
        url_values = urllib.parse.urlencode(values)
        full_url = url + '?' + url_values
        print('Full_url is:',full_url)
        req = urllib.request.urlopen(full_url)
        response = self.openURL(req)
        return response

    def openSession(self):
        response = self.openURL(self.baseURL)
        responseHTML=response.read()
        formURL = response.geturl()
        #Here check that this looks like a properly crafted login form
        try:
            parser = MyParser()
        except:
            raise Exception("Unable to instantiate MyParser ", response.geturl())
        try:
            parser.feed(responseHTML)
        except:
            raise Exception("Unable to parse ", responseHTML)
        if len(parser.nameInputElements) == 0:
            raise Exception("No username field in login page retrieved from ", response.geturl())
        if len(parser.passwordInputElements) == 0:
            raise Exception("No password field in login page retrieved from ", response.geturl())
        else:
            print('Retrieved login form with url', formURL)
        #Now post to this url, using the provided csrftoken and our stored credentials
        values = {'username':self.username,'password':self.password,'csrfmiddlewaretoken':self.cookieLookup['csrftoken'],next:self.baseURL}
        secondResponse = self.postURLWithValues(formURL, values)
        html = secondResponse.read()
        #print html

    def multiPartForm(self, baseURL):
        form = DjangoMultiPartForm(baseURL)
        form.add_field('csrfmiddlewaretoken',self.cookieLookup['csrftoken'])
        return form

class MyParser(HTMLParser):
    
    nameInputElements = set()
    passwordInputElements = set()
    
    def handle_startendtag(self, tag, attrs):
        if tag == "input":
            for att in attrs:
                if att[0] == "name":
                    if att[1] == 'username': self.nameInputElements.add(att[1])
                    if att[1] == 'password': self.passwordInputElements.add(att[1])

    def handle_starttag(self, tag, attrs):
        self.handle_startendtag(tag, attrs)

class CCP4i2DjangoSession(DjangoSession):
    def __init__(self, *args, **kws):
        super(CCP4i2DjangoSession, self).__init__(*args, **kws)
        self._pm = kws.get('projectsManager', None)

    def projectsManager(self):
        if self._pm is None: self._pm = self.myStartProjectsManager()
        return self._pm

    def myStartProjectsManager(self):
        import os
        import CCP4Utils
        CCP4I2_TOP = CCP4Utils.getCCP4I2Dir()
        sys.path.append(os.path.join(CCP4I2_TOP,'utils'))
        from startup import setupEnvironment, setupPythonpath, startProjectsManager
        setupEnvironment(path=CCP4I2_TOP)
        setupPythonpath(top=CCP4I2_TOP,mode='qtgui')
        pm = startProjectsManager()
        return pm
    
    def projectIdForName(self, projectName=None, strict=False):
        try:
            if strict: projectIdList = [projectTuple[0] for projectTuple in self.projectsManager().db().listProjects() if projectName == projectTuple[1]]
            else: projectIdList = [projectTuple[0] for projectTuple in self.projectsManager().db().listProjects() if projectName in projectTuple[1]]
            return projectIdList[0]
        except:
            return None
    
    def slugify(self, value):
        """
            Normalizes string, converts to lowercase, removes non-alpha characters,
            and converts spaces to hyphens.
            """
        import unicodedata, re
        value = unicodedata.normalize('NFKD', value.decode('unicode-escape')).encode('ascii', 'ignore')
        value = unicode(re.sub(r'[^\w\s-]', '', value).strip().lower())
        result = re.sub(r'[-\s]+', '-', value)
        return result

    def pushProjectWithName(self,projectName):
        projectId = self.projectIdForName(projectName=projectName, strict=False)
        self.pushProjectWithId(projectId)
        
    def pushProjectWithId(self, projectId, jobList=None):
        from core import CCP4NonGuiProjectUtils
        
        projectList = self.projectsManager().db().listProjects()
        projectNameList = [projectTuple[1] for projectTuple in projectList if projectId in projectTuple[0]]
        projectName = self.slugify(projectNameList[-1])

        usingStringIO = True
        if usingStringIO:
            zippedStringIO = StringIO() 
            self.projectsManager().compressProject(projectId,after=None,jobList=jobList, excludeI2files=False,fileName=zippedStringIO,blocking=True)
            form = self.multiPartForm(self.baseURL+'/importProject')
            zippedStringIO.seek(0)
            form.add_file('file', str(projectName+".ccp4_project.zip"), zippedStringIO)
        
        else:
            import tempfile
            with tempfile.NamedTemporaryFile(suffix='.ccp4_project.zip',delete=False) as tmpArchive:
                self.projectsManager().compressProject(projectId,after=None,jobList=jobList, excludeI2files=False,fileName=tmpArchive,blocking=True)
                form = self.multiPartForm(self.baseURL+'/importProject')
                tmpArchive.flush()
                import os
                os.fsync(tmpArchive.fileno())
                tmpArchive.seek(0)
                form.add_file('file', str(projectName+".ccp4_project.zip"), tmpArchive)
        
        request = form.get_request()
        response = self.openRequest(request)
        print(response.read())

    def pushZip(self,zipName):
        form = self.multiPartForm(self.baseURL+'/importProject')
        import os
        with open(zipName,"r") as fileHandle:
            form.add_file('file', os.path.split(zipName)[1], fileHandle)
        
        request = form.get_request()
        response = self.openRequest(request)
        print(response.read())

    def fetchRemoteProjects(self):
        response = self.getURLWithValues(self.baseURL+"?listProjects",{})
        responseText = response.read()
        import json
        remoteProjects = json.loads(responseText)
        return remoteProjects

    def fetchProjectWithName(self, projectName, jobList=None):
        remoteProjects = self.fetchRemoteProjects()
        matchingProjects = [project for project in remoteProjects if project[1] == projectName]
        if len(matchingProjects) == 1:
            self.fetchProjectWithId(matchingProjects[0][0], jobList=jobList)

    def fetchProjectWithId(self, projectId, jobList=None):
        response = self.getURLWithValues(self.baseURL+"/ProjectZipForProjectId/"+projectId,{"jobList":jobList})
        import tempfile
        tmpArchive = tempfile.NamedTemporaryFile(suffix='.ccp4_project.zip',delete=False)
        print(tmpArchive.name)
        CHUNK = 16 * 1024
        while True:
            chunk = response.read(CHUNK)
            if not chunk: break
            tmpArchive.write(chunk)
        tmpArchive.close()
        from core import CCP4NonGuiProjectUtils
        importer = CCP4NonGuiProjectUtils.CCP4NonGuiProjectUtils(tmpArchive.name)
        self.projectsManager().db().commit()
        return tmpArchive.name

    def lastRefmacExportForProjectId(self, projectId):
        response = self.getURLWithValues(self.baseURL+"/LastRefmacExportForProjectId/"+projectId,{})
        import tempfile
        tmpArchive = tempfile.NamedTemporaryFile(suffix='.zip',delete=False)
        print(tmpArchive.name)
        CHUNK = 16 * 1024
        while True:
            chunk = response.read(CHUNK)
            if not chunk: break
            tmpArchive.write(chunk)
        tmpArchive.close()
        return tmpArchive.name

    def lastLigandPipelineExportForProjectId(self, projectId):
        response = self.getURLWithValues(self.baseURL+"/LastLigandPipelineExportForProjectId/"+projectId,{})
        import tempfile
        tmpArchive = tempfile.NamedTemporaryFile(suffix='.zip',delete=False)
        print(tmpArchive.name)
        CHUNK = 16 * 1024
        while True:
            chunk = response.read(CHUNK)
            if not chunk: break
            tmpArchive.write(chunk)
        tmpArchive.close()
        return tmpArchive.name

    def syncProject(self, projectId):
        #Retrieve local job List
        localProjects = self.projectsManager().db().listProjects()
        # convert project Ids to unicode
        localProjectIds = [unicode(project[0],"utf-8") for project in localProjects]
        
        #Retrieve remote project list
        remoteProjects = self.fetchRemoteProjects()
        # Remote project Ids will be unicode already
        remoteProjectIds = [project[0] for project in remoteProjects]
        
        #Handle simple cases where entire project has to be pulled or pushed from archive
        #Convert projectId to unicode
        try:
            uProjectId = unicode(projectId,"utf-8")
        except:
            uProjectId = projectId
        
        if (uProjectId in localProjectIds) and (not uProjectId in remoteProjectIds):
            self.pushProjectWithId(uProjectId)
            return
        elif (not uProjectId in localProjectIds) and  (uProjectId in remoteProjectIds):
            self.fetchProjectWithId(uProjectId)
            return

        #Handle more complicated case where ony want to push/pull selected jobs
        localProjectJobList = self.projectsManager().db().getProjectJobListInfo(projectId=projectId, topLevelOnly=False)
        localProjectJobIds = [unicode(jobInfo['jobid'],"utf-8") for jobInfo in localProjectJobList]
        #print "Local job Ids", localProjectJobIds
        
        #Retrieve remote job List
        response = self.getURLWithValues(self.baseURL+"?getProjectJobListInfo",{"projectId":str(projectId),"topLevelOnly":"False"})
        responseText = response.read()
        import json
        remoteProjectJobList = json.loads(responseText)
        remoteProjectJobIds = [jobInfo['jobid'] for jobInfo in remoteProjectJobList]
        
        jobIdsMissingFromRemote = [jobId for jobId in localProjectJobIds if jobId not in remoteProjectJobIds]
        #print "jobIdsMissingFromRemote", jobIdsMissingFromRemote
        jobIdsMissingFromLocal = [jobId for jobId in remoteProjectJobIds if jobId not in localProjectJobIds]
        #print "jobIdsMissingFromLocal", jobIdsMissingFromLocal
        
        #print "jobIdsMissingFromRemote",jobIdsMissingFromRemote
        #print "jobIdsMissingFromLocal",jobIdsMissingFromLocal
        if len(jobIdsMissingFromRemote)>0: self.pushProjectWithId(projectId, jobList=jobIdsMissingFromRemote)
        if len(jobIdsMissingFromLocal)>0: self.fetchProjectWithId(projectId, jobList=jobIdsMissingFromLocal)


