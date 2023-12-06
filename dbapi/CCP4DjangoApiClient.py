import unittest
import requests
from pathlib import Path
from urllib import parse
import collections.abc

CCP4i2UUID = str
CCP4i2ParamsPath = str


class CCP4DjangoApiClient:

    """A python class providing utility methods to access the CCP4i2 server API"""

    def __init__(self, urlBase: str = "http://127.0.0.1:43434/database"):
        """
        Args:
            urlBase: the root url for the connection to the CCP4i2 server, typically `http://127.0.0.1:43434/database`
        """
        #: urlBase stores the root url for the connection to the CCP4i2 server
        self.urlBase = urlBase

    def setJobParameterByXML(self, jobId: CCP4i2UUID, objectPath: CCP4i2ParamsPath, valueXmlText: str):
        """
        Update the params.xml of the CCP4i2 task by specifying xml describing the new parameter
        Args:

            jobId: the primary key of the CCP4i2 job to be patched
            objectPath: the XML path to the parameter (e.g. `controlParameters.NCYCLES`)
            valueXmlText: a string representation of the XML for the updated element (e.g. `<NCYCLES>10</NCYCLES>`)
        """
        files = []
        data = {"jobId": jobId,
                "objectPath": objectPath,
                "valueXMLText": valueXmlText}
        url = f"{self.urlBase}/setJobParameterByXML"
        response = requests.post(url=url,
                                 data=data,
                                 files=files)
        jsonData = response.json()
        return jsonData

    def setJobParameter(self, jobId: CCP4i2UUID, objectPath: CCP4i2ParamsPath, value):
        """
        Update the params.xml of the CCP4i2 task by specifying the value of a "simple" control element
        Args:

            jobId: the primary key of the CCP4i2 job to be patched
            objectPath: the XML path to the parameter (e.g. `controlParameters.NCYCLES`)
            value: the new value to be applied (e.g. 10)
        """
        tag = objectPath.split(".")[-1]
        valueXmlText = f"<{tag}>{value}</{tag}>"
        return self.setJobParameterByXML(jobId, objectPath, valueXmlText)

    def uploadFileForJobObject(self, filePath: str, objectPath: CCP4i2ParamsPath, jobId: CCP4i2UUID = None, projectId: CCP4i2UUID = None,
                               projectName: str = None, jobNumber: str = None):
        """
        Upload an input file to a specified CCP4i2 task
        Args:
            filePath: path the the file to be uploaded
            objectPath: the XML path specifying which input the file corresponds to (e.g. `inputData.XYZIN`)

            Job specifier:

                Either:

                    jobId: jobid of the job for which the file is input

                Or:

                    Either: 

                        projectId: projectid of the project containg in the job for which the file is input
                    Or:

                        projectName: projectname of the project containg in the job for which the file is input (e.g. `MDM2`)

                    And:
                        jobNumber: number of the job for which the file is input (e.g. `10`)
        """
        filePathAsPath = Path(filePath)
        url = f"{self.urlBase}/uploadFileForJobObject"
        files = [("file", open(filePathAsPath, "rb"))]
        data = {"objectPath": objectPath,
                "fileName": filePathAsPath.name,
                }
        if jobId is not None:
            data['jobId'] = jobId
        if jobNumber is not None:
            data['jobNumber'] = jobNumber
        if projectId is not None:
            data['projectId'] = projectId
        if projectName is not None:
            data['projectName'] = projectName
        response = requests.post(url=url, data=data, files=files)
        uploadResult = response.json()
        return uploadResult

    def createJob(self, projectId: CCP4i2UUID = None, projectName: str = None, taskName: str = None):
        """
        Create a new CCP4i2 job
        Args:

            taskName: name of the CCP4i2 plugin to be created (e.g. `prosmart_refmac`)

            Project specifier:

                Either: 

                    projectId: projectid of the project in which to createjob

                Or:

                    projectName: projectname of the project in which to createjob (e.g. `MDM2`)

        """
        url = f"{self.urlBase}/createJob"
        files = []
        data = {"taskName": taskName}
        if projectId is not None:
            data["projectId"] = projectId
        if projectName is not None:
            data["projectName"] = projectName

        response = requests.post(url=url, data=data, files=files)
        jsonData = response.json()
        return jsonData

    def modelValues(self, type: str, values: list = None, order_by: list = None, **predicate) -> dict:
        """
        Method to flexibly interrogate the CCP4i2 database
        Args:

            type: Name of the CCP4i2 database object type about which this query applies (e.g. `Projects`, `Jobs`, or `Files`)
            values: List of the fields to be fetched (e.g. `["projectname", "projectdirectory", "projectid"]`).  

                Note that fields are not limited to the table of the object itself. For example, each `Jobs` object has a 
                ForeignKey (stored in the `projectid` field) to a corresponding `Projects` object. Using this, one can 
                retrieve in one call properties of the Job itself, and of the Job's Project using a "double underscore" syntax 
                (e.g. to retrieve a job's jobid, number, task name and the name of the job's project, specify
                  `["jobid", "jobnumber", "taskname", "projectid__projectname"]`).

                If no value is provided for `values`, then the default is to retrieve all of the fields in the object's table

            **predicate: List of filters to be applied in selecting objects from the database. Examples of calls exploiting filters:

                `modelValues("Jobs", taskname="prosmart_refmac")`: 
                Returns the default fields of all prosmart_refmac jobs in the database. Note 
                that all modelfields used in filters or values are in *lower case*

                `modelValues("Jobs", taskname="prosmart_refmac", parentjobid__isnull=True)`: 
                Returns the default fields of all "top-level" prosmart_refmac jobs in the database. Note 
                use of the "__isnull" modifier to test whether the field value is null.

                `modelValues("Jobs", taskname="prosmart_refmac", status__statusname="Finished")`: 
                Returns the default fields of all successfully completed prosmart_refmac jobs in the database

                `modelValues("Jobs", taskname="prosmart_refmac", status__statusname__in=["Finished", "Failed", "Unsatisfactory"])`: 
                Returns the default fields of all prosmart_refmac jobs that have stopped because completed, failed, or terminated with status 
                "unsatisfactory".  Note the use of the "__in" modifier for specifying multiple options

                For a general introduction to django's filtering predicates, see https://docs.djangoproject.com/en/4.2/topics/db/queries/#retrieving-specific-objects-with-filters
                
        Returns:

            A dict with entries:

                status - should have value "Success"
                results - A list of dictionaries containing entries for each of the requested values for each object matching the predicate
                exception - A textual summary of an exception if one was raised in executing the call

        """
        url = f"{self.urlBase}/ModelValues"
        params = {"__type__": type}
        if values is not None:
            params["__values__[]"] = values
        if order_by is not None:
            params["__order_by__[]"] = values
        for key in predicate:
            if isinstance(predicate[key], collections.abc.Sequence) and not isinstance(predicate[key], str):
                params[f"{key}[]"] = predicate[key]
            else:
                params[key] = predicate[key]
        url = f"{url}?{parse.urlencode(params, doseq=True)}"
        response = requests.get(url=url)
        jsonData = response.json()
        return jsonData

    def getProjectJobFile(self, projectId: CCP4i2UUID = None, projectName: str = None, fileName: str = None, jobNumber: str = "1")  -> requests.Response:
        """
        fetch a unicode-encoded file with name fileName from the directory of the specified job
        Args:

            Either:

                projectId: projectid of the project of the job in question

            Or:

                projectName: name of the project of the job in question

            fileName: path of the file relative to the job's working directory
            jobNumber: period delimited job "path" of the job in question (e.g. subjob 1 of subjob 1 of job 32 would have `jobNumber="32.1.1"` )

        Returns:
            Requests.response (see https://requests.readthedocs.io/en/latest/user/quickstart/#response-status-codes )         
        """
        url = f"{self.urlBase}/getProjectJobFile"
        params = {"fileName": fileName, "jobNumber": jobNumber}
        if projectId is not None: params["projectId"] = projectId
        if projectName is not None: params["projectName"] = projectName
        url = f"{url}?{parse.urlencode(params, doseq=True)}"
        response = requests.get(url=url)
        return response

    def getFileWithPredicate(self, **predicate) -> requests.Response:
        """
        fetch a file known to the CCP4 database, specified uniquely by predicate
        Args:

            predicate: any set of filters that will identify a single file.

                Examples: 
                
                    getFileWithPredicate(fileid=fileId): retrieve a file based on fileid 

                    getFileWithPredicate(
                        jobid__projectid__projectname="PhasingTest", 
                        jobid__jobnumber="104",
                        jobparamname="XYZOUT"): retrieve the file generated as XYZOUT on job number 104 of project with name PhasingTest

        Returns:
            Requests.response (see https://requests.readthedocs.io/en/latest/user/quickstart/#response-status-codes )                 
        """
        url = f"{self.urlBase}/getFileWithPredicate"
        url = f"{url}?{parse.urlencode(predicate, doseq=True)}"
        response = requests.get(url=url)
        print (type(response))
        return response
    
class MyTest(unittest.TestCase):

    def test_method1(self):
        client = CCP4DjangoApiClient("http://127.0.0.1:43434/database")

    def test_retrieve_status_name(self):
        client = CCP4DjangoApiClient("http://127.0.0.1:43434/database")
        createJobResult = client.createJob(
            projectName="CDK2CyclinACDC25", taskName="MakeLink")
        self.assertEqual(createJobResult["status"], "Success")
        fetchStatusResult = client.modelValues(
            "Jobs", values=["status__statustext"], jobid=createJobResult["jobId"])
        self.assertEqual(
            fetchStatusResult["results"][0]["status__statustext"], "Pending")

    def test_retrieve_finished_jobs(self):
        client = CCP4DjangoApiClient("http://127.0.0.1:43434/database")
        finishedJobs = client.modelValues("Jobs", values=[
                                          "projectid__projectname", "jobnumber", "taskname"], status__statustext="Finished")
        self.assertEqual(finishedJobs, "banana")

    def test_createMakeLink(self):
        client = CCP4DjangoApiClient("http://127.0.0.1:43434/database")
        result = client.createJob(
            projectName="CDK2CyclinACDC25", taskName="MakeLink")
        jobId = result["jobId"]
        result = client.uploadFileForJobObject(
            "TestData/THR-PO4-prelink.pdb", "inputData.XYZIN", jobId=jobId)
        result = client.setJobParameter(
            jobId, "controlParameters.TOGGLE_LINK", "True")


if __name__ == "__main__":

    if False:
        client = CCP4DjangoApiClient("http://127.0.0.1:43434/database")
        createJobResult = client.createJob(
            projectName="CDK2CyclinACDC25", taskName="MakeLink")
        uploadFileResult = client.uploadFileForJobObject(
            "TestData/THR-PO4-prelink.pdb", "inputData.XYZIN", jobId=createJobResult["jobId"])
        setParameterResult = client.setJobParameter(
            createJobResult["jobId"], "controlParameters.TOGGLE_LINK", "True")

    if True:
        client = CCP4DjangoApiClient("http://127.0.0.1:43434/database")
        response = client.getFileWithPredicate(
                        jobid__projectid__projectname="PhasingTest", 
                        jobid__jobnumber="104",
                        jobparamname="XYZOUT")
        