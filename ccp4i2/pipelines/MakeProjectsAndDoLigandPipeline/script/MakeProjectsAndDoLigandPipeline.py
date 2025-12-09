"""
    MakeProjectsAndDoLigandPipeline.py: CCP4 GUI Project
    
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
    """

from __future__ import print_function

import os
from ccp4i2.core.CCP4PluginScript import CPluginScript
from ccp4i2.core import CCP4ErrorHandling
from ccp4i2.core import CCP4Utils

class MakeProjectsAndDoLigandPipeline(CPluginScript):
    TASKNAME = 'MakeProjectsAndDoLigandPipeline'   # Task name - should be same as class name and match pluginTitle in the .def.xml file
    TASKVERSION= 0.1               # Version of this plugin
    MAINTAINER = 'martin.noble@ncl.ac.uk'
    ERROR_CODES = {201 : {'severity':CCP4ErrorHandling.SEVERITY_WARNING, 'description' : 'Failed to create project' },
        202 : {'severity':CCP4ErrorHandling.SEVERITY_WARNING, 'description' : 'Failed to create job' },
        203 : {'severity':CCP4ErrorHandling.SEVERITY_WARNING, 'description' : 'Failed to determine job number' },
        204 : {'severity':CCP4ErrorHandling.SEVERITY_WARNING, 'description' : 'Failed to generate dbHandler' },
        205 : {'severity':CCP4ErrorHandling.SEVERITY_WARNING, 'description' : 'Failed to create job directory' },
        206 : {'severity':CCP4ErrorHandling.SEVERITY_WARNING, 'description' : 'Failed to retrieve plugin' },
        207 : {'severity':CCP4ErrorHandling.SEVERITY_WARNING, 'description' : 'Failed to instantiate plugin' },
        208 : {'severity':CCP4ErrorHandling.SEVERITY_WARNING, 'description' : 'Failed to attach db data to plugin' },
        209 : {'severity':CCP4ErrorHandling.SEVERITY_WARNING, 'description' : 'Failed to launch subprocess' },
        210 : {'severity':CCP4ErrorHandling.SEVERITY_WARNING, 'description' : 'Failed to reparent project' },
        211 : {'severity':CCP4ErrorHandling.SEVERITY_WARNING, 'description' : 'Failed to open and /or parse XML' },
        212 : {'severity':CCP4ErrorHandling.SEVERITY_ERROR, 'description' : 'Failed to determine projectId or Directory of presumed existing project' },
                    }
    PURGESEARCHLIST = [ [ 'hklin.mtz' , 0 ],
                        ['log_mtzjoin.txt', 0]
                       ]
    RUNEXTERNALPROCESS = False
    ASYNCHRONOUS=False

    #Uncomment the following if this pipeline will run plugins asynchronously
    #ASYNCHRONOUS=True

    def __init__(self, *args, **kws):
        super(MakeProjectsAndDoLigandPipeline, self).__init__(*args, **kws)
        from lxml import etree
        self.xmlroot = etree.Element("MakeProjectsAndDoLigandPipeline")
        etree.SubElement(self.xmlroot,"Message").text = "Starting processing"
        etree.SubElement(self.xmlroot,"Warnings").text = "No warnings"
        self.jobIds = []
        self.runningJobs = []
        self.processes = {}
        self.pluginsStarted = {}
        self.datasetElements = {}
        self.dumpXml()

    #The startProcess method is where you build in the pipeline logic
    def startProcess(self, command, **kws):
        from ccp4i2.core.CCP4Modules import PROJECTSMANAGER, JOBCONTROLLER
        from dbapi import CCP4DbApi
        import sys, os
        import functools
        from datetime import datetime
        from lxml import etree
        from ccp4i2.baselayer import QtCore
        pm = PROJECTSMANAGER()
        
        for iLigand, projectName in enumerate(self.container.inputData.PROJECTNAME_LIST):
            print(iLigand, projectName)
            messageNode = self.xmlroot.xpath("Message")[0]
            messageNode.text = "Starting on project " + projectName.__str__() + " which is " + str(iLigand+1) + "/" + str(len(self.container.inputData.PROJECTNAME_LIST))
            pluginName = "SubstituteLigand"
            pluginTitle = "SubstituteLigand"
            projectNameStr = projectName.__str__()
            
            datasetElement = etree.SubElement(self.xmlroot,"Dataset")
            etree.SubElement(datasetElement,"ProjectName").text = projectNameStr
            etree.SubElement(datasetElement,"SMILES").text = self.container.inputData.SMILES_LIST[iLigand].__str__()
            etree.SubElement(datasetElement,"MaximumResolution").text = "-"
            etree.SubElement(datasetElement,"Rfree").text = "-"
            etree.SubElement(datasetElement,"Rfactor").text = "-"
            self.dumpXml()

            #Scriptmatically create a new project with specified name and location
            projectId = None
            projectPath = os.path.join(CCP4Utils.getProjectDirectory(),projectNameStr)
            try:
                projectId = pm.createProject(projectName=projectNameStr, projectPath=projectPath)
            except:
                self.appendErrorReport(201, projectNameStr)
                self.xmlroot.xpath("//Warnings")[0].text = self._errorReport.report(user=False, ifStack=False)
                self.dumpXml()
                try:
                    projectId = pm.db().getProjectId(projectNameStr)
                    projectPath = pm.db().getProjectDirectory(projectId=projectId)
                except:
                    self.appendErrorReport(212, projectNameStr)
                    self.xmlroot.xpath("//Warnings")[0].text = self._errorReport.report(user=False, ifStack=False)
                    return CPluginScript.FAILED

            #Scriptmatically reparent project to be child of this one
            try:
                pm.db().updateProject(projectId,'parentprojectid',self.projectId())
                print("Reparented", projectId)
            except:
                self.appendErrorReport(210, projectNameStr+" "+projectId+" to "+str(self.projectId()))
                self.xmlroot.xpath("//Warnings")[0].text = self._errorReport.report(user=False, ifStack=False)
                continue

            #Create a placeholder job
            try:
                jobId = pm.db().createJob(projectId, pluginName, jobTitle=pluginTitle)
                print("jobId is", jobId)
            except CException as err:
                self.appendErrorReport(202, projectId+" "+pluginName)
                self.xmlroot.xpath("//Warnings")[0].text = self._errorReport.report(user=False, ifStack=False)
                continue
            
            #Record the mapping between jobId and etree Element
            self.jobIds.append(jobId)
            self.datasetElements[jobId] = datasetElement

            #Determine the "number" of the job within the project
            jobNumber = None
            try:
                jobNumber = pm.db().getJobInfo(jobId,mode='jobnumber')
                print("jobNumber is", jobNumber)
            except:
                self.appendErrorReport(203, str(jobId))
                self.xmlroot.xpath("//Warnings")[0].text = self._errorReport.report(user=False, ifStack=False)
                continue

            #Create a dbHandler for the plugin
            from ccp4i2.core.CCP4PluginScript import CDatabaseHandler
            try:
                dbHandler = CDatabaseHandler(projectId= projectId,jobNumber = jobNumber, projectName=projectNameStr)
                # Db is already running so _dbHandler should pick up that one
                dbOk = dbHandler.openDb()
                print("dbHandler opened")
            except CException as e:
                self.appendErrorReport(204, "")
                self.xmlroot.xpath("//Warnings")[0].text = self._errorReport.report(user=False, ifStack=False)
                continue

            #Create job's working directory
            workDir=None
            try:
                workDir = os.path.join(projectPath,'CCP4_JOBS','job_'+str(jobNumber))
                if not os.path.exists(workDir): os.mkdir(workDir)
                print("Able to create", workDir)
            except:
                print("Unable to create workDir ",workDirectory)
                self.appendErrorReport(205, "")
                self.xmlroot.xpath("//Warnings")[0].text = self._errorReport.report(user=False, ifStack=False)
                continue
            name = projectNameStr+'_'+str(jobNumber)

            #Retrieve plugin
            from ccp4i2.core import CCP4TaskManager
            cls = CCP4TaskManager.TASKMANAGER().getPluginScriptClass(pluginName)
            if cls is None:
                self.appendErrorReport(206, pluginName)
                self.xmlroot.xpath("//Warnings")[0].text = self._errorReport.report(user=False, ifStack=False)
                continue

            #Instantiate plugin
            try:
                plugin = cls(parent=None,name=name,workDirectory=workDir)
                print("able to instantiate plugin")
            except:
                self.appendErrorReport(207, pluginName)
                self.xmlroot.xpath("//Warnings")[0].text = self._errorReport.report(user=False, ifStack=False)
                continue

            if pluginTitle is not None and plugin.container is not None: plugin.container.header.pluginTitle = pluginTitle

            #And link the plugin back to the dbHandler
            try:
                plugin.setDbData(handler=dbHandler,projectName=projectNameStr,projectId=projectId,jobNumber=jobNumber,jobId=jobId)
                print("Able to link to dbHandler")
            except CException as e:
                self.appendErrorReport(208, pluginName + e.__str__())
                self.xmlroot.xpath("//Warnings")[0].text = self._errorReport.report(user=False, ifStack=False)
                continue

            #Configure the pipeline inputs and controlParameters
            plugin.container.controlParameters.LIGANDAS="SMILES"
            plugin.container.inputData.SMILESIN=self.container.inputData.SMILES_LIST[iLigand]
            plugin.container.inputData.XYZIN=self.container.inputData.XYZIN
            if self.container.inputData.FREERFLAG.isSet():
                plugin.container.inputData.FREERFLAG_IN=self.container.inputData.FREERFLAG
            plugin.container.inputData.PIPELINE=self.container.inputData.PIPELINE
            
            unmergedInput = plugin.container.inputData.UNMERGEDFILES[-1]
            unmergedPath = os.path.join(self.container.inputData.ROOT_DIRECTORY.fullPath.__str__(),
                                        self.container.inputData.PATH_LIST[iLigand].__str__()[1:])
            unmergedInput.file.fullPath = unmergedPath
            unmergedInput.crystalName = "Xtal1"
            unmergedInput.dataset = "DS1"
            print('fullUnmergedPath',unmergedInput.file.fullPath.__str__())

            plugin.saveParams()
            paramsFileName = plugin.makeFileName('PARAMS')
            inputParamsFileName = plugin.makeFileName('JOB_INPUT')
            import shutil
            shutil.copyfile(paramsFileName, inputParamsFileName)
            print('From to',paramsFileName, inputParamsFileName)

            self.processes[jobId] = JOBCONTROLLER().runTask(jobId=jobId)
            print('JOBCONTROLLER', JOBCONTROLLER())

            self.pluginsStarted[jobId] = plugin
            self.runningJobs.append(jobId)
            #self.connectSignal(self.pluginsStarted[jobId],'finished',
            #functools.partial(self.pluginFinished, jobId))
            while len(self.runningJobs) > 10:
                print(str(datetime.now())+"::Scanning for finished ")
                sys.stdout.flush()
                for jobId in self.runningJobs:
                    if self.processes[jobId].waitForFinished(1000):
                        self.pluginFinished(jobId)
                        break
        
        while len(self.runningJobs) > 0:
            print(str(datetime.now())+"::Scanning for finished ")
            sys.stdout.flush()
            for jobId in self.runningJobs:
                if self.processes[jobId].waitForFinished(1000):
                    self.pluginFinished(jobId)
                    break

        return CPluginScript.SUCCEEDED

    #This method will be called as each plugin completes
    def pluginFinished(self, jobId):
        print("Have heard from ", jobId)
        from lxml import etree
        whichPlugin = self.pluginsStarted[jobId]
        self.runningJobs.remove(jobId)
        datasetElement = self.datasetElements[jobId]
        datasetElement.xpath("MaximumResolution")[0].text = "N/D"
        datasetElement.xpath("Rfree")[0].text = "N/D"
        datasetElement.xpath("Rfactor")[0].text = "N/D"
        try:
            print("name of correspnding PROGRAMXML", whichPlugin.makeFileName("PROGRAMXML"))
            pluginEtree = CCP4Utils.openFileToEtree(whichPlugin.makeFileName("PROGRAMXML"))
            try:
                rnodeText = pluginEtree.xpath('//ResolutionLimitEstimate[@type="CChalf" and ./Direction[text()="Overall"]]/MaximumResolution/text()')[-1]
                datasetElement.xpath("MaximumResolution")[0].text = rnodeText
            except:
                pass
            
            try:
                rfreeNodeText = pluginEtree.xpath('//Overall_stats/stats_vs_cycle/new_cycle[last()]/r_free/text()')[-1]
                datasetElement.xpath("Rfree")[0].text = rfreeNodeText
            except:
                pass
                    
            try:
                rfactorNodeText = pluginEtree.xpath('//Overall_stats/stats_vs_cycle/new_cycle[last()]/r_factor/text()')[-1]
                datasetElement.xpath("Rfactor")[0].text = rfreeNodeText
            except:
                pass
        except Exception as e:
            self.appendErrorReport(211, jobId)
            self.xmlroot.xpath("//Warnings")[0].text = self._errorReport.report(user=False, ifStack=False)

        self.dumpXml()
        
        return CPluginScript.SUCCEEDED

    def dumpXml(self):
        from lxml import etree
        with open(self.makeFileName("PROGRAMXML"),"w") as programXmlFile:
            CCP4Utils.writeXML(programXmlFile,etree.tostring(self.xmlroot))
