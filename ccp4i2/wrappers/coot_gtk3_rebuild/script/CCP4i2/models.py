from __future__ import print_function    # (at top of module)

# This is an auto-generated Django model module.
# You'll have to do the following manually to clean this up:
#   * Rearrange models' order
#   * Make sure each model has one field with primary_key=True
#   * Make sure each ForeignKey has `on_delete` set to the desired behavior.
#   * Remove `managed = True` lines if you wish to allow Django to create, modify, and delete the table
# Feel free to rename the models, but don't rename db_table values or field names.
from __future__ import unicode_literals
from asyncio.log import logger

from django.db import models
import uuid
import os
from django.db.models import Max
from django.conf import settings
from django.db.models.signals import pre_delete
from django.dispatch import receiver
from django.utils.text import slugify
import pathlib
import shutil
import time


def generate_ccp4_uuid():
    return uuid.uuid1().hex

class Comments(models.Model):
    commentid = models.CharField(db_column='CommentID', unique=True, max_length=40,
                                 primary_key=True, default=generate_ccp4_uuid)  # Field name made lowercase.
    # Field name made lowercase.
    jobid = models.ForeignKey('Jobs', on_delete=models.CASCADE,
                              db_column='JobID', related_name='commentsofjob', null=True)
    userid = models.ForeignKey('Users', on_delete=models.SET_NULL, db_column='UserID',
                               null=True,  related_name='commentsofuser')  # Field name made lowercase.
    # Field name made lowercase. This field type is a guess.
    timeofcomment = models.TextField(db_column='TimeOfComment')
    # Field name made lowercase.
    comment = models.TextField(db_column='Comment', blank=True, null=True)

    class Meta:
        managed = True
        db_table = 'Comments'
        unique_together = (('jobid', 'commentid'),)
        app_label = 'CCP4i2'


class Databases(models.Model):
    databaseid = models.CharField(db_column='DatabaseId', unique=True, max_length=40,
                                  primary_key=True, default=generate_ccp4_uuid)  # Field name made lowercase.
    # Field name made lowercase. This field type is a guess.
    thisisme = models.SmallIntegerField(db_column='ThisIsMe')
    # Field name made lowercase.
    schemaversion = models.CharField(
        db_column='SchemaVersion', max_length=50, blank=True, null=True)
    # Field name made lowercase.
    schemadate = models.CharField(
        db_column='SchemaDate', max_length=50, blank=True, null=True)
    # Field name made lowercase.
    creationtime = models.FloatField(
        db_column='CreationTime', blank=True, null=True)
    # Field name made lowercase.
    creationhostname = models.CharField(
        db_column='CreationHostName', max_length=255, blank=True, null=True)
    # Field name made lowercase.
    creationusername = models.CharField(
        db_column='CreationUserName', max_length=100, blank=True, null=True)

    class Meta:
        managed = True
        db_table = 'Databases'
        app_label = 'CCP4i2'


class Directoryaliases(models.Model):
    directoryalias = models.CharField(db_column='DirectoryAlias', unique=True,  max_length=100,
                                      primary_key=True, default=generate_ccp4_uuid)  # Field name made lowercase.
    # Field name made lowercase.
    directory = models.TextField(db_column='Directory', blank=True, null=True)

    class Meta:
        managed = True
        db_table = 'DirectoryAliases'
        app_label = 'CCP4i2'


class Exportfiles(models.Model):
    exportid = models.CharField(db_column='ExportId', unique=True, max_length=40,
                                primary_key=True, default=generate_ccp4_uuid)  # Field name made lowercase.
    # Field name made lowercase.
    fileid = models.ForeignKey('Files', on_delete=models.CASCADE,
                               db_column='FileID', related_name='exportfiles', null=True)
    # Field name made lowercase.
    exportfilename = models.TextField(
        db_column='ExportFilename', blank=True, null=True)
    # Field name made lowercase.
    annotation = models.TextField(
        db_column='Annotation', blank=True, null=True)
    # Field name made lowercase.
    creationtime = models.FloatField(
        db_column='CreationTime', blank=True, null=True)

    class Meta:
        managed = True
        db_table = 'ExportFiles'
        app_label = 'CCP4i2'


def nextFileassociationmembers():
    if Fileassociationmembers.objects.count() == 0:
        return 1
    return Fileassociationmembers.objects.only('fileassociationmembersid').aggregate(Max('fileassociationmembersid'))['fileassociationmembersid__max']+1


# Here seeing if rowid can serve as pk
class Fileassociationmembers(models.Model):
    # Field name made lowercase.
    fileassociationmembersid = models.IntegerField(
        db_column='rowid', unique=True, primary_key=True, default=nextFileassociationmembers)
    fileassociationid = models.ForeignKey('Fileassociations', on_delete=models.CASCADE, db_column='FileAssociationID',
                                          related_name='associationgroups', null=True)  # Field name made lowercase.
    fileid = models.ForeignKey('Files', on_delete=models.CASCADE, db_column='FileID',
                               related_name='associationgroupsoffile', null=True)  # Field name made lowercase.
    # Field name made lowercase. This field type is a guess.
    fileassociationroleid = models.ForeignKey('Fileassociationroles', on_delete=models.CASCADE,
                                              db_column='FileAssociationRoleID', related_name='associationgroupsofrole', null=True)

    class Meta:
        managed = True
        db_table = 'FileAssociationMembers'
        unique_together = (
            ('fileid', 'fileassociationid', 'fileassociationroleid'),)
        app_label = 'CCP4i2'


class Fileassociationroles(models.Model):
    # Field name made lowercase. This field type is a guess.
    fileassociationroleid = models.IntegerField(
        db_column='FileAssociationRoleID', unique=True, primary_key=True)
    # Field name made lowercase.
    fileassociationtypeid = models.ForeignKey(
        'Fileassociationtypes', on_delete=models.CASCADE, db_column='FileAssociationTypeID', null=True)
    # Field name made lowercase.
    fileassociationrolename = models.CharField(
        db_column='FileAssociationRoleName', unique=True, max_length=50)
    fileassociationroledescription = models.TextField(
        db_column='FileAssociationRoleDescription', blank=True, null=True)  # Field name made lowercase.

    class Meta:
        managed = True
        db_table = 'FileAssociationRoles'
        ordering = ['fileassociationroleid']
        app_label = 'CCP4i2'


class Fileassociationtypes(models.Model):
    # Field name made lowercase. This field type is a guess.
    fileassociationtypeid = models.IntegerField(
        db_column='FileAssociationTypeID', unique=True, primary_key=True)
    # Field name made lowercase.
    fileassociationtypename = models.CharField(
        db_column='FileAssociationTypeName', unique=True, max_length=50)
    fileassociationtypedescription = models.TextField(
        db_column='FileAssociationTypeDescription', blank=True, null=True)  # Field name made lowercase.

    class Meta:
        managed = True
        db_table = 'FileAssociationTypes'
        ordering = ['fileassociationtypeid']
        app_label = 'CCP4i2'


class Fileassociations(models.Model):
    # Field name made lowercase.
    fileassociationid = models.CharField(
        db_column='FileAssociationID', unique=True, max_length=40, primary_key=True, default=generate_ccp4_uuid)
    # Field name made lowercase. This field type is a guess.
    fileassociationtypeid = models.ForeignKey(
        'Fileassociationtypes', on_delete=models.CASCADE, db_column='FileAssociationTypeId', null=True)
    # Field name made lowercase.
    fileassociationname = models.CharField(
        db_column='FileAssociationName', max_length=50, blank=True, null=True)

    class Meta:
        managed = True
        db_table = 'FileAssociations'
        app_label = 'CCP4i2'


class Fileroles(models.Model):
    # Field name made lowercase. This field type is a guess.
    roleid = models.IntegerField(
        db_column='RoleID', unique=True, primary_key=True)
    # Field name made lowercase.
    roletext = models.CharField(
        db_column='RoleText', unique=True, max_length=50)

    class Meta:
        managed = True
        db_table = 'FileRoles'
        ordering = ['roleid']
        app_label = 'CCP4i2'

    def __unicode__(self):
        return '<djangoDB::{} {} {}>'.format(self.__class__.__name__, self.roleid, self.roletext)


class Filetypes(models.Model):
    # Field name made lowercase. This field type is a guess.
    filetypeid = models.IntegerField(
        db_column='FileTypeID', unique=True, primary_key=True)
    # Field name made lowercase.
    filetypename = models.CharField(
        db_column='FileTypeName', unique=True, max_length=50)
    # Field name made lowercase.
    filetypedescription = models.TextField(
        db_column='FileTypeDescription', blank=True, null=True)

    class Meta:
        managed = True
        db_table = 'FileTypes'
        ordering = ['filetypeid']
        app_label = 'CCP4i2'

    def __unicode__(self):
        return '<djangoDB::{} {} {} "{}">'.format(self.__class__.__name__, self.filetypeid, self.filetypename, self.filetypedescription)


def nextFileuses():
    if Fileuses.objects.count() == 0:
        return 1
    return Fileuses.objects.only('fileuseid').aggregate(Max('fileuseid'))['fileuseid__max']+1


class Fileuses(models.Model):
    # Field name made lowercase.
    fileuseid = models.IntegerField(
        db_column='rowid', unique=True, primary_key=True, default=nextFileuses)
    # Field name made lowercase.
    fileid = models.ForeignKey('Files', on_delete=models.CASCADE,
                               db_column='FileID', related_name='usesoffile', null=True)
    jobid = models.ForeignKey('Jobs', on_delete=models.CASCADE, db_column='JobID',
                              related_name='fileusesofjob', null=True, blank=True)  # Field name made lowercase.
    roleid = models.ForeignKey(Fileroles, on_delete=models.CASCADE, db_column='RoleID',
                               related_name='fileusesofrole', null=True)  # Field name made lowercase.
    # Field name made lowercase.
    jobparamname = models.CharField(
        db_column='JobParamName', max_length=40, blank=True, null=True)

    class Meta:
        managed = True
        db_table = 'FileUses'
        unique_together = (('fileid', 'jobid', 'roleid', 'jobparamname'),)
        ordering = ['fileuseid']
        app_label = 'CCP4i2'

    def __unicode__(self):
        return '<djangoDB::{} {} {} "{}">'.format(self.__class__.__name__, self.jobid.projectid.projectname, self.jobid.jobnumber, self.jobparamname)

    def __str__(self):
        return '<djangoDB::{} {} {} "{}">'.format(self.__class__.__name__, self.jobid.projectid.projectname, self.jobid.jobnumber, self.jobparamname)


class Files(models.Model):
    fileid = models.CharField(db_column='FileID', unique=True, primary_key=True,
                              max_length=40, default=generate_ccp4_uuid)  # Field name made lowercase.
    # Field name made lowercase.
    filename = models.TextField(db_column='Filename')
    # Field name made lowercase.
    annotation = models.TextField(
        db_column='Annotation', blank=True, null=True)
    # Field name made lowercase.
    filetypeid = models.ForeignKey('Filetypes', on_delete=models.CASCADE,
                                   db_column='FiletypeID', null=True, related_name='filesoftype')
    # Field name made lowercase. This field type is a guess.
    filesubtype = models.IntegerField(
        db_column='FileSubType', blank=True, null=True)
    # Field name made lowercase. This field type is a guess.
    filecontent = models.IntegerField(
        db_column='FileContent', blank=True, null=True)
    jobid = models.ForeignKey('Jobs', on_delete=models.CASCADE, db_column='JobID', max_length=40,
                              related_name='filesofjob', null=True, blank=True)  # Field name made lowercase.
    # Field name made lowercase.
    jobparamname = models.CharField(
        db_column='JobParamName', max_length=40, blank=True, null=True)
    # Field name made lowercase. This field type is a guess.
    pathflag = models.IntegerField(db_column='PathFlag', blank=True, null=True)

    class Meta:
        managed = True
        db_table = 'Files'
        app_label = 'CCP4i2'

    @property
    def relPath(self):
        import os
        # print 'In filePath'
        try:
            if self.pathflag == 1:
                jobSubPath = os.path.sep.join(
                    "job_"+digit for digit in self.jobid.jobnumber.split("."))
                filesDir = os.path.join("CCP4_JOBS", jobSubPath)
            elif self.pathflag == 2 or self.importid is not None:
                filesDir = os.path.join("CCP4_IMPORTED_FILES")
            else:
                print(
                    'Unexpected pathflag: assuming this is an imported file', self.pathflag)
                filesDir = os.path.join("CCP4_IMPORTED_FILES")
        except Exception as err:
            print('Exception', err)
        return filesDir

    @property
    def filePath(self):
        filePath = os.path.abspath(os.path.join(
            self.jobid.projectid.projectdirectory, self.relPath, self.filename))
        return filePath
'''
    @property
    def digest(self):
        from CCP4i2.shim import CCP4i2SharedRoutines
        from dbapi import CCP4DbApi
        from core import CCP4DataManager
        fileClassName = f'C{CCP4DbApi.FILETYPES_CLASS[int(self.filetypeid.filetypeid)]}'
        fileClass = CCP4DataManager.DATAMANAGER().getClass(fileClassName)
        instance = fileClass()
        instance.setFullPath(self.filePath)
        return CCP4i2SharedRoutines.fileAnalysis(instance)
'''

class Importfiles(models.Model):
    # Field name made lowercase.
    importid = models.CharField(db_column='ImportId', unique=True,
                                primary_key=True, max_length=40, default=generate_ccp4_uuid)
    # Field name made lowercase.
    fileid = models.OneToOneField(
        Files, on_delete=models.CASCADE, db_column='FileID', related_name='importid', null=True)
    # Field name made lowercase.
    sourcefilename = models.TextField(
        db_column='SourceFilename', blank=True, null=True)
    sourcefileid = models.ForeignKey(Files, on_delete=models.CASCADE, db_column='SourceFileID',
                                     blank=True, related_name='importsourcefileids', null=True)  # Field name made lowercase.
    # Field name made lowercase.
    exportfileid = models.ForeignKey(
        Exportfiles, on_delete=models.CASCADE, db_column='ExportFileID', blank=True, null=True)
    # Field name made lowercase.
    annotation = models.TextField(
        db_column='Annotation', blank=True, null=True)
    # Field name made lowercase.
    creationtime = models.FloatField(
        db_column='CreationTime', blank=True, null=True)
    # Field name made lowercase.
    lastmodifiedtime = models.FloatField(
        db_column='LastModifiedTime', blank=True, null=True)
    # Field name made lowercase.
    checksum = models.CharField(
        db_column='Checksum', max_length=40, blank=True, null=True)
    # Field name made lowercase. This field type is a guess.
    importnumber = models.TextField(
        db_column='ImportNumber', blank=True, null=True)
    # Field name made lowercase.
    reference = models.TextField(db_column='Reference', blank=True, null=True)

    class Meta:
        managed = True
        db_table = 'ImportFiles'
        app_label = 'CCP4i2'


class Jobevaluation(models.Model):
    # Field name made lowercase. This field type is a guess.
    evaluationid = models.IntegerField(
        db_column='EvaluationID', unique=True, primary_key=True)
    # Field name made lowercase.
    evaluationtext = models.CharField(
        db_column='EvaluationText', max_length=255)

    class Meta:
        managed = True
        db_table = 'JobEvaluation'
        ordering = ['evaluationid']
        app_label = 'CCP4i2'

    def __unicode__(self):
        return '<djangoDB::{} {} "{}">'.format(self.__class__.__name__, self.evaluationid, self.evaluationtext)


class Keytypes(models.Model):
    # Field name made lowercase. This field type is a guess.
    keytypeid = models.IntegerField(
        db_column='KeyTypeID', unique=True, primary_key=True)
    # Field name made lowercase.
    keytypename = models.CharField(
        db_column='KeyTypeName', unique=True, max_length=50)
    # Field name made lowercase.
    keytypedescription = models.TextField(
        db_column='KeyTypeDescription', blank=True, null=True)

    class Meta:
        managed = True
        db_table = 'KeyTypes'
        ordering = ['keytypeid']
        app_label = 'CCP4i2'

    def __unicode__(self):
        return '<djangoDB::{} {} {} "{}">'.format(self.__class__.__name__, self.keytypeid, self.keytypename,  self.keytypedescription)


def nextJobkeycharvalues():
    if Jobkeycharvalues.objects.count() == 0:
        return 1
    return Jobkeycharvalues.objects.only('jobkeycharvaluesid').aggregate(Max('jobkeycharvaluesid'))['jobkeycharvaluesid__max']+1


class Jobkeycharvalues(models.Model):  # Here seeing if rowid can serve as pk
    # Field name made lowercase.
    jobkeycharvaluesid = models.IntegerField(
        db_column='rowid', unique=True, primary_key=True, default=nextJobkeycharvalues)
    jobid = models.ForeignKey('Jobs', on_delete=models.CASCADE, db_column='JobID',
                              related_name='keycharvaluesofjob', null=True, blank=True)  # Field name made lowercase.
    # Field name made lowercase. This field type is a guess.
    keytypeid = models.ForeignKey('Keytypes', on_delete=models.CASCADE,
                                  db_column='KeyTypeID', related_name='keycharvaluesoftype', null=True)
    # Field name made lowercase.
    value = models.CharField(
        db_column='Value', max_length=255, blank=True, null=True)

    class Meta:
        managed = True
        db_table = 'JobKeyCharValues'
        unique_together = (('jobid', 'keytypeid'),)
        ordering = ['jobkeycharvaluesid']
        app_label = 'CCP4i2'


def nextJobkeyvalues():
    if Jobkeyvalues.objects.count() == 0:
        return 1
    return Jobkeyvalues.objects.only('jobkeyvaluesid').aggregate(Max('jobkeyvaluesid'))['jobkeyvaluesid__max']+1


class Jobkeyvalues(models.Model):  # Here seeing if rowid can serve as pk
    # Field name made lowercase.
    jobkeyvaluesid = models.IntegerField(
        db_column='rowid', unique=True, primary_key=True, default=nextJobkeyvalues)
    jobid = models.ForeignKey('Jobs', on_delete=models.CASCADE, db_column='JobID',
                              related_name='keyvaluesofjob', null=True, blank=True)  # Field name made lowercase.
    # Field name made lowercase. This field type is a guess.
    keytypeid = models.ForeignKey(Keytypes, on_delete=models.CASCADE,
                                  db_column='KeyTypeID', related_name='keyvaluesoftype', null=True)
    # Field name made lowercase.
    value = models.FloatField(db_column='Value', blank=True, null=True)

    class Meta:
        managed = True
        db_table = 'JobKeyValues'
        unique_together = (('jobid', 'keytypeid'),)
        ordering = ['jobkeyvaluesid']
        app_label = 'CCP4i2'


class Jobstatus(models.Model):
    # Field name made lowercase. This field type is a guess.
    statusid = models.IntegerField(
        db_column='StatusID', unique=True, primary_key=True)
    # Field name made lowercase.
    statustext = models.CharField(
        db_column='StatusText', unique=True, max_length=255)

    class Meta:
        managed = True
        db_table = 'JobStatus'
        ordering = ['statusid']
        app_label = 'CCP4i2'

    def __unicode__(self):
        return u'<djangoDB::{} {} {}>'.format(self.__class__.__name__, self.statusid, self.statustext)

    def __str__(self):
        return '<djangoDB::{} {} {}>'.format(self.__class__.__name__, self.statusid, self.statustext)


class Jobs(models.Model):
    jobid = models.CharField(db_column='JobID', unique=True, primary_key=True,
                             default=generate_ccp4_uuid, max_length=40)  # Field name made lowercase.
    # Field name made lowercase.
    jobnumber = models.CharField(db_column='JobNumber', max_length=50)
    # Field name made lowercase.
    creationtime = models.FloatField(
        db_column='CreationTime', blank=True, null=True, default=time.time)
    # Field name made lowercase.
    finishtime = models.FloatField(
        db_column='FinishTime', blank=True, null=True)
    # Field name made lowercase. This field type is a guess.
    status = models.ForeignKey(Jobstatus, on_delete=models.SET_NULL,
                               db_column='Status', related_name='jobsofstatus', null=True)
    # Field name made lowercase. This field type is a guess.
    evaluation = models.ForeignKey(Jobevaluation, on_delete=models.SET_NULL,
                                   db_column='Evaluation', related_name='jobsofevaluation', null=True, blank=True)
    # Field name made lowercase. This field type is a guess.
    useragent = models.ForeignKey('Useragents', on_delete=models.SET_NULL,
                                  db_column='UserAgent', related_name='jobsofuseragent', null=True)
    # Field name made lowercase.
    jobtitle = models.CharField(
        db_column='JobTitle', max_length=255, blank=True, null=True)
    projectid = models.ForeignKey('Projects', on_delete=models.CASCADE, db_column='ProjectID',
                                  related_name='jobsofproject', null=True)  # Field name made lowercase.
    # Field name made lowercase.
    taskname = models.CharField(db_column='TaskName', max_length=100)
    # Field name made lowercase.
    taskversion = models.CharField(
        db_column='TaskVersion', max_length=40, blank=True, null=True)
    parentjobid = models.ForeignKey('self', on_delete=models.CASCADE, db_column='ParentJobID',
                                    related_name='childrenofjob', null=True, blank=True)  # Field name made lowercase.
    preceedingjobid = models.ForeignKey('self', on_delete=models.SET_NULL,   db_column='PreceedingJobID',
                                        related_name='succeedingjobs', null=True, blank=True)  # Field name made lowercase.
    # Field name made lowercase.
    treeleft = models.IntegerField(db_column='TreeLeft', blank=True, null=True)
    # Field name made lowercase.
    treeright = models.IntegerField(
        db_column='TreeRight', blank=True, null=True)
    userid = models.ForeignKey('Users', on_delete=models.SET_NULL, null=True,  db_column='UserId',
                               related_name='jobsofuser', blank=True)  # Field name made lowercase.
    # Field name made lowercase.
    processid = models.IntegerField(
        db_column='ProcessId', blank=True, null=True)
    # Field name made lowercase.
    # MN Comment this field out for legacy database access
    #processcharid = models.CharField(
    #    db_column='ProcessCharId', blank=True, null=True, max_length=128)

    class Meta:
        managed = True
        db_table = 'Jobs'
        ordering = ['-creationtime']
        app_label = 'CCP4i2'

    def __unicode__(self):
        return '<djangoDB::{} {} {} "{}" {} {}>'.format(self.__class__.__name__, self.projectid.projectname, self.jobnumber, self.status.statustext, self.jobtitle, self.taskname)

    def __str__(self):
        return '<djangoDB::{} {} {} "{}" {} {}>'.format(self.__class__.__name__, self.projectid.projectname, self.jobnumber, self.status.statustext, self.jobtitle, self.taskname)

    @property
    def jobDirectory(self):
        pathElements = [os.path.normpath(self.projectid.projectdirectory), 'CCP4_JOBS'] + [
            'job_'+digit for digit in self.jobnumber.split('.')]
        return os.path.sep.join(pathElements)

    def descendents(self, growingList=None):
        if growingList is None:
            growingList = []
        for descendent in self.childrenofjob.all():
            growingList.append(descendent)
            descendent.descendents(growingList)
        return growingList

    def gracefulDelete(self):
        if self.status.statustext not in ['Finished', 'Failed', 'Unsatisfactory', 'Pending', 'Interrupted']:
            raise Exception(
                f"Unable to gracefullly delete because self.status.statustext is {self.status.statustext}")
        jobDirectory = self.jobDirectory
        normAbsPath = os.path.normpath(jobDirectory)
        if not normAbsPath.startswith(os.path.split(settings.CCP4I2_PROJECTS_DIR)[0]):
            raise Exception("normAbsPath ({}) does not start with {}".format(
                normAbsPath, os.path.split(settings.CCP4I2_PROJECTS_DIR)[0]))
        if not 'CCP4_JOBS' in normAbsPath:
            raise Exception("'CCP4_JOBS' not in normAbsPath ")
        if not os.path.sep.join(['job_{}'.format(element) for element in self.jobnumber.split('.')]) in normAbsPath:
            raise Exception("'{}' not in normAbsPath ".format(os.path.sep.join(
                ['job_{}'.format(element) for element in self.jobnumber.split('.')])))

        import shutil
        for descendent in self.childrenofjob.all():
            if descendent.status.statustext not in ['Finished', 'Failed', 'Unsatisfactory', 'Pending', 'Interrupted']:
                raise Exception(
                    f"Unable to gracefullly delete because descendent.status.statustext is {descendent.status.statustext}")
        for descendent in self.childrenofjob.all():
            descendent.gracefulDelete()
        try:
            shutil.rmtree(self.jobDirectory)
        except FileNotFoundError as err:
            logger.warning(f"Missing file when trying to delete {self.projectid.projectname} Job {self.jobnumber}")
            
        self.delete()

    def dependents(self, growingList=None):
        if growingList is None:
            growingList = []
        filesOfJob = self.filesofjob.all()
        fileUses = Fileuses.objects.filter(
            jobid__parentjobid__isnull=True).filter(fileid__in=filesOfJob)
        newJobs = list(
            set([aFileUse.jobid for aFileUse in fileUses]).difference(set(growingList)))
        for newJob in newJobs:
            growingList.append(newJob)
            newJob.dependents(growingList)
        return growingList

    def fail(self):
        self.status = Jobstatus.objects.get(
            statustext='Failed')
        self.save()


class Lastuniqueids(models.Model):
    # Field name made lowercase.
    tablename = models.CharField(
        db_column='TableName', unique=True, max_length=50, primary_key=True)
    # Field name made lowercase. This field type is a guess.
    lastid = models.TextField(db_column='LastId', blank=True, null=True)

    class Meta:
        managed = True
        db_table = 'LastUniqueIds'
        app_label = 'CCP4i2'


class Projectcomments(models.Model):
    # Field name made lowercase.
    projectcommentid = models.CharField(
        db_column='ProjectCommentID', unique=True, primary_key=True,  default=generate_ccp4_uuid, max_length=40)
    # Field name made lowercase.
    projectid = models.ForeignKey('Projects', on_delete=models.CASCADE,
                                  db_column='ProjectID', related_name='comments', null=True)
    userid = models.ForeignKey('Users', on_delete=models.SET_NULL, null=True,  db_column='UserID',
                               related_name='projectcommentsofuser')  # Field name made lowercase.
    # Field name made lowercase. This field type is a guess.
    timeofcomment = models.TextField(db_column='TimeOfComment')
    # Field name made lowercase.
    comment = models.TextField(db_column='Comment', blank=True, null=True)

    class Meta:
        managed = True
        db_table = 'ProjectComments'
        unique_together = (('projectid', 'projectcommentid'),)
        app_label = 'CCP4i2'


class Projectexports(models.Model):
    # Field name made lowercase.
    projectexportid = models.CharField(
        db_column='ProjectExportId', unique=True, primary_key=True,  default=generate_ccp4_uuid, max_length=40)
    projectid = models.ForeignKey('Projects', on_delete=models.CASCADE, db_column='ProjectID',
                                  related_name='exportsofproject', null=True)  # Field name made lowercase.
    # Field name made lowercase. This field type is a guess.
    projectexporttime = models.TextField(db_column='ProjectExportTime')
    # Field name made lowercase. This field type is a guess.
    projectexportafter = models.TextField(
        db_column='ProjectExportAfter', blank=True, null=True)

    class Meta:
        managed = True
        db_table = 'ProjectExports'
        app_label = 'CCP4i2'


class Projectimports(models.Model):
    # Field name made lowercase.
    projectimportid = models.CharField(
        db_column='ProjectImportId', unique=True, primary_key=True,  default=generate_ccp4_uuid, max_length=40)
    projectid = models.ForeignKey('Projects', on_delete=models.CASCADE, db_column='ProjectID',
                                  related_name='importsofproject', null=True)  # Field name made lowercase.
    # Field name made lowercase. This field type is a guess.
    projectimporttime = models.TextField(db_column='ProjectImportTime')
    projectexportid = models.ForeignKey(Projectexports, on_delete=models.SET_NULL, null=True,
                                        db_column='ProjectExportId', related_name='projectimportsofexport')  # Field name made lowercase.
    # Field name made lowercase.
    projectexportdatabaseid = models.CharField(
        db_column='ProjectExportDatabaseId', max_length=40)
    # Field name made lowercase. This field type is a guess.
    projectexporttime = models.TextField(
        db_column='ProjectExportTime', blank=True, null=True)
    # Field name made lowercase. This field type is a guess.
    projectexportafter = models.TextField(
        db_column='ProjectExportAfter', blank=True, null=True)

    class Meta:
        managed = True
        db_table = 'ProjectImports'
        app_label = 'CCP4i2'


def nextProjecttags():
    if Projecttags.objects.count() == 0:
        return 1
    return Projecttags.objects.only('projecttagsid').aggregate(Max('projecttagsid'))['projecttagsid__max']+1


class Projecttags(models.Model):  # Here seeing if rowid can serve as pk
    # Field name made lowercase.
    projecttagsid = models.IntegerField(
        db_column='rowid', unique=True, primary_key=True, default=nextProjecttags)
    tagid = models.ForeignKey('Tags', on_delete=models.CASCADE, db_column='TagID',
                              related_name='projecttagsoftag', null=True)  # Field name made lowercase.
    # Field name made lowercase.
    projectid = models.ForeignKey('Projects', on_delete=models.CASCADE,
                                  db_column='ProjectID', related_name='projecttagsofproject')

    class Meta:
        managed = True
        db_table = 'ProjectTags'
        unique_together = (('projectid', 'tagid'),)
        ordering = ['projecttagsid']
        app_label = 'CCP4i2'


class Projects(models.Model):
    projectid = models.CharField(db_column='ProjectID', unique=True, primary_key=True,
                                 default=generate_ccp4_uuid, max_length=40)  # Field name made lowercase.
    # Field name made lowercase.
    projectname = models.CharField(
        db_column='ProjectName', unique=True, max_length=100)
    # Field name made lowercase.
    projectcreated = models.FloatField(
        db_column='ProjectCreated', default=time.time)
    # Field name made lowercase.
    userid = models.ForeignKey('Users', on_delete=models.SET_NULL,
                               null=True, db_column='UserID', related_name='projectsofuser')
    parentprojectid = models.ForeignKey('self', on_delete=models.SET_NULL, null=True,
                                        db_column='ParentProjectID', related_name='childprojects')  # Field name made lowercase.
    # Field name made lowercase.
    projectdirectory = models.TextField(
        db_column='ProjectDirectory', blank=True, null=True)
    # Field name made lowercase. This field type is a guess.
    lastjobnumber = models.TextField(
        db_column='LastJobNumber', blank=True, null=True)
    followfromjobid = models.ForeignKey(Jobs, on_delete=models.SET_NULL, null=True,  blank=True,
                                        db_column='FollowFromJobID', related_name='followingjobs')  # Field name made lowercase.
    # Field name made lowercase.
    lastcleanuptime = models.FloatField(
        db_column='LastCleanupTime', blank=True, null=True)
    # Field name made lowercase.
    i1projectname = models.CharField(
        db_column='I1ProjectName', max_length=100, blank=True, null=True)
    # Field name made lowercase.
    i1projectdirectory = models.TextField(
        db_column='I1ProjectDirectory', blank=True, null=True)
    # Field name made lowercase.
    lastaccess = models.FloatField(
        db_column='LastAccess', blank=True, null=True)

    class Meta:
        managed = True
        db_table = 'Projects'
        app_label = 'CCP4i2'

    def __unicode__(self):
        return '<djangoDB::{} {} {}>'.format(self.__class__.__name__, self.projectid, self.projectname)

    def __str__(self):
        return '<djangoDB::{} {} {}>'.format(self.__class__.__name__, self.projectid, self.projectname)


@receiver(pre_delete, sender=Projects)
def pre_ProjectDelete(sender, instance, **kwargs):
    try:
        print("Would have to delete", instance.projectdirectory)
        if instance.projectdirectory is not None:
            projectDirectory = pathlib.Path(instance.projectdirectory)
            assert projectDirectory.parent.name == "CCP4I2_PROJECTS"
            assert projectDirectory.exists()
            assert (projectDirectory.name == slugify(instance.projectname)
                    ) or projectDirectory.name == instance.projectname
            shutil.rmtree(projectDirectory)
    except AssertionError as err:
        print(f'Unsafe to delete project {instance.projectname} Because {err}')
        raise err


def nextProjectsuserspermissions():
    if Projectsuserspermissions.objects.count() == 0:
        return 1
    return Projectsuserspermissions.objects.only('projectsuserspermissionsid').aggregate(Max('projectsuserspermissionsid'))['projectsuserspermissionsid__max']+1


# Here seeing if rowid can serve as pk
class Projectsuserspermissions(models.Model):
    projectsuserspermissionsid = models.IntegerField(
        db_column='rowid', unique=True, primary_key=True, default=nextProjectsuserspermissions)  # Field name made lowercase.
    projectid = models.ForeignKey(Projects, on_delete=models.CASCADE, db_column='ProjectID',
                                  related_name='projectsuserspermissionsofproject', null=True)  # Field name made lowercase.
    userid = models.ForeignKey('Users', on_delete=models.CASCADE, db_column='UserID',
                               related_name='projectsuserspermissionsofuser', null=True)  # Field name made lowercase.
    privilegeid = models.ForeignKey('Userprivileges', on_delete=models.CASCADE, db_column='PrivilegeID',
                                    related_name='projectsuserspermissionsofprivilege', null=True)  # Field name made lowercase.
    creatorid = models.ForeignKey('Users', on_delete=models.SET_NULL, null=True,  db_column='CreatorID',
                                  blank=True, related_name='projectsuserspermissionscreated')  # Field name made lowercase.

    class Meta:
        managed = True
        db_table = 'ProjectsUsersPermissions'
        unique_together = (('projectid', 'userid'),)
        ordering = ['projectsuserspermissionsid']
        app_label = 'CCP4i2'


class Serverjobs(models.Model):
    # Field name made lowercase.
    jobid = models.OneToOneField(
        Jobs, models.CASCADE, db_column='JobId', primary_key=True, null=False, blank=True)
    # Field name made lowercase.
    serverprocessid = models.IntegerField(
        db_column='ServerProcessId', blank=True, null=True)
    # Field name made lowercase.
    machine = models.CharField(
        db_column='Machine', max_length=255, blank=True, null=True)
    # Field name made lowercase.
    username = models.CharField(
        db_column='Username', max_length=100, blank=True, null=True)
    # Field name made lowercase.
    mechanism = models.CharField(
        db_column='Mechanism', max_length=40, blank=True, null=True)
    # Field name made lowercase.
    remotepath = models.CharField(
        db_column='RemotePath', max_length=255, blank=True, null=True)
    # Field name made lowercase.
    customcodefile = models.CharField(
        db_column='CustomCodeFile', max_length=255, blank=True, null=True)
    # Field name made lowercase.
    validate = models.CharField(
        db_column='Validate', max_length=40, blank=True, null=True)
    # Field name made lowercase.
    keyfilename = models.CharField(
        db_column='KeyFilename', max_length=255, blank=True, null=True)
    # Field name made lowercase.
    servergroup = models.CharField(
        db_column='ServerGroup', max_length=40, blank=True, null=True)

    class Meta:
        managed = True
        db_table = 'ServerJobs'
        app_label = 'CCP4i2'


class Tags(models.Model):
    tagid = models.CharField(db_column='TagID', unique=True, primary_key=True,
                             default=generate_ccp4_uuid, max_length=40)  # Field name made lowercase.
    # Field name made lowercase.
    parenttagid = models.ForeignKey(
        'self', on_delete=models.CASCADE, db_column='ParentTagID', null=True, related_name='childtags')
    # Field name made lowercase.
    text = models.CharField(db_column='Text', max_length=50, unique=True)

    class Meta:
        managed = True
        db_table = 'Tags'
        app_label = 'CCP4i2'


class Useragents(models.Model):
    # Field name made lowercase. This field type is a guess.
    useragentid = models.IntegerField(
        db_column='UserAgentID', unique=True, primary_key=True)
    # Field name made lowercase.
    useragentname = models.CharField(
        db_column='UserAgentName', unique=True, max_length=255)
    # Field name made lowercase.
    version = models.CharField(
        db_column='Version', max_length=50, blank=True, null=True)

    class Meta:
        managed = True
        db_table = 'UserAgents'
        ordering = ['useragentid']
        app_label = 'CCP4i2'

    def __unicode__(self):
        return '<djangoDB::{} {} {} {}>'.format(self.__class__.__name__, self.useragentid, self.useragentname, self.version)


class Userprivileges(models.Model):
    # Field name made lowercase. This field type is a guess.
    privilegeid = models.IntegerField(
        db_column='PrivilegeID', unique=True, primary_key=True)
    # Field name made lowercase.
    privilegetext = models.CharField(db_column='PrivilegeText', max_length=100)

    class Meta:
        managed = True
        db_table = 'UserPrivileges'
        ordering = ['privilegeid']
        app_label = 'CCP4i2'

    def __unicode__(self):
        return '<djangoDB::{} {} {}>'.format(self.__class__.__name__, self.privilegeid, self.privilegetext)


class Userroles(models.Model):
    # Field name made lowercase. This field type is a guess.
    userroleid = models.IntegerField(
        db_column='UserRoleID', unique=True, primary_key=True)
    # Field name made lowercase.
    userrole = models.CharField(
        db_column='UserRole', max_length=255, blank=True, null=True)

    class Meta:
        managed = True
        db_table = 'UserRoles'
        ordering = ['userroleid']
        app_label = 'CCP4i2'

    def __unicode__(self):
        return '<djangoDB::{} {} "{}">'.format(self.__class__.__name__, self.userroleid, self.userrole)


class Users(models.Model):
    userid = models.CharField(db_column='UserId', unique=True, primary_key=True,
                              default=generate_ccp4_uuid, max_length=40)  # Field name made lowercase.
    # Field name made lowercase.
    username = models.CharField(
        db_column='UserName', unique=True, max_length=100)
    # Field name made lowercase.
    userpassword = models.CharField(
        db_column='UserPassword', max_length=100, blank=True, null=True)
    userroleid = models.ForeignKey(Userroles, on_delete=models.SET_NULL, null=True,
                                   db_column='UserRoleID', related_name='userrolesofuser')  # Field name made lowercase.

    class Meta:
        managed = True
        db_table = 'Users'
        app_label = 'CCP4i2'

    def __unicode__(self):
        return '<djangoDB::{} {} {} role:{}>'.format(self.__class__.__name__, self.userid, self.username, self.userroleid.userrole)

    @classmethod
    def ownerUser(cls, username=None, password=""):
        return cls(username=username, userpassword="", userroleid=Userroles.objects.get(userrole="Owner"))


class Xdata(models.Model):
    xdataid = models.CharField(db_column='XDataID', unique=True, primary_key=True,
                               default=generate_ccp4_uuid, max_length=40)  # Field name made lowercase.
    # Field name made lowercase.
    xdataclass = models.TextField(db_column='XDataClass')
    # Field name made lowercase.
    xdataxml = models.TextField(db_column='XDataXml')
    jobid = models.ForeignKey(Jobs, on_delete=models.CASCADE, db_column='JobID',
                              related_name='xdatasofjob', null=True, blank=True)  # Field name made lowercase.

    class Meta:
        managed = True
        db_table = 'XData'
        app_label = 'CCP4i2'
