# Database documentation

## Mapping from the old i2 SQL schema to the new Django models

| Old                                                 | New                          |
| --------------------------------------------------- | ---------------------------- |
| Comments.Comment                                    | Job.comments                 |
| Comments.CommentID                                  |                              |
| Comments.JobID                                      |                              |
| Comments.TimeOfComment                              |                              |
| Comments.UserID                                     |                              |
| Databases.CreationHostName                          | Project.creation_host        |
| Databases.CreationTime                              |                              |
| Databases.CreationUserName                          | Project.creation_user        |
| Databases.DatabaseId                                |                              |
| Databases.SchemaDate                                |                              |
| Databases.SchemaVersion                             |                              |
| Databases.ThisIsMe                                  |                              |
| DirectoryAliases.Directory                          |                              |
| DirectoryAliases.DirectoryAlias                     |                              |
| ExportFiles.Annotation                              |                              |
| ExportFiles.CreationTime                            | FileExport.time              |
| ExportFiles.ExportFilename                          | FileExport.name              |
| ExportFiles.ExportId                                |                              |
| ExportFiles.FileID                                  | FileExport.file              |
| FileAssociationMembers.FileAssociationID            |                              |
| FileAssociationMembers.FileAssociationRoleID        |                              |
| FileAssociationMembers.FileID                       |                              |
| FileAssociationRoles.FileAssociationRoleDescription |                              |
| FileAssociationRoles.FileAssociationRoleID          |                              |
| FileAssociationRoles.FileAssociationRoleName        |                              |
| FileAssociationRoles.FileAssociationTypeID          |                              |
| FileAssociations.FileAssociationID                  |                              |
| FileAssociations.FileAssociationName                |                              |
| FileAssociations.FileAssociationTypeId              |                              |
| FileAssociationTypes.FileAssociationTypeDescription |                              |
| FileAssociationTypes.FileAssociationTypeID          |                              |
| FileAssociationTypes.FileAssociationTypeName        |                              |
| FileRoles.RoleID                                    |                              |
| FileRoles.RoleText                                  |                              |
| Files.Annotation                                    | File.annotation              |
| Files.FileContent                                   | File.content                 |
| Files.FileID                                        | File.uuid                    |
| Files.Filename                                      | File.name                    |
| Files.FilePath                                      |                              |
| Files.FileSubType                                   | File.sub_type                |
| Files.FiletypeID                                    | File.type                    |
| Files.JobID                                         | File.job                     |
| Files.JobParamName                                  | File.job_param_name          |
| Files.PathFlag                                      | File.directory               |
| FileTypes.FileTypeDescription                       | FileType.description         |
| FileTypes.FileTypeID                                |                              |
| FileTypes.FileTypeName                              | FileType.name                |
| FileUses.FileID                                     | FileUse.file                 |
| FileUses.JobID                                      | FileUse.job                  |
| FileUses.JobParamName                               | FileUse.job_param_name       |
| FileUses.RoleID                                     | FileUse.role                 |
| ImportFiles.Annotation                              |                              |
| ImportFiles.Checksum                                | FileImport.checksum          |
| ImportFiles.CreationTime                            | FileImport.time              |
| ImportFiles.ExportFileID                            |                              |
| ImportFiles.FileID                                  | FileImport.file              |
| ImportFiles.ImportId                                |                              |
| ImportFiles.ImportNumber                            |                              |
| ImportFiles.LastModifiedTime                        | FileImport.last_modified     |
| ImportFiles.Reference                               |                              |
| ImportFiles.SourceFileID                            |                              |
| ImportFiles.SourceFilename                          | FileImport.name              |
| JobEvaluation.EvaluationID                          |                              |
| JobEvaluation.EvaluationText                        |                              |
| JobKeyCharValues.JobID                              | JobCharValue.job             |
| JobKeyCharValues.KeyTypeID                          | JobCharValue.key             |
| JobKeyCharValues.Value                              | JobCharValue.value           |
| JobKeyValues.JobID                                  | JobFloatValue.job            |
| JobKeyValues.KeyTypeID                              | JobFloatValue.key            |
| JobKeyValues.Value                                  | JobFloatValue.value          |
| Jobs.CreationTime                                   | Job.creation_time            |
| Jobs.Evaluation                                     | Job.evaluation               |
| Jobs.FinishTime                                     | Job.finish_time              |
| Jobs.JobID                                          | Job.uuid                     |
| Jobs.JobNumber                                      | Job.number                   |
| Jobs.JobTitle                                       | Job.title                    |
| Jobs.ParentJobID                                    | Job.parent                   |
| Jobs.PreceedingJobID                                |                              |
| Jobs.processId                                      | Job.process_id               |
| Jobs.ProjectID                                      | Job.project                  |
| Jobs.Status                                         | Job.status                   |
| Jobs.TaskName                                       | Job.task_name                |
| Jobs.TaskVersion                                    |                              |
| Jobs.TreeLeft                                       |                              |
| Jobs.TreeRight                                      |                              |
| Jobs.UserAgent                                      |                              |
| Jobs.UserId                                         |                              |
| JobStatus.StatusID                                  |                              |
| JobStatus.StatusText                                |                              |
| KeyTypes.KeyTypeDescription                         | JobValueKeys.description     |
| KeyTypes.KeyTypeID                                  |                              |
| KeyTypes.KeyTypeName                                | JobValueKeys.name            |
| LastUniqueIds.LastId                                |                              |
| LastUniqueIds.TableName                             |                              |
| ProjectComments.Comment                             | Project.description          |
| ProjectComments.ProjectCommentID                    |                              |
| ProjectComments.ProjectID                           |                              |
| ProjectComments.TimeOfComment                       |                              |
| ProjectComments.UserID                              |                              |
| ProjectExports.ProjectExportAfter                   |                              |
| ProjectExports.ProjectExportId                      |                              |
| ProjectExports.ProjectExportTime                    | ProjectExport.time           |
| ProjectExports.ProjectID                            | ProjectExport.project        |
| ProjectImports.ProjectExportAfter                   |                              |
| ProjectImports.ProjectExportDatabaseId              |                              |
| ProjectImports.ProjectExportId                      |                              |
| ProjectImports.ProjectExportTime                    |                              |
| ProjectImports.ProjectID                            | ProjectImport.project        |
| ProjectImports.ProjectImportId                      |                              |
| ProjectImports.ProjectImportTime                    | ProjectImport.time           |
| Projects.FollowFromJobID                            | Project.follow_from_job      |
| Projects.I1ProjectDirectory                         | Project.i1_project_directory |
| Projects.I1ProjectName                              | Project.i1_project_name      |
| Projects.LastAccess                                 | Project.last_access          |
| Projects.LastCleanupTime                            |                              |
| Projects.LastJobNumber                              | Project.last_job_number      |
| Projects.ParentProjectID                            |                              |
| Projects.ProjectCreated                             | Project.creation_time        |
| Projects.ProjectDirectory                           | Project.directory            |
| Projects.ProjectID                                  | Project.uuid                 |
| Projects.ProjectName                                | Project.name                 |
| Projects.UserID                                     |                              |
| ProjectsUsersPermissions.CreatorID                  |                              |
| ProjectsUsersPermissions.PrivilegeID                |                              |
| ProjectsUsersPermissions.ProjectID                  |                              |
| ProjectsUsersPermissions.UserID                     |                              |
| ProjectTags.ProjectID                               | ProjectTag.projects          |
| ProjectTags.TagID                                   |                              |
| ServerJobs.CustomCodeFile                           | ServerJob.custom_code_file   |
| ServerJobs.JobId                                    | ServerJob.job                |
| ServerJobs.KeyFilename                              | ServerJob.key_file_name      |
| ServerJobs.Machine                                  | ServerJob.machine            |
| ServerJobs.Mechanism                                | ServerJob.mechanism          |
| ServerJobs.RemotePath                               | ServerJob.remote_path        |
| ServerJobs.ServerGroup                              | ServerJob.server_group       |
| ServerJobs.ServerProcessId                          | ServerJob.server_process_id  |
| ServerJobs.Username                                 | ServerJob.username           |
| ServerJobs.Validate                                 | ServerJob.validate           |
| Tags.ParentTagID                                    | ProjectTag.parent            |
| Tags.TagID                                          |                              |
| Tags.Text                                           | ProjectTag.text              |
| UserAgents.UserAgentID                              |                              |
| UserAgents.UserAgentName                            |                              |
| UserAgents.Version                                  |                              |
| UserPrivileges.PrivilegeID                          |                              |
| UserPrivileges.PrivilegeText                        |                              |
| UserRoles.UserRole                                  |                              |
| UserRoles.UserRoleID                                |                              |
| Users.UserID                                        |                              |
| Users.UserName                                      |                              |
| Users.UserPassword                                  |                              |
| Users.UserRoleID                                    |                              |
| XData.JobID                                         | XData.job                    |
| XData.XDataClass                                    | XData.class                  |
| XData.XDataID                                       |                              |
| XData.XDataXml                                      | XData.xml                    |
