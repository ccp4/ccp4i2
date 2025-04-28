-- CCP4I2 SCHEMA
-- DO NOT REMOVE OR CHANGE FORMAT OF VERSION/DATE LINES!
-- THE VERSION/DATE SHOULD BE INCREMENTED WHENEVER THE SCHEMA IS CHANGED
-- VERSION 0.1.22
-- DATE 23-09-2016

-- !!!!BEWARE CHANGES TO MANY TABLES MUST BE REPLICATED IN temp_datase_schema.sql!!!!

-- 0.1.1 02-03-2012 Added ImportFiles table
-- 0.1.2 08-03-2012 Added ExportFiles table
-- 0.1.3 20-04-2012 Changes to ImportFiles/ExportFiles
-- 0.1.4 14-09-2012 Change to use UUIDs as VARCHAR(32)
-- 0.1.5 10-10-2012 Remove Files.ImportId, add ImportFiles.FileId, add 'ON DELETE' keywords
-- 0.1.6 25-10-2012 Add Files.FileSubType and Files.FileContent. Remove Files.RelPath
-- 0.1.7 21-11-2012 Add refmac-dictionary to FileTypes table
-- 0.1.8 02-01-2013 Add lastCleanupTime to project table
-- 0.1.9 06-04-2013 Remove ON DELETE CASCADE. Remove Jobs.machine. Add Jobs.TaskVersion
-- 0.1.10 17-05-2013 Convert DbVersion table to Databases table and add ImportProject and ExportProject table
-- 0.1.11 23-07-2013 Add JobParamName to Files  and Fileuses table
-- 0.1.12 14-12-2013 Add ModifiedTime to ImportFiles table
-- 0.1.13 06-02-2014 AddChecksum to ImportFiles table
-- 0.1.14 07-02-2014 Add ImportNumber/ImportReference to ImportFiles table
-- 0.1.15 15-05-2014 Add UserId to Jobs table, Key/value table and KeyValueTypes table, 
--                    association table and association members table
-- 0.1.16 19-08-2014 Add JobKeyCharValue table
-- 0.1.17 30-01-2015 Add Files.FilePath
-- 0.1.18 01-03-2016 Add Projects.I1ProjectName and Projects.I1ProjectDirectory
-- 0.1.19 09-06-2016 Add Projects.LastAccess 
-- 0.1.20 02-08-2016 Add Jobs.processId and ServerParams table
-- 0.1.21 18-08-2016 Add  ServerParams.ServerGroup
-- 0.1.22 23-09-2016 Add Tags, ProjectTags and ProjectComments tables

CREATE TABLE Databases
(
	DatabaseId VARCHAR(32) NOT NULL,
	ThisIsMe SHORTINT NOT NULL DEFAULT 0,
	SchemaVersion VARCHAR(50),
	SchemaDate VARCHAR(50),
	CreationTime REAL,
	CreationHostName VARCHAR(255),
	CreationUserName VARCHAR(100),
	PRIMARY KEY (DatabaseId),
	UNIQUE (DatabaseId)

);

CREATE TABLE UserAgents
(
	UserAgentID MEDIUMINT NOT NULL,
	UserAgentName VARCHAR(255) NOT NULL,
	Version VARCHAR(50),
	PRIMARY KEY (UserAgentID),
	UNIQUE (UserAgentID),
	UNIQUE (UserAgentName)
) ;

CREATE TABLE UserRoles
(
	UserRoleID MEDIUMINT NOT NULL,
	UserRole VARCHAR(255),
	UNIQUE (UserRoleID),
	PRIMARY KEY (UserRoleID)
);

CREATE TABLE Users
(
	UserID VARCHAR(32) NOT NULL,
	UserName VARCHAR(100) NOT NULL,
        UserPassword VARCHAR(100),
        UserRoleID MEDIUMINT NOT NULL,
	PRIMARY KEY (UserID),
	UNIQUE (UserID),
	UNIQUE (UserName),
        FOREIGN KEY (UserRoleID) REFERENCES UserRoles (UserRoleID)
);

CREATE TABLE UserPrivileges
(
	PrivilegeID MEDIUMINT NOT NULL,
	PrivilegeText VARCHAR(100) NOT NULL,
	PRIMARY KEY (PrivilegeID),
	UNIQUE (PrivilegeID)
) ;


CREATE TABLE Projects
(
	ProjectID VARCHAR(32) NOT NULL,
	ProjectName VARCHAR(100) NOT NULL,
	ProjectCreated REAL NOT NULL,
	UserID MEDIUMINT NOT NULL,
	ParentProjectID VARCHAR(32),
        ProjectDirectory TEXT,
	LastJobNumber MEDIUMINT,
	FollowFromJobID VARCHAR(32),
        LastCleanupTime REAL,
	I1ProjectName  VARCHAR(100),
	I1ProjectDirectory TEXT,
	LastAccess REAL,
	PRIMARY KEY (ProjectID),
	UNIQUE (ProjectID),
	UNIQUE (ProjectName),
	FOREIGN KEY (UserID) REFERENCES Users (UserID),
	FOREIGN KEY (ParentProjectID) REFERENCES Projects (ProjectID),
	FOREIGN KEY (FollowFromJobID) REFERENCES Jobs (JobID) ON DELETE SET NULL
) ;

CREATE TABLE DirectoryAliases
(
	DirectoryAlias VARCHAR(100) NOT NULL,
        Directory TEXT,
	PRIMARY KEY (DirectoryAlias),
	UNIQUE (DirectoryAlias)
) ;


CREATE TABLE FileTypes
(
	FileTypeID MEDIUMINT NOT NULL,
	FileTypeName VARCHAR(50) NOT NULL,
	FileTypeDescription TEXT,
	PRIMARY KEY (FileTypeID),
	UNIQUE (FileTypeID),
	UNIQUE (FileTypeName)
) ;

CREATE TABLE Files
(
	FileID VARCHAR(32) NOT NULL,
	Filename TEXT NOT NULL,
	Annotation TEXT,
	FiletypeID MEDIUMINT NOT NULL,
        FileSubType MEDIUMINT,
	FileContent MEDIUMINT,
	FilePath MEDIUMINT,
	JobID VARCHAR(32),
	JobParamName VARCHAR(32),
	PathFlag MEDIUMINT,
	PRIMARY KEY (FileID),
	UNIQUE (FileID),
        FOREIGN KEY (FiletypeID) REFERENCES Filetypes (FiletypeID),
--        FOREIGN KEY (JobID) REFERENCES Jobs (JobID) ON DELETE CASCADE
        FOREIGN KEY (JobID) REFERENCES Jobs (JobID)
) ;

CREATE TABLE ImportFiles
(
	ImportId VARCHAR(32) NOT NULL,
	FileID VARCHAR(32),
	SourceFilename TEXT,
	SourceFileID VARCHAR(32),
	ExportFileID VARCHAR(32),
	Annotation TEXT,
	CreationTime REAL,
	LastModifiedTime REAL,
	Checksum VARCHAR(32),
	ImportNumber MEDIUMINT,
	Reference TEXT,
	PRIMARY KEY (ImportId),
	UNIQUE (ImportId),
	FOREIGN KEY (FileID) REFERENCES Files (FileID),
	FOREIGN KEY (SourceFileID) REFERENCES Files (FileID),
	FOREIGN KEY (ExportFileID) REFERENCES ExportFiles (ExportID)
); 

CREATE TABLE ExportFiles
(
	ExportId VARCHAR(32) NOT NULL,
	FileID VARCHAR(32) NOT NULL,
	ExportFilename TEXT,
	Annotation TEXT,
	CreationTime REAL,
	PRIMARY KEY (ExportId),	
	UNIQUE (ExportId),
	FOREIGN KEY (FileID) REFERENCES Files (FileID)
); 

CREATE TABLE XData
(
	XDataID VARCHAR(32) NOT NULL,
	XDataClass TEXT NOT NULL,
	XDataXml  TEXT NOT NULL,
	JobID  VARCHAR(32),
	PRIMARY KEY (XDataID),
	UNIQUE (XDataID),
--	FOREIGN KEY (JobID) REFERENCES Jobs (JobID) ON DELETE CASCADE
	FOREIGN KEY (JobID) REFERENCES Jobs (JobID)
);

CREATE TABLE JobStatus
(
	StatusID MEDIUMINT NOT NULL,
	StatusText VARCHAR(255) NOT NULL,
	PRIMARY KEY (StatusID),
	UNIQUE (StatusID),
	UNIQUE (StatusText)
);

CREATE TABLE JobEvaluation
(
        EvaluationID MEDIUMINT NOT NULL,
        EvaluationText VARCHAR(255) NOT NULL,
        PRIMARY KEY (EvaluationID),
        UNIQUE (EvaluationID)
);

CREATE TABLE Jobs
(
	JobID VARCHAR(32) NOT NULL,
	JobNumber VARCHAR(50) NOT NULL,
	CreationTime REAL NOT NULL,
	FinishTime REAL,
	Status MEDIUMINT NOT NULL,
        Evaluation MEDIUMINT,
        UserId VARCHAR(32),
	UserAgent MEDIUMINT NOT NULL,
	JobTitle VARCHAR(255),
	ProjectID VARCHAR(32) NOT NULL,
	TaskName VARCHAR(100) NOT NULL,
	TaskVersion VARCHAR(32),
	ParentJobID VARCHAR(32),
	PreceedingJobID VARCHAR(32),
        TreeLeft INT,
	TreeRight INT,
	processId INT,
	PRIMARY KEY (JobID),
	UNIQUE (JobID),
	FOREIGN KEY (ProjectID) REFERENCES Projects (ProjectID),
--	FOREIGN KEY (ParentJobID) REFERENCES Jobs (JobID) ON DELETE CASCADE,
	FOREIGN KEY (ParentJobID) REFERENCES Jobs (JobID),
	FOREIGN KEY (PreceedingJobID) REFERENCES Jobs (JobID) ON DELETE SET NULL,
	FOREIGN KEY (Status) REFERENCES JobStatus  (StatusID),
	FOREIGN KEY (Evaluation) REFERENCES JobEvaluation  (EvaluationID),
        FOREIGN KEY (UserId) REFERENCES Users (UserId),
	FOREIGN KEY (UserAgent) REFERENCES UserAgents  (UserAgentID)
) ;


CREATE TABLE ServerJobs
(
	JobId VARCHAR(32) NOT NULL,
	ServerProcessId INT,
	Machine VARCHAR(255),
	Username  VARCHAR(100),
	Mechanism VARCHAR(32),
	RemotePath VARCHAR(255),
	CustomCodeFile  VARCHAR(255),
	Validate VARCHAR(32),
	KeyFilename  VARCHAR(255),
	ServerGroup VARCHAR(32),
	UNIQUE (JobId),
	PRIMARY KEY (JobId),
        FOREIGN KEY (JobId) REFERENCES Jobs (JobId)

) ;


CREATE TABLE Comments
(
	CommentID VARCHAR(32) NOT NULL,
	JobID VARCHAR(32) NOT NULL,
	UserID VARCHAR(255),
	TimeOfComment FLOAT NOT NULL,
	Comment TEXT,
	PRIMARY KEY (JobID, CommentID),
	UNIQUE  (CommentID),
--        FOREIGN KEY (JobID) REFERENCES Jobs (JobID) ON DELETE CASCADE
        FOREIGN KEY (JobID) REFERENCES Jobs (JobID)
) ;

CREATE TABLE FileRoles
(
	RoleID MEDIUMINT NOT NULL,
	RoleText VARCHAR(50) NOT NULL,
	PRIMARY KEY (RoleID),
	UNIQUE (RoleID),
	UNIQUE (RoleText)
) ;

CREATE TABLE FileUses
(
	FileID VARCHAR(32) NOT NULL,
	JobID VARCHAR(32) NOT NULL,
	RoleID MEDIUMINT NOT NULL,
        JobParamName VARCHAR(32),
	PRIMARY KEY (FileID, JobID, RoleID),
--	FOREIGN KEY (JobID) REFERENCES Jobs (JobID) ON DELETE CASCADE,
	FOREIGN KEY (JobID) REFERENCES Jobs (JobID),
	FOREIGN KEY (FileID) REFERENCES Files (FileID),
	FOREIGN KEY (RoleID) REFERENCES FileRoles (RoleID)
) ;


CREATE TABLE ProjectsUsersPermissions
(
	ProjectID VARCHAR(32) NOT NULL,
	UserID VARCHAR(32) NOT NULL,
	PrivilegeID MEDIUMINT NOT NULL,
	CreatorID VARCHAR(32),
	PRIMARY KEY (ProjectID, UserID),
	FOREIGN KEY (ProjectID) REFERENCES Projects (ProjectID),
	FOREIGN KEY (UserID) REFERENCES Users (UserID),
	FOREIGN KEY (CreatorID) REFERENCES Users (UserID),
	FOREIGN KEY (PrivilegeID) REFERENCES UserPrivileges (PrivilegeID)
) ;


CREATE TABLE LastUniqueIds
(
	TableName VARCHAR(50),
	LastId MEDIUMINT 
);


CREATE TABLE ProjectExports
(
	ProjectExportId VARCHAR(32) NOT NULL,
	ProjectID VARCHAR(32) NOT NULL,
	ProjectExportTime FLOAT NOT NULL,
	ProjectExportAfter FLOAT,
	--ProjectExportSelection -- need to handle this
	PRIMARY KEY (ProjectExportId),
	UNIQUE (ProjectExportId),
	FOREIGN KEY (ProjectID) REFERENCES Projects (ProjectID)
);

CREATE TABLE ProjectImports
(
	ProjectImportId VARCHAR(32) NOT NULL,
	ProjectID VARCHAR(32) NOT NULL,
	ProjectImportTime FLOAT NOT NULL,
	ProjectExportId VARCHAR(32),
	ProjectExportDatabaseId VARCHAR(32) NOT NULL,
	ProjectExportTime FLOAT,
	ProjectExportAfter FLOAT,
	PRIMARY KEY (ProjectImportId),
	UNIQUE (ProjectImportId),
	FOREIGN KEY (ProjectID) REFERENCES Projects (ProjectID),
	FOREIGN KEY (ProjectExportDatabaseId) REFERENCES Databases (DatabaseID)
);

CREATE TABLE KeyTypes (
        KeyTypeID MEDIUMINT NOT NULL,
        KeyTypeName VARCHAR(50) NOT NULL,
	KeyTypeDescription TEXT,
	PRIMARY KEY (KeyTypeID),
	UNIQUE (KeyTypeID),
	UNIQUE (KeyTypeName) );

CREATE TABLE JobKeyValues (
        JobID VARCHAR(32) NOT NULL,
        KeyTypeID MEDIUMINT NOT NULL,
        Value REAL,
        PRIMARY KEY (JobID,KeyTypeID),
        FOREIGN KEY (JobId)  REFERENCES Jobs (JobId),
        FOREIGN KEY (KeyTypeId)  REFERENCES KeyTypes (KeyTypeId) );

CREATE TABLE JobKeyCharValues (
        JobID VARCHAR(32) NOT NULL,
        KeyTypeID MEDIUMINT NOT NULL,
        Value  VARCHAR(255),
        PRIMARY KEY (JobID,KeyTypeID),
        FOREIGN KEY (JobId)  REFERENCES Jobs (JobId),
        FOREIGN KEY (KeyTypeId)  REFERENCES KeyTypes (KeyTypeId) );

CREATE TABLE FileAssociationTypes (
        FileAssociationTypeID MEDIUMINT NOT NULL,
        FileAssociationTypeName VARCHAR(50) NOT NULL,
	FileAssociationTypeDescription TEXT,
	PRIMARY KEY (FileAssociationTypeID),
	UNIQUE (FileAssociationTypeID),
	UNIQUE (FileAssociationTypeName) );

CREATE TABLE FileAssociationRoles (
        FileAssociationRoleID MEDIUMINT NOT NULL,
	FileAssociationTypeID MEDIUMINT NOT NULL,
        FileAssociationRoleName VARCHAR(50) NOT NULL,
	FileAssociationRoleDescription TEXT,
	PRIMARY KEY (FileAssociationRoleID),
	FOREIGN KEY (FileAssociationTypeID) REFERENCES FileAssociationTypes (FileAssociationTypeID),
	UNIQUE (FileAssociationRoleID),
	UNIQUE (FileAssociationRoleName) );

CREATE TABLE FileAssociations (
       FileAssociationID VARCHAR(32) NOT NULL,
       FileAssociationTypeId MEDIUMINT NOT NULL,
       FileAssociationName VARCHAR(50),
       PRIMARY KEY (FileAssociationId),
       UNIQUE (FileAssociationId) );

CREATE TABLE FileAssociationMembers (
       FileAssociationID VARCHAR(32) NOT NULL,
       FileID VARCHAR(32) NOT NULL,
       FileAssociationRoleID MEDIUMINT NOT NULL,
       PRIMARY KEY (FileID,FileAssociationID,FileAssociationRoleID),
       FOREIGN KEY (FileID) REFERENCES Files (FileID),
       FOREIGN KEY (FileAssociationID) REFERENCES FileAssociations (FileAssociationID) );

CREATE TABLE Tags (
       TagID VARCHAR(32) NOT NULL,
       ParentTagID VARCHAR(32),
       Text VARCHAR(50) NOT NULL,
       PRIMARY KEY (TagID),
       FOREIGN KEY (ParentTagID) REFERENCES Tags (TagID) );

CREATE TABLE ProjectTags (
       TagID VARCHAR(32) NOT NULL,
       ProjectID VARCHAR(32) NOT NULL,
       PRIMARY KEY ( ProjectID, TagID ),
       FOREIGN KEY (TagID) REFERENCES Tags (TagID) ,
       FOREIGN KEY (ProjectID) REFERENCES Projects (ProjectID) );
       
CREATE TABLE ProjectComments
(
	ProjectCommentID VARCHAR(32) NOT NULL,
	ProjectID VARCHAR(32) NOT NULL,
	UserID VARCHAR(255),
	TimeOfComment FLOAT NOT NULL,
	Comment TEXT,
	PRIMARY KEY (ProjectID, ProjectCommentID),
        FOREIGN KEY (ProjectID) REFERENCES Projects (ProjectID)
) ;

INSERT INTO UserAgents (UserAgentID,UserAgentName) VALUES (0, 'Unknown');
INSERT INTO UserAgents (UserAgentID,UserAgentName) VALUES (1, 'CCP4i2');
INSERT INTO UserAgents (UserAgentID,UserAgentName) VALUES (2, 'CCP4mg');
INSERT INTO UserAgents (UserAgentID,UserAgentName) VALUES (3, 'Coot');

INSERT INTO UserPrivileges (PrivilegeID,PrivilegeText) VALUES( 0, 'No access' );
INSERT INTO UserPrivileges (PrivilegeID,PrivilegeText) VALUES( 1, 'Read access' );
INSERT INTO UserPrivileges (PrivilegeID,PrivilegeText) VALUES( 2, 'Read and comment' );
INSERT INTO UserPrivileges (PrivilegeID,PrivilegeText) VALUES( 3, 'Read and write' );
INSERT INTO UserPrivileges (PrivilegeID,PrivilegeText) VALUES( 4, 'Read, write, delete');
INSERT INTO UserPrivileges (PrivilegeID,PrivilegeText) VALUES( 5, 'Create & delete projects');
INSERT INTO UserPrivileges (PrivilegeID,PrivilegeText) VALUES( 6, 'Extend access' );

INSERT INTO JobStatus (StatusID,StatusText) VALUES(0, 'Unknown' );
INSERT INTO JobStatus (StatusID,StatusText) VALUES(1, 'Pending' );
INSERT INTO JobStatus (StatusID,StatusText) VALUES(2, 'Queued' );
INSERT INTO JobStatus (StatusID,StatusText) VALUES(3, 'Running' );
INSERT INTO JobStatus (StatusID,StatusText) VALUES(4, 'Interrupted');
INSERT INTO JobStatus (StatusID,StatusText) VALUES(5, 'Failed' );
INSERT INTO JobStatus (StatusID,StatusText) VALUES(6, 'Finished' );
INSERT INTO JobStatus (StatusID,StatusText) VALUES(7, 'Running remotely' );

INSERT INTO JobEvaluation (EvaluationID,EvaluationText) VALUES(1, 'Best' );
INSERT INTO JobEvaluation (EvaluationID,EvaluationText) VALUES(2, 'Good' );
INSERT INTO JobEvaluation (EvaluationID,EvaluationText) VALUES(3, 'Rejected' );

INSERT INTO UserRoles (UserRoleID,UserRole) VALUES(0, 'Unknown' );
INSERT INTO UserRoles (UserRoleID,UserRole) VALUES(1, 'Manager' );
INSERT INTO UserRoles (UserRoleID,UserRole) VALUES(2, 'Owner' );
INSERT INTO UserRoles (UserRoleID,UserRole) VALUES(3, 'User' );

INSERT INTO FileRoles (RoleID,RoleText) VALUES(0, 'out' );
INSERT INTO FileRoles (RoleID,RoleText) VALUES(1, 'in' );

--INSERT INTO Filetypes (FiletypeID,FiletypeName,FiletypeDescription) VALUES(0,'Unknown','File type unknown');
--INSERT INTO Filetypes (FiletypeID,FiletypeName,FiletypeDescription) VALUES(1,'application/CCP4-seq','Sequence file');
--INSERT INTO Filetypes (FiletypeID,FiletypeName,FiletypeDescription) VALUES(2,'chemical/x-pdb','PDB or equivalent coordinate file');
--INSERT INTO Filetypes (FiletypeID,FiletypeName,FiletypeDescription) VALUES(3,'MultiPDB','Multiple PDB files');
--INSERT INTO Filetypes (FiletypeID,FiletypeName,FiletypeDescription) VALUES(4,'application/CCP4-mtz','MTZ of merged experimental data');
--INSERT INTO Filetypes (FiletypeID,FiletypeName,FiletypeDescription) VALUES(5,'application/CCP4-mtz-unmerged','MTZ of unmerged experimental data');
--INSERT INTO Filetypes (FiletypeID,FiletypeName,FiletypeDescription) VALUES(6,'application/CCP4-unmerged-experimental','Unmerged experimental data - not MTZ');
--INSERT INTO Filetypes (FiletypeID,FiletypeName,FiletypeDescription) VALUES(7,'application/CCP4-map','Electron density map');
--INSERT INTO Filetypes (FiletypeID,FiletypeName,FiletypeDescription) VALUES(8,'application/refmac-dictionary','Refmac dictionary');
--INSERT INTO Filetypes (FiletypeID,FiletypeName,FiletypeDescription) VALUES(9,'application/refmac-TLS','Refmac TLS');
--INSERT INTO Filetypes (FiletypeID,FiletypeName,FiletypeDescription) VALUES(10,'application/CCP4-mtz-freerflag','FreeR flag');
--INSERT INTO Filetypes (FiletypeID,FiletypeName,FiletypeDescription) VALUES(11,'application/CCP4-mtz-observed','Observed intensities and structure factors');
--INSERT INTO Filetypes (FiletypeID,FiletypeName,FiletypeDescription) VALUES(12,'application/CCP4-mtz-phases','Phases');
--INSERT INTO Filetypes (FiletypeID,FiletypeName,FiletypeDescription) VALUES(13,'application/CCP4-mtz-map','Map coefficients');
