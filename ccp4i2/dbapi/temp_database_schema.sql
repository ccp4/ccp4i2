-- CCP4I2 SCHEMA FOR TEMPORARY TABLES TO HANDLE IMPORTING DATA
-- DO NOT REMOVE OR CHANGE FORMAT OF VERSION/DATE LINES!
-- THE VERSION/DATE SHOULD BE INCREMENTED WHENEVER THE SCHEMA IS CHANGED
-- THIS SHOULD BE IN SYNC WITH DATABASE_SCHEMA.SQL VERSION
-- VERSION 0.1.10
-- DATE 17-05-2013


CREATE TABLE IF NOT EXISTS  TempFiles
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
	Excl INT DEFAULT 0,
	PRIMARY KEY (FileID),
	UNIQUE (FileID),
        FOREIGN KEY (FiletypeID) REFERENCES Filetypes (FiletypeID)
        --FOREIGN KEY (JobID) REFERENCES Jobs (JobID)
) ;

CREATE TABLE IF NOT EXISTS  TempImportFiles
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
	Excl INT DEFAULT 0,
	PRIMARY KEY (ImportId),
	UNIQUE (ImportId)
	--FOREIGN KEY (FileID) REFERENCES Files (FileID),
	--FOREIGN KEY (SourceFileID) REFERENCES Files (FileID),
	--FOREIGN KEY (ExportFileID) REFERENCES ExportFiles (ExportID)
); 

CREATE TABLE IF NOT EXISTS  TempExportFiles
(
	ExportId VARCHAR(32) NOT NULL,
	FileID VARCHAR(32) NOT NULL,
	ExportFilename TEXT,
	Annotation TEXT,
	CreationTime REAL,
	Excl INT DEFAULT 0,
	PRIMARY KEY (ExportId),	
	UNIQUE (ExportId)
	--FOREIGN KEY (FileID) REFERENCES Files (FileID)
); 

CREATE TABLE IF NOT EXISTS  TempXData
(
	XDataID VARCHAR(32) NOT NULL,
	XDataClass TEXT NOT NULL,
	XDataXml  TEXT NOT NULL,
	JobID  VARCHAR(32),
	Excl INT DEFAULT 0,
	PRIMARY KEY (XDataID),
	UNIQUE (XDataID)
	--FOREIGN KEY (JobID) REFERENCES Jobs (JobID)
);


CREATE TABLE  IF NOT EXISTS TempJobs
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
	processId INT,
	Excl INT DEFAULT 0,
        NewJobNumber  VARCHAR(50),
	PRIMARY KEY (JobID),
	UNIQUE (JobID),
	FOREIGN KEY (ProjectID) REFERENCES Projects (ProjectID),
	--FOREIGN KEY (ParentJobID) REFERENCES Jobs (JobID),
	--FOREIGN KEY (PreceedingJobID) REFERENCES Jobs (JobID) ON DELETE SET NULL,
	FOREIGN KEY (Status) REFERENCES JobStatus  (StatusID),
	FOREIGN KEY (Evaluation) REFERENCES JobEvaluation  (EvaluationID),
        FOREIGN KEY (UserId) REFERENCES Users (UserId),
	FOREIGN KEY (UserAgent) REFERENCES UserAgents  (UserAgentID)
) ;

CREATE TABLE IF NOT EXISTS  TempComments
(
	CommentID VARCHAR(32) NOT NULL,
	JobID VARCHAR(32) NOT NULL,
	UserID VARCHAR(255),
	TimeOfComment FLOAT NOT NULL,
	Comment TEXT,
	Excl INT DEFAULT 0,
	PRIMARY KEY (JobID, CommentID),
	UNIQUE  (CommentID)
        --FOREIGN KEY (JobID) REFERENCES Jobs (JobID)
) ;


CREATE TABLE IF NOT EXISTS  TempFileUses
(
	FileID VARCHAR(32) NOT NULL,
	JobID VARCHAR(32) NOT NULL,
	RoleID MEDIUMINT NOT NULL,
        JobParamName VARCHAR(32),
	Excl INT DEFAULT 0,
	PRIMARY KEY (FileID, JobID, RoleID)
	--FOREIGN KEY (JobID) REFERENCES Jobs (JobID),
	--FOREIGN KEY (FileID) REFERENCES Files (FileID),
	--FOREIGN KEY (RoleID) REFERENCES FileRoles (RoleID)
) ;


CREATE TABLE IF NOT EXISTS  TempJobKeyValues (
        JobID VARCHAR(32) NOT NULL,
        KeyTypeID MEDIUMINT NOT NULL,
        Value REAL,
	Excl INT DEFAULT 0,
        PRIMARY KEY (JobID,KeyTypeID)
        --FOREIGN KEY (JobId)  REFERENCES Jobs (JobId),
        --FOREIGN KEY (KeyTypeId)  REFERENCES KeyTypes (KeyTypeId) )
	);

CREATE TABLE IF NOT EXISTS  TempJobKeyCharValues (
        JobID VARCHAR(32) NOT NULL,
        KeyTypeID MEDIUMINT NOT NULL,
        Value VARCHAR(255),
	Excl INT DEFAULT 0,
        PRIMARY KEY (JobID,KeyTypeID)
        --FOREIGN KEY (JobId)  REFERENCES Jobs (JobId),
        --FOREIGN KEY (KeyTypeId)  REFERENCES KeyTypes (KeyTypeId) )
	);

CREATE TABLE IF NOT EXISTS TempFileAssociations (
       FileAssociationID VARCHAR(32) NOT NULL,
       FileAssociationTypeId MEDIUMINT NOT NULL,
       FileAssociationName VARCHAR(50),
       Excl INT DEFAULT 0,
       PRIMARY KEY (FileAssociationId),
       UNIQUE (FileAssociationId) );

CREATE TABLE IF NOT EXISTS TempFileAssociationMembers (
       FileAssociationID VARCHAR(32) NOT NULL,
       FileID VARCHAR(32) NOT NULL,
       FileAssociationRoleID MEDIUMINT NOT NULL,
       Excl INT DEFAULT 0,
       PRIMARY KEY (FileID,FileAssociationID,FileAssociationRoleID)
       --FOREIGN KEY (FileID) REFERENCES Files (FileID),
       --FOREIGN KEY (FileAssociationID) REFERENCES FileAssociations (FileAssociationID) 
       );

CREATE TABLE IF NOT EXISTS TempTags (
       TagID VARCHAR(32) NOT NULL,
       ParentTagID VARCHAR(32),
       Text VARCHAR(50) NOT NULL,
       PRIMARY KEY (TagID)
       --FOREIGN KEY (ParentTagID) REFERENCES Tags (TagID)
       );

CREATE TABLE IF NOT EXISTS TempProjectTags (
       TagID VARCHAR(32) NOT NULL,
       ProjectID VARCHAR(32) NOT NULL,
       PRIMARY KEY ( ProjectID, TagID )
       --FOREIGN KEY (TagID) REFERENCES Tags (TagID) ,
       --FOREIGN KEY (ProjectID) REFERENCES Projects (ProjectID)
       );
       
CREATE TABLE IF NOT EXISTS TempProjectComments (
	ProjectCommentID VARCHAR(32) NOT NULL,
	ProjectID VARCHAR(32) NOT NULL,
	UserID VARCHAR(255),
	TimeOfComment FLOAT NOT NULL,
	Comment TEXT,
	PRIMARY KEY (ProjectID, ProjectCommentID) ) ;
        --FOREIGN KEY (ProjectID) REFERENCES Projects (ProjectID) 

DELETE FROM TempJobs;
DELETE FROM TempFiles;
DELETE FROM TempFileUses;
DELETE FROM TempImportFiles;
DELETE FROM TempExportFiles;
DELETE FROM TempComments;
DELETE FROM TempXData;
DELETE FROM TempJobKeyValues;
DELETE FROM TempJobKeyCharValues;
DELETE FROM TempFileAssociations;
DELETE FROM TempFileAssociationMembers;
DELETE FROM TempTags;
DELETE FROM TempProjectTags;
DELETE FROM TempProjectComments;
