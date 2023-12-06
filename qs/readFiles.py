from __future__ import print_function

import sys
import sqlite3

def readFiles(f,projId):

    com, args = "SELECT Jobs.TaskName,Jobs.JobNumber, Files.FileId, Files.Annotation,ImportFiles.SourceFilename,Files.FiletypeID,Files.FileSubType,Files.JobParamName, Jobs.CreationTime FROM Jobs INNER JOIN Files ON Files.JobId=Jobs.JobId LEFT OUTER JOIN ImportFiles ON Files.FileId=ImportFiles.FileId WHERE Jobs.ProjectId=?  AND ImportFiles.ImportId IS NULL AND Jobs.ParentJobId IS NULL ORDER BY Jobs.CreationTime DESC;", [projId]

    conn = sqlite3.connect(f)
    cursor = conn.cursor()
    cursor.execute(com,args)
    files = cursor.fetchall()
    conn.close()

    print("----------------------------")
    print("Got",len(files),"files")
    print("----------------------------")
    return files

if __name__ == "__main__":
    f = sys.argv[1]
    files = readFiles(f,'8cb7c294c25711e4803360f81db2874e')

    print(len(files))
    for f in files:
        print(f)

