import sys
import traceback
import os
import shutil
import sqlite3
import getopt

project = None
oldRoot = None
newRoot = None
dbFile = None

def usage():
    print("CCP4i2 Project Root Renamer")
    print()
    print("Usage:")
    print()
    print("python "+sys.argv[0] + " -i|--inputRoot -o|--outputRoot -f fname|--dbFile=fname [-p name|--project=name]")
    print()
    print("This program changes the project directories root in a CCP4i2 database file")
    print("and also updates any project DATABASE.db.xml files.")
    print()
    print("The main purpose is to enable the moving of an entire heirarchy of CCP4I2 projects from one")
    print("filesystem location to another, even on a different computer.")
    print()
    print(" -f fname|--dbFile=fname          : Specify path to the database file 'fname'")
    print(" -i oldRoot|--inputRoot=oldRoot   : The old projects root, e.g. '/Users/stuart/CCP4I2_PROJECTS'")
    print(" -o newRoot|--outputRoot=newRoot  : The new projects root, e.g. '/home/mcnicholas/CCP4I2_PROJECTS'")
    print(" -p name|--project=pname          : Only edit the location of project 'pname' (Optional)")
    print()
    print("e.g.:")
    print()
    print("python "+sys.argv[0] + " -i /Users/stuart/CCP4I2_PROJECTS -o /home/mcnicholas/CCP4I2_PROJECTS -f database.sqlite")
    print()
    print("will edit the database file database.sqlite, changing all ProjectDirectory")
    print("entries containing  /Users/stuart/CCP4I2_PROJECTS instead contain /home/mcnicholas/CCP4I2_PROJECTS" )
    print()
    print("python "+sys.argv[0] + " -i /Users/stuart/CCP4I2_PROJECTS -o /home/mcnicholas/CCP4I2_PROJECTS -f database.sqlite -p refine_insulin")
    print()
    print("will edit the database file database.sqlite, changing ProjectDirectory for project 'refine_insulin'")
    print("from containing  /Users/stuart/CCP4I2_PROJECTS to instead contain /home/mcnicholas/CCP4I2_PROJECTS" )

try:
    opts, args = getopt.getopt(sys.argv[1:], "i:o:p:f:", ["inputRoot=","outputRoot=","project=","dbFile="])
except  getopt.GetoptError as err:
    print(str(err))
    print()
    usage()
    sys.exit(2)

for o, a in opts:
    if o in ("-f", "--dbFile"):
        dbFile = a
    elif o in ("-i", "--inputRoot"):
        oldRoot = os.path.normpath(a)
    elif o in ("-o", "--outputRoot"):
        newRoot = os.path.normpath(a)
    elif o in ("-p", "--project"):
        project = a

if not oldRoot or not newRoot or not dbFile:
    usage()
    sys.exit(2)

conn = sqlite3.connect(dbFile)
cursor = conn.cursor()

cursor.execute("SELECT * FROM Projects ORDER BY ProjectName ASC")
rows = cursor.fetchall()

for row in rows:

    projectID   = row[0]
    projectName = row[1]

    if project is not None and project != projectName:
        continue

    newDir = row[5].replace(oldRoot,newRoot)

    if newDir == row[5]:
        print()
        print('\033[47m\033[35m' + "Warning: did not change project directory "+row[5]+"\033[0m")
        print('\033[47m\033[35m' + "because it does not contain path to be replaced '"+oldRoot+"'\033[0m")
        print()

    else:

        try:
            args = (newDir,projectID,)
            cursor.execute("UPDATE Projects SET ProjectDirectory= ? WHERE ProjectID= ?",args)
            print("Project directory for",row[1],"successfully changed to ",newDir)
            dbxml_filename = os.path.join(newDir,"DATABASE.db.xml")
            if os.path.isfile(dbxml_filename):
                dbxml_filename_new = dbxml_filename + "-new"
                with open(dbxml_filename) as f_in:
                    lines = f_in.readlines()
                    with open(dbxml_filename_new,"w+") as f_out:
                        for l in lines:
                            f_out.write(l.replace(oldRoot,newRoot))
                shutil.move(dbxml_filename_new,dbxml_filename)
                print(dbxml_filename,"successfully updated")
        except:
            print()
            print('\033[47m\033[31m' + "Error: failed to change project directory for "+row[5]+"\033[0m")
            print(traceback.format_exc())
            print()

conn.commit()
conn.close()
