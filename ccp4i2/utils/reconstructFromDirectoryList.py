import argparse
import os
import shutil
import sqlite3
import sys
import tempfile
import xml.etree.ElementTree as ET

import reconstructDBFromXML
if __name__ == "__main__":
    sys.path.append(os.path.join(os.path.dirname(__file__),".."))
import importDir

if __name__ == "__main__":

    parser = argparse.ArgumentParser(usage="A script to convert a list of projectdirectories in a file to CCP4I2 format project exprot files (.ccp4_export)")
    parser.add_argument('--inputFileName', help='Filename listing input projects',required=True)
    parser.add_argument('--dbFile', help='Output database filename',required=True)
    args = parser.parse_args()
    if args.dbFile is not None and args.inputFileName is not None:
        fileName = args.dbFile
        inputFileName = args.inputFileName
    else:
        #argparse should stop us ever getting here...
        print(parser.usage)
        sys.exit()

    with open(inputFileName, 'r') as infile: 
        dirs = []
        for line in infile:
            dirs.append(line.strip())
            print("Adding",line.strip(),"to list of directories to process.")

        tfile = tempfile.NamedTemporaryFile(delete=False)
        tfn = tfile.name+".sqlite"
        tfile2 = tempfile.NamedTemporaryFile(delete=False)
        tfn2 = tfile2.name
        tfile.close()
        tfile2.close()

        for d in dirs:
            print("Generating database from project files. Please wait.\nCurrently processing project directory: "+str(d))
            print("------------------------------------------------------------")
            print(d)
            project_tree = reconstructDBFromXML.generate_xml_from_project_directory(str(d))
            if len(project_tree.xpath("//ccp4i2_body/jobTable")) == 0:
                continue
            if len(project_tree.xpath("//ccp4i2_body/jobTable/job")) == 0:
                continue
            ET.indent(project_tree)
            outl = ET.tostring(project_tree).decode()
            dbxmlout = os.path.join(str(d),"DATABASE.db.xml")
            with open(dbxmlout,"w+") as outfd:
                outfd.write(outl)

            importDir.importFilesFromDirXML(tfn,str(d))

        #Now we use sqlite api to make fresh, "clean" hopefully copy of new db...
        con = sqlite3.connect(tfn)
        sql = "".join([s+"\n" for s in con.iterdump()])

        f2trunc = open(tfn2,"w")
        f2trunc.close()
        print("Writing to",tfn2)

        conbak = sqlite3.connect(tfn2)
        cur = conbak.cursor()

        try:
            cur.executescript(sql)
        except:
            print("Fail",com)
            conbak.close()
            raise

        conbak.commit()
        
        #... and then copy to desited location.
        shutil.copy(tfn2,fileName)

        print()
        print("############################################################")
        print("Saved backup database to",fileName)
        print("############################################################")
        conbak.close()

