import argparse
import shutil
import sqlite3
import sys
import tempfile
import xml.etree.ElementTree as ET

from . import importDir


def parse_from_unicode(unicode_str):
    return ET.fromstring(unicode_str)


def reconstruct(inputFileName,fileName):
    with open(inputFileName, 'r') as infile: 
        dirs = []
        if inputFileName.endswith(".xml"):
            s = infile.read()
            tree = parse_from_unicode(s)
            projs = tree.xpath("//project")
            for proj in projs:
                print("Adding",proj.text.strip(),"to list of directories to process.")
                dirs.append(proj.text.strip())
        else:
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
            print("Generating database from DATABASE.db.xml files. Please wait.\nCurrently processing project directory: "+str(d))
            print("------------------------------------------------------------")
            print(d)
            importDir.importFilesFromDirXML(tfn,str(d),importProjectComments=True)

        #Now we use sqlite api to make fresh, "clean" hopefully copy of new db...
        print("Create new DB",tfn)
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
        print("Written new DB",tfn)
        conbak.close()
        
        #... and then copy to desired location.
        print("Copying new DB",tfn2,fileName)
        shutil.copy(tfn2,fileName)
        print("Copied new DB",fileName)

        print()
        print("############################################################")
        print("Saved backup database to",fileName)
        print("############################################################")

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

    reconstruct(inputFileName,fileName)
