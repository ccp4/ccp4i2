"""
Liz Potterton Feb 2010 - Some minimal functionality to bootstrap program startup
"""

import atexit
import os
import shutil
import sys
import tempfile
import time

from PySide2 import QtWidgets

from ccp4i2 import I2_TOP


class DatabaseFailException(Exception):
    pass


def getCCP4I2Dir():
    return str(I2_TOP)


def setupEnvironment(path=''):
    os.environ["CCP4I2_TOP"] = path or getCCP4I2Dir()


def setupPythonpath():
    sys.path.insert(0, getCCP4I2Dir())


def createMissingDATABASEdbXML():
    from lxml import etree
    from core import CCP4Modules, CCP4Utils
    proj_dir_list0=CCP4Modules.PROJECTSMANAGER().db().getProjectDirectoryList()

    for proj in proj_dir_list0:
        updateDBXML = False
        projectid = proj[0]
        d = proj[2]
        dbxml = os.path.join(d,"DATABASE.db.xml")
        if not os.path.exists(dbxml):
#DATABASE.db.xml does not exist, so we create one.
            updateDBXML = True
        else:
            dbxml_mtime = os.stat(dbxml).st_mtime
#Check mtime of CCP4_JOBS and subdirs with mtime of DATABASE.db.xml
            jobsDir = os.path.join(d,"CCP4_JOBS")
            if os.path.exists(jobsDir):
                jobsDir_mtime = os.stat(jobsDir).st_mtime
                if jobsDir_mtime > dbxml_mtime:
                    print(jobsDir,"newer than",dbxml)
                    updateDBXML = True
                else:
                    jobDirs = os.listdir(jobsDir)
                    for jobDir in jobDirs:
                        thisJobDir = os.path.join(d,"CCP4_JOBS",jobDir)
                        job_mtime = os.stat(thisJobDir).st_mtime
                        if jobsDir_mtime > dbxml_mtime:
                            print(thisJobDir,"newer than",dbxml)
                            updateDBXML = True
                            break
        if updateDBXML:
            print("Rebuilding DATABASE.db.xml for project",d)
            CCP4Modules.PROJECTSMANAGER().db().exportProjectXml(projectid,fileName=dbxml)

#TODO At this point we should probably be checking that proj_dir_list0 is consistent 
#     with projectList-backup.xml
#
#     If projectList-backup.xml contains stuff not in proj_dir_list0 we probably need to abort and might need to rebuild database, or something else ...?
#     Also (not here), if DB opening fails, we might want to reconstruct it (automatically)


    dbListBackupName = os.path.join(CCP4Utils.getDotDirectory(),'projectList-backup.xml')

    if len(proj_dir_list0) > 0 and not os.path.exists(dbListBackupName):
        print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        print("Project directory list does not exist, creating")
        CCP4Modules.PROJECTSMANAGER().backupDBXML()
        print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")

    elif os.path.exists(dbListBackupName):
        utf8_parser = etree.XMLParser(encoding='utf-8')

        def parse_from_unicode(unicode_str):
            s = unicode_str.encode('utf-8')
            return etree.fromstring(s, parser=utf8_parser)
    
        with open(dbListBackupName, 'r') as infile: 
            s = infile.read()
            tree = parse_from_unicode(s)
            projs = tree.xpath("//project")
            setXML = set([str(x.text) for x in projs])
            setDB = set([str(x[2]) for x in proj_dir_list0])
            if (setXML!=setDB):
                if setXML.issubset(setDB):
                    print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
                    print("Project directory list is not up to date, recreating")
                    CCP4Modules.PROJECTSMANAGER().backupDBXML()
                    print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
                else:
                   dirListStr =  "<br/>".join([str(x) for x in (setXML - setDB)])
                   print("Backup project list is inconsistent with database contents. This might be a serious problem.")
                   res = QtWidgets.QMessageBox.warning(None,"Database file / Project list backup","Backup project list file "+dbListBackupName+" is inconsistent with database contents. This might be a serious problem.<br/><br/>The following project directories contained in "+dbListBackupName+" do not correspond to anything in the database:<br/><br/>"+dirListStr+"<br/><br/>If you think all is well, then click 'OK', otherwise click 'Abort' and investigate the problem.<br/><br/>Do not click 'OK' if you are not sure what is wrong.",QtWidgets.QMessageBox.Ok|QtWidgets.QMessageBox.Abort)
                   if res == QtWidgets.QMessageBox.Abort:
                       sys.exit()


def startBrowser(args, app=None, splash=None):
    from core import CCP4Modules
    from qtgui import CCP4WebBrowser
    from qtgui import CCP4StyleSheet, CCP4WebView
    # KJS: Removed the linux check. Unclear why it's in here.
    # if sys.platform == "linux2": win = QtWidgets.QWidget(); splash.finish(win); splash.show()
    config = loadConfig()
    config.set('graphical', True)
    config.set('qt', True)
    kw = {'graphical' : True}
    ii = 0
    while ii < len(args):
        if args[ii] in ['-db']:
            if ii+1 < len(args):
                ii += 1
                kw['dbFileName'] = os.path.abspath(args[ii])
        elif args[ii] in ['-dbLocal']:
            tFile = tempfile.NamedTemporaryFile()
            tFileName = tFile.name
            tFile.close()
            kw['dbFileName'] = tFileName
            kw['useLocalDBFile'] = tFileName
            print("Using local database file",tFileName)
            print()
            print("*********************************************************************")
            print("To avoid possible database corruption issues, you are running CCP4i2 ")
            print("in 'Local Database Mode'.")
            print()
            print("Please be aware that to avoid data loss:")
            print()
            print("You must only have one instance of CCP4i2 running at a time.")
            print("You must wait for jobs to finish before closing CCP4i2.")
            print("   If you do quit CCP4i2 whilst jobs are still running, then any")
            print("   information (output files/reports) about these jobs will be lost")
            print("*********************************************************************")
            print()
            QtWidgets.QMessageBox.information(None,"Using local copy of database","To avoid possible database corruption issues, you are running CCP4i2 in 'Local Database Mode'.<br/><br>Please be aware that to avoid data loss:<br/><br/>You must only have one instance of CCP4i2 running at a time.<br/>You must wait for jobs to finish before closing CCP4i2.<br/><br/>It is intended that these will be temporary restrictions.")
        elif args[ii] in ['-u', '-user', '-username']:
            if ii + 1 < len(args):
                ii += 1
                kw['userName'] = args[ii]
        ii += 1
    if app is None:
        app = CCP4Modules.QTAPPLICATION()
    startPrintHandler(app)
    #CCP4WebBrowser.setupWebkit()
    proj_man = startProjectsManager(**kw)
    proj_man.startCheckForFinishedJobs()
    job_cont = startJobController()
    if kw.get('dbFileName',None) is not None:
        job_cont.setDbFile(kw['dbFileName'])
    launcher = CCP4Modules.LAUNCHER()
    #print 'startup.startBrowser app',app
    createMissingDATABASEdbXML()
    def pushLocalDB():
        from core import CCP4Config, CCP4Utils
        try:
            print()
        except ValueError:
            sys.stdout = sys.__stdout__
            print("Restoring stdout")
        if 'dbFileName' in kw and 'useLocalDBFile' in kw and kw['dbFileName'] is not None and kw['useLocalDBFile'] is not None and kw['dbFileName'] == kw['useLocalDBFile']:
            print("Need to push back temporary database file ......")
            origFileName = CCP4Config.DBFILE()
            if origFileName is None:
                origFileName = os.path.join(CCP4Utils.getDotDirectory(), 'db', 'database.sqlite')
            pm = CCP4Modules.PROJECTSMANAGER()
            tFile = tempfile.NamedTemporaryFile(delete=False)
            tFileName = tFile.name
            tFile.close()
            pm.db().backupDB(tFileName)
            try:
                os.remove(origFileName+"-wal")
            except:
                pass
            try:
                os.remove(origFileName+"-shm")
            except:
                pass
            print("Copying",kw['dbFileName'],"back to",origFileName,"original DB")
            try:
                shutil.copy(tFileName,origFileName)
            except:
                time.sleep(1)
                try:
                    shutil.copy(tFileName,origFileName)
                except:
                    print("Second attempt to copy database file",kw['dbFileName'],"back to",origFileName,"failed.")
                    print("This is a serious error.")
                    QtWidgets.QMessageBox.critical(None,"Copying back database failed","Second attempt to copy database file "+kw['dbFileName']+" back to "+origFileName+" failed.<br/> This is a <b>serious</b> error. You need to type:<br/><br/>cp "+kw['dbFileName']+" "+origFileName+"<br/><br/>on the command line (terminal) to fix this.")
        else:
            print("No need to push back local db file.")
    #----------------------------------------------------------------------
    def backupListOfProjects():
        CCP4Modules.PROJECTSMANAGER().backupDBXML()
 
    #----------------------------------------------------------------------
    def closeHTTPServer(parent=None):
        # Called on exit to find and close all CHTTPServerThread threads
        from qtcore import CCP4HTTPServerThread
        if CCP4HTTPServerThread.CHTTPServerThread.insts is not None:
            CCP4HTTPServerThread.CHTTPServerThread.insts.quitServer()   # KJS : Should be ok.
    #----------------------------------------------------------------------
    proj_man.doCheckForFinishedJobs.connect(proj_man.checkForFinishedJobs)
    # the aboutToQuit() signal does not work atexit is too late to save status - 
    for func in (backupListOfProjects, pushLocalDB, proj_man.Exit, job_cont.Exit, launcher.Exit, closeHTTPServer, CCP4Modules.PRINTHANDLER().exit):
        atexit.register(func)
    # KJS : Moved the style & font changes to make up some time. Should be fine.
    CCP4StyleSheet.setStyleSheet()
    #CCP4WebView.setGlobalSettings()
    startHTTPServer(parent=app,**kw)
    CCP4WebBrowser.restoreStatus()
    CCP4WebBrowser.applyCommandLine(args)
    if hasattr(splash, "close"):
        splash.close()
    showTipsOfTheDay = CCP4Modules.PREFERENCES().SHOW_TIPS_OF_THE_DAY
    if showTipsOfTheDay:
        from qtgui.CCP4TipsOfTheDay import CTipsOfTheDay
        tipsOfTheDay = CTipsOfTheDay()
        tipsOfTheDay.exec_()


def startJobController():
    from core import CCP4Modules
    print('Starting Job Controller')
    jc = CCP4Modules.JOBCONTROLLER()
    jc.startTimer()
    print('Starting Job Controller - DONE')
    jc.restoreRunningJobs()
    return jc
 

def copyDB(origFileName,f2):
    import sqlite3
    print(origFileName,f2)
    con = sqlite3.connect(origFileName,5.0,1, check_same_thread = False)
    sql = "".join([s+"\n" for s in con.iterdump()])
    f2trunc = open(f2,"w")
    f2trunc.close()
    conbak = sqlite3.connect(f2)
    cur = conbak.cursor()
    #I am doing a simpler thing. Splitting based on semi-colon + newline is not safe. Strings may contain
    #this sequence. e.g. in comments. The simpler way may be (theoretically?) slower, but should be safer.
    try:
        cur.executescript(sql)
    except:
        print("Fail cur.executescript(sql) in startup copyDB")# ,com KJS. Not sure what this was meant to print.
        conbak.close()
        raise
    """p = re.compile(";[\r\n]+")
    sql_commands = p.split(sql)
    for com in sql_commands:
        if com == "COMMIT":
            continue
        try:
            cur.execute(com)
        except:
            print("Fail",com)
            conbak.close()
            raise"""
    conbak.commit()
    print("Saved backup database to",f2)
    conbak.close()
    con.close()


def startProjectsManager(dbFileName=None, checkForFinishedJobs=False, useLocalDBFile=None, **kw):
    from core import CCP4Config
    from core import CCP4Utils
    from core import CCP4Modules,CCP4ProjectsManager, CCP4Utils
    import string
    print(dbFileName)
    print(useLocalDBFile)
    if CCP4ProjectsManager.CProjectsManager.insts is not None:
        print('Attempting to start Project Manager for second time - should not do this!')
        return CCP4ProjectsManager.CProjectsManager.insts
    else:
        print('Starting Project Manager')
        if dbFileName is not None:
            print('Using database file: ', dbFileName)
        pm = CCP4Modules.PROJECTSMANAGER()
        try:
            if dbFileName is not None and useLocalDBFile is not None and dbFileName == useLocalDBFile:
                origFileName = CCP4Config.DBFILE()
                if origFileName is None:
                    origFileName = os.path.join(CCP4Utils.getDotDirectory(), 'db', 'database.sqlite')
                if os.path.exists(origFileName):
                    print("Copying",origFileName,"to",dbFileName,"and using this as DB")
                    t1 = time.time()
                    copyDB(origFileName,dbFileName)
                    t2 = time.time()
                    print("Copying database to local file took",t2-t1,"seconds")
                else:
                    print("Using local file",dbFileName,"as DB")

            nullHome = (CCP4Utils.getHOME() == "")
            if nullHome:
                message = "CCP4i2 can not determine your home directory's name correctly. This is probably because your username contains special (non-ASCII) characters. CCP4i2 and other CCP4 programs are unlikely to work correctly in these circumstances. Please try creating a new user without such characters for running CCP4i2. <br/><br><b>Any errors after this are likely caused by this problem.</b>"
                print(message)
                if CCP4Config.GRAPHICAL():
                    QtWidgets.QMessageBox.warning(None,"Filename warning",message)

            db = startDb(pm, mode='sqlite', fileName=dbFileName,**kw)

            if sys.version_info > (3,6):
                isAscii = db._fileName.isascii()
            else:
                isAscii = all(ord(char) < 128 for char in db._fileName)
            hasSpace = True in [c in db._fileName for c in string.whitespace]

            if not isAscii or hasSpace:
                message = 'The database filename:\n\n'+db._fileName+'\n\ncontains characters that mean ccp4i2, or other CCP4\nprogram is unlikely to work. If your UserName contains spaces or special characters and you run into problems, please try creating a new user without such characters for running CCP4i2.'
                print(message)
                if CCP4Config.GRAPHICAL():
                    QtWidgets.QMessageBox.warning(None,"Filename warning",message)

            pm.setDatabase(db)
        except DatabaseFailException:
            print("Failed to open database file!!!!!!!!")
            if dbFileName is None and CCP4Config.DBFILE() is None:
                #This is manual backup
                from qtgui import CCP4BackupDBBrowser
                dbOrigName = os.path.join(CCP4Utils.getDotDirectory(), 'db', 'database.sqlite')
                win = CCP4BackupDBBrowser.CBackupDBBrowser()
                win.setDBDirectory(os.path.join(CCP4Utils.getDotDirectory(), 'db'))
                if win.exec_():
                    f = win.selectedFile
                    print("Recovering from file",f)
                    #The journal and shm files are not going to be valid with backup, so we remove them. The backups should be 'complete'
                    try:
                        os.remove(dbOrigName + "-wal")
                    except:
                        pass
                    try:
                        os.remove(dbOrigName + "-shm")
                    except:
                        pass
                    shutil.copyfile(f, dbOrigName)
                    #FIXME - Should auto restart here.
                    QtWidgets.QMessageBox.information(None,"Restart CCP4i2","CCP4i2 will now quit. You will need to restart CCP4i2 to continue.")
                    sys.exit()
                else:
                    sys.exit()
            elif dbFileName is not None:
                res = QtWidgets.QMessageBox.warning(None,"Database file corrupted?","The specified database file "+dbFileName+" is corrupted or not a valid CCP4i2 database. As specifying the db parameter on the command line is non-standard behaviour, automatic fixing of the problem is not possible. You must specify a valid CCP4i2 database file  or seek expert assistance.")
                print(dbFileName,"corrupted")
                sys.exit()
            else:
                res = QtWidgets.QMessageBox.warning(None,"Database file corrupted?","The specified database file "+CCP4Config.DBFILE()+" is corrupted or not a valid CCP4i2 database. As specifying the dbFile parameter in the config file is non-standard behaviour, automatic fixing of the problem is not possible. You must specify a valid CCP4i2 database file in this parameter or seek expert assistance.")
                print(CCP4Config.DBFILE(),"corrupted")
                sys.exit()
        except SystemExit:
            raise
        if checkForFinishedJobs:
            pm.startCheckForFinishedJobs()
        print('Starting Project Manager - DONE')
        return pm


def startDb(parent=None, fileName=None, mode='sqlite', userName=None, userPassword=None,**kw):
    from core import CCP4Utils
    from core import CCP4ErrorHandling
    from dbapi import CCP4DbApi
    from qtgui import CCP4DbManagerGui
    try:
        db = CCP4DbApi.CDbApi(parent=parent, fileName=fileName, mode=mode, createDb=True, userName=userName,
                              userPassword=userPassword, loadDiagnostic=kw.get('loadDiagnostic',True))
    except CCP4ErrorHandling.CException as e:
        print('startDb', e[0]['code'], kw.get('graphical', False))
        if e[0]['code'] == 121 and kw.get('graphical', False):
            try:
                # Basically 121 is a username miss-match in the db.
                p_users = e[0]['details']['p_users']
                p_user = p_users[0][1] # The original username (most likely the first one listed in db)
                print('p_user', p_user)
                currentUserName = CCP4Utils.getUserId()
                dialog = CCP4DbManagerGui.CChallengeUnknownUser(currentUserName=currentUserName, origUserName=p_user)
                rv = dialog.exec_()
            except:
                rv = False
            if rv:
                ownerName = dialog.getOwner()
                reset = dialog.getReset()
                try:
                    db = CCP4DbApi.CDbApi(parent=parent, fileName=fileName, mode=mode, createDb=True,
                                          userName=ownerName, userPassword=userPassword)
                except:
                    pass
                else:
                    try:
                        db.createUser(currentUserName, userRole=CCP4DbApi.USER_ROLE_OWNER)
                    except:
                        print('Failed adding current user to database user list')
                    if reset:
                        try:
                            db.updateUser(ownerName, newUserName=currentUserName)
                        except CCP4ErrorHandling.CException as e:
                            e.warningMessage(parent=parent)
                    return db
        if kw.get('graphical', True):
            try:
                e.warningMessage(parent=parent)
            except:
                pass
            print(e.report())
            print('Please try removing $HOME/.CCP4I2/db directory so new database is created')
            if e[0]['code'] == 222:
                print('In future updates to database schema will be handled better!')
            else:
                print('Please report to ccp4@ccp4.ac.uk')
            raise DatabaseFailException
    except Exception as e:
        print(e)
        print('Please try removing $HOME/.CCP4I2/db directory so new database is created')
        print('Please report to ccp4@ccp4.ac.uk')
        raise DatabaseFailException

    return db


def loadConfig():
    from core import CCP4Config, CCP4Utils
    config = CCP4Config.CONFIG(os.path.join(CCP4Utils.getDotDirectory(), 'configs', 'ccp4i2_config.params.xml'))
    return config


def  startHTTPServer(parent=None, port=None, **kw):
    from qtcore import CCP4HTTPServerThread
    from core import CCP4Utils
    if port is None:
        port = CCP4HTTPServerThread.DEFAULT_PORT
    diry = os.path.join(CCP4Utils.getCCP4I2Dir(), 'docs')
    try:
        fileName = kw['dbFileName']
    except:
        fileName = None
    t = CCP4HTTPServerThread.CHTTPServerThread(parent=parent, parentDir=diry, port=port, fileName=fileName)
    t.start()
    return


def startPrintHandler(app):
    from core import CCP4PrintHandler
    printHandler = CCP4PrintHandler.CPrintHandler(parent=app)
    try:
        printHandler.cleanupLogs()
    except Exception as e:
        print('ERROR cleaning up print logs')
        print(e)
    # Do a getFileObject() here to ensure the filename is main_thread
    f = printHandler.getFileObject(thread=None, name='main_thread')
    sys.stdout = printHandler
    sys.stderr = printHandler
    app.aboutToQuit.connect(printHandler.exit)
    # Output version info to the print log
