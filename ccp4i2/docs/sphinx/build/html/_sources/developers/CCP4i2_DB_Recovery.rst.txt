Database Recovery Guide
=======================

This document describes four ways that the CCP4i2 database can be
recovered.

#. This method assumes that nothing is wrong with actual databse file
   and that there is some sqlite/OS problem.
#. This method assumes that the working DB is broken but that meta-data
   XML files can be used to reconstruct the database.
#. This method assumes that the working DB is broken and depends on
   CCP4i2 recognising this and that it can solve it itself
   semi-automatically.
#. This method is for when working DB is broken, there is a working
   backup, but CCP4i2 cannot solve it automatically.
#. This method is when there is no working DB or backups. *This method
   is also useful or transferring many projects from one computer to
   another without having to export/import them all.*

Method 1 - Reboot
-----------------

Shut down all instances of CCP4i2 and then reboot computer.

Method 2 - Reconstruct from XML
-------------------------------

**NB. This method requires that you have been running a very recent
version of CCP4i2, which writes project metadata and project list to
disk periodiccaly.**

-  Windows

   ::

       set CCP4I2=%CCP4%\share\ccp4i2
       ccp4-python %CCP4I2%\utils\reconstructFromDBXML.py

-  Linux/Mac (csh)

   ::

       setenv CCP4I2 $CCP4/share/ccp4i2
       ccp4-python utils/reconstructFromDBXML.py --inputFileName $HOME/projectList-backup.xml --dbFile ccp4i2_backup_db.sqlite

-  Linux/Mac (sh)

   ::

       CCP4I2=$CCP4/share/ccp4i2 ccp4-python utils/reconstructFromDBXML.py --inputFileName $HOME/projectList-backup.xml --dbFile ccp4i2_backup_db.sqlite

A file will be written in your home directory called

::

    ccp4i2_backup_db.sqlite

This file is a new database and needs to be put in the location of your
"working" database.

Start ccp4i2, then quit i2 to recreate dummy empty ``~/.CCP4I2`` (if one
does not already exist).

Then move the reconstructed database into place.

-  Windows
   Move ``$HOME/ccp4i2_backup_db.sqlite`` to
   ``$HOME/CCP4I2/db/database.sqlite`` either with a command prompt or
   Explorer window.
-  Linux/Mac

   ::

       mv $HOME/ccp4i2_backup_db.sqlite $HOME/.CCP4I2/db/database.sqlite

Remove the wal, shm files from the dummy db.

::

    rm $HOME/.CCP4I2/db/database.sqlite-shm $HOME/.CCP4I2/db/database.sqlite-wal





    Method 3 - Database Backup GUI

    If a GUI pops-up saying that database is corrupted, simply follow the instructions to pick a backup file (usually the
    most recent) and restart CCP4i2.


    Method 3 - Copying Database Backup Files by hand
    Following applies to Linux/Mac, for Windows, the appropriate directory is %USERPROFILE%\CCP4I2\db 
    or %HOMEDRIVE%\%HOMEPATH%\CCP4I2\db.

    cd ~/.CCPI2/db
    ls -l 

    You should see something like this:
    fasolt:db stuart$ ls -l
    total 26824
    -rw-r--r--  1 stuart  staff   536576  8 Nov 09:02 database.sqlite
    -rw-r--r--  1 stuart  staff    32768  5 Nov 14:16 database.sqlite-shm
    -rw-r--r--  1 stuart  staff  1722192  8 Nov 09:19 database.sqlite-wal
    -rw-------  1 stuart  staff   327680  5 Nov 14:18 database_sqlite_backup-05112018-141812
    -rw-------  1 stuart  staff   327680  5 Nov 15:07 database_sqlite_backup-05112018-150733
    -rw-------  1 stuart  staff   327680  5 Nov 15:26 database_sqlite_backup-05112018-152607
    -rw-------  1 stuart  staff   327680  5 Nov 15:27 database_sqlite_backup-05112018-152703
    -rw-------  1 stuart  staff   327680  5 Nov 15:27 database_sqlite_backup-05112018-152745
    -rw-------  1 stuart  staff   344064  5 Nov 15:30 database_sqlite_backup-05112018-153016
    -rw-------  1 stuart  staff   348160  5 Nov 15:31 database_sqlite_backup-05112018-153152
    -rw-------  1 stuart  staff   348160  7 Nov 11:27 database_sqlite_backup-07112018-112718
    -rw-------  1 stuart  staff   348160  7 Nov 11:28 database_sqlite_backup-07112018-112829
    -rw-------  1 stuart  staff   352256  7 Nov 11:33 database_sqlite_backup-07112018-113332
    -rw-------  1 stuart  staff   356352  7 Nov 11:34 database_sqlite_backup-07112018-113438
    -rw-------  1 stuart  staff   356352  7 Nov 11:41 database_sqlite_backup-07112018-114118
    -rw-------  1 stuart  staff   364544  7 Nov 11:44 database_sqlite_backup-07112018-114423
    -rw-------  1 stuart  staff   389120  7 Nov 14:18 database_sqlite_backup-07112018-141833
    -rw-------  1 stuart  staff   389120  7 Nov 14:19 database_sqlite_backup-07112018-141923
    -rw-------  1 stuart  staff   397312  7 Nov 14:19 database_sqlite_backup-07112018-141949
    -rw-------  1 stuart  staff   413696  7 Nov 14:21 database_sqlite_backup-07112018-142114
    -rw-------  1 stuart  staff   430080  7 Nov 14:24 database_sqlite_backup-07112018-142401
    -rw-------  1 stuart  staff   430080  8 Nov 08:55 database_sqlite_backup-08112018-085500
    -rw-------  1 stuart  staff   430080  8 Nov 08:57 database_sqlite_backup-08112018-085709
    -rw-------  1 stuart  staff   438272  8 Nov 08:57 database_sqlite_backup-08112018-085737
    -rw-------  1 stuart  staff   442368  8 Nov 09:00 database_sqlite_backup-08112018-090001
    -rw-------  1 stuart  staff   450560  8 Nov 09:03 database_sqlite_backup-08112018-090339
    -rw-------  1 stuart  staff   450560  8 Nov 09:06 database_sqlite_backup-08112018-090652
    -rw-------  1 stuart  staff   454656  8 Nov 09:07 database_sqlite_backup-08112018-090737
    -rw-------  1 stuart  staff   454656  8 Nov 09:08 database_sqlite_backup-08112018-090828
    -rw-------  1 stuart  staff   327680 17 Oct 16:11 database_sqlite_backup-17102018-161157
    -rw-------  1 stuart  staff   327680 17 Oct 16:12 database_sqlite_backup-17102018-161246
    -rw-------  1 stuart  staff   327680 17 Oct 16:15 database_sqlite_backup-17102018-161520

| **If there are lots of backup files like above**, then you can copy
  one of these to your "main" database file:
| e.g.

::

    cp database_sqlite_backup-17102018-161520 database.sqlite
    rm database.sqlite-shm database.sqlite-wal

| The second command gets rid of temporary files which would be invalid
  with the backup file.
| Restart ccp4i2.

Method 4 - Reconstruct GUI
--------------------------

-  Windows

   ::

       set CCP4I2=%CCP4%\share\ccp4i2
       ccp4-python %CCP4I2%\utils\reconstructBrowser.py

-  Linux/Mac

   ::

       sh $CCP4/share/ccp4i2/bin/reconstructDB

Click OK. After a while (depending on number of projects/jobs/disk
speed, etc.) a file will be written in your home directory called

::

    ccp4i2_backup_db.sqlite

This file is a new database and needs to be put in the location of your
"working" database.

Start ccp4i2, then quit i2 to recreate dummy empty ``~/.CCP4I2`` (if one
does not already exist).

Then move the reconstructed database into place.

-  Windows
   Move ``$HOME/ccp4i2_backup_db.sqlite`` to
   ``$HOME/CCP4I2/db/database.sqlite`` either with a command prompt or
   Explorer window.
-  Linux/Mac

   ::

       mv $HOME/ccp4i2_backup_db.sqlite $HOME/.CCP4I2/db/database.sqlite

Remove the wal, shm files from the dummy db.

::

    rm $HOME/.CCP4I2/db/database.sqlite-shm $HOME/.CCP4I2/db/database.sqlite-wal


