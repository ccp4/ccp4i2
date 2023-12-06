##########################################
CCP4I2 Developers: Demo data for tutorials
##########################################

Introduction
~~~~~~~~~~~~

CCP4i2 has a standard mechanism to simplify user access to tutorial
documentation and data files for practical workshops using CCP4i2. The
user access is via *Utilities->Copy demo data to project*. This gives a
menu of demo data that are available as part of the CCP4i2 release. Once
selected, a set of demo data is copied to the user's project directory
for easy access. The *More info* menu item opens a web page listing a
brief description of each set of demo data and links to further
documentation. The *Download data* menu option opens a dialog box that
enables downloading data from a remote site - the user must enter the
URL of a web page. This document describes the conventions for the demo
data and how to set up a site that CCP4i2 can download new demo data
from.

The demo data are in the *$CCP4I2/demo_data* directory and more demo
data, downloaded from the web, can be placed there if the user has write
permission. Alternatively downloaded demo data is placed in the user's
'dot' directory, ($HOME/.CCP4I2 on Linux and MacOSX and CCP4I2 on
Windows).

Each set of demo data is in a single directory with a short specific
name. Besides data files the directory should contain a file called
*README.html* that describes the data files in the directory and
explains the tutorial/workshop. *README.html* can include links to other
documents - either in the directory, in sister demo data directories or
on the web. There is also an *INFO.txt* file that contains three lines:

#. a very short name for the demo data that will appear on the menu
#. a single line description that will appear in the 'More info' page
#. the sub title for a section of the 'More info' page. Please use
   existing subtitles if appropriate or something like 'York MX
   Workshop, 2016' might also be appropriate.

Templates for creating these two files can be found in
*$CCP4I2/demo_data/resources*.

Creating downloadable demo data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you want to make additional demo data sets available to CCP4i2
users,probably for use in a classroom situation, then:

#. Assemble data files in a directory with short, informative name. Give
   the files informative names. Note that your tutorial instructions
   could also include downloading data from protein database sites via
   the options on the file browser.
#. Copy the template files *$CCP4I2/demo_data/resources/README.html* and
   *$CCP4I2/demo_data/resources/INFO.txt* to your directory and edit
   them. You may want to add other tutorial files to your directory - be
   sure to link them from *README.html*.
#. Create a compressed zip file of the directory in CCP4i2 using
   *Utilities->Developers Tools->Compress Demo Data*. (Doing this via
   CCP4i2 should avoid any possible problems with zip formats.)
#. Create a directory for CCP4i2 demo downloads on your website and move
   all compressed demo files to that directory.
#. Copy the *$CCP4I2/demo_data/resources/downloadPage.html* to that
   directory and rename it *index.html*. Edit this file to create links
   to all the zip files. Be careful to retain the structure of the links
   as CCP4i2 will parse this file in order to automatically download the
   zip files.
#. The download web page can be accessed via a normal web browser but
   can be downloaded quickly to the correct place by CCP4i2: go via
   *Utilities->Copy demo data to project->Download data* to a *Download
   demo data* window and enter the URL of the download directory. Beware
   classroom PCs often have restrictive firewalls, download, write
   access etc. and to use this in a classroom situation it would be wise
   to test beforehand and have a copy of data on a USB stick!

Last modified: Wed Sep 30 18:05:45 BST 2015
