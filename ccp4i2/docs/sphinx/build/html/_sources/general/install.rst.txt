###############################
CCP4i2 Installation for Testers
###############################

The present mechanism for installing CCP4i2 relies on having a recent
CCP4mg installed and the libraries in CCP4mg being used by i2. The
Python and other (mostly) uncompiled stuff is accessed via the source
control at CCPForge. The advantage of this is that it should be possible
to get updates quickly and easily. The disadvantage is that you are at
the bleeding edge and may see some (very temporary) bugs.

Using CCPForge to Access CCP4i2 Code
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You need to register with CCPForge and join the CCP4i2 project. To do
that:

-  Go to `ccpforge.cse.rl.ac.uk <http://ccpforge.cse.rl.ac.uk>`__
-  Click **Register New Account** (currently top right of home page).
-  Fill in form, submit and wait for email response
-  Once you are registered login in and go to the CCP4i2 page.
   (Currently this is click on the **Projects** tab at the top of the
   CCP4Forge home page and then scroll down and go onto the **Next
   page** and you should see CCP4i2 listed). On the CCP4i2 page there is
   a **Request to Join Project** button on the right hand side.
-  After your request I will have to reply by email to confirm OK so
   there may be some delay

To access the SVN you will need the environment variables set:

::

   SVN_EDITOR=vi
   svnroot=http://ccpforge.cse.rl.ac.uk/svn/ccp4i2

Then the command

::

    svn co $svnroot

should do the trick. After this you should update your installation by
going into the ccp4i2 directory and typing the command:

::

   svn update

If you get any odd behaviour with this try typing before retrying

::

   svn cleanup

If the update stops to warn you of conflicts then accept 'their' version
(that will be my updated version). The ccpforge site is usually pretty
reliable but may occasionally be down.

| The checked out project is all Python or uncompiled files (so no
  building required!). You also need to have installed a recent nightly
  build of CCP4mg - this is used by i2 to display results but also
  provides the Qt and other libraries for CCP4i2. This can be got from:
| http://www.ysbl.york.ac.uk/~ccp4mg/nightly/
| If you are running on Linux you will need to edit the file
  *ccp4i2/bin/browser* to set the path to CCP4MG (search the file for
  'CCP4MG'). If you are running on MacOSX and have installed CCP4mg in
  the default /Applications this should not be necessary. You will also
  need to setup an up-to-date version of the CCP4 suite before starting
  CCP4i2 with the command:

::

   ccp4i2/bin/browser

For Developers Only
~~~~~~~~~~~~~~~~~~~

Developers may find the file *ccp4i2/utils/setup.sh* useful to set up a
developer environment. This will need some intelligent editing for your
system. I source this file in my .bashrc.

Other Programs
~~~~~~~~~~~~~~

CCP4i2 currently may launch CCP4mg,Coot or viewhkl. These can be
downloaded from:

+---------+---------+------------------------------------------------+
| Coot    | Linux   | http://lmb.b                                   |
|         |         | ioch.ox.ac.uk/coot/software/binaries/releases/ |
+---------+---------+------------------------------------------------+
|         | Mac     | http://sage.ucsc.edu                           |
|         |         | /~wgscott/xtal/wiki/index.php/Stand-Alone_Coot |
+---------+---------+------------------------------------------------+
|         | Windows | http://                                        |
|         |         | www.ysbl.york.ac.uk/~lohkamp/coot/wincoot.html |
+---------+---------+------------------------------------------------+
| hklview | All     |                                                |
+---------+---------+------------------------------------------------+

| CCP4i2 needs to know the path of the executable for each package. This
  is set in the config file which is (in Linux speak):
| $HOME/.CCP4I2/configs/ccp4i2_config.params.xml

This is an XMl file with a **searchPath** element with sub-elements that
look like:

::

       
         <ccp4mg>
            <macosx>
              <path>/Applications/QtMG.app/Contents/ccp4mg*/bin/ccp4mg </path>
              <path>/usr/local/ccp4mg*/bin/ccp4mg </path>
            </macosx>
            <linux>
               ....
            <windows>
               .....
          </ccp4mg>

For each program and each OS there is one or more **path** elements.
These may include '*' characters to wildcard features such as version
numbers. I2 will test these in order take the first file that exisits.
If a wildcard results in two or more hits then the one with the latest
modification date is used. You may need to edit this file (though should
have default install paths set for Mac and Windows).

If your config file does not contain the **searchPath** element delete
the file and rerun CCP4i2 browser.

--------------

Last modified: Fri Dec 14 12:56:41 GMT 2012
