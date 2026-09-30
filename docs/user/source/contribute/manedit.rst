###############################################################
How to Contribute to the Documentation (using the command line)
###############################################################
Sometimes users may wish to use the command line for editing the documentation. Though this may be desirable for more signifcant editing of the docs, in general it is advised to use the Gitlab interface (in general this is a more user friendly way of editing, especially for minor editing).


How to use Git from the Command Line (checking out the docs)
============================================================
The first step is to obtain a cloned copy of the documentation on your local machine. This can be done with the following command::

    git clone https://gitlab.com/ccp4i2/rstdocs.git/

which will create a local copy of the documentation for editing and testing.

If you wish to update this local copy (to match what is currently in Gitlab) use::

    git pull origin master

If you wish to see a record of changes::

    git log
  
If you are unfamiliar with using Git, it is strongly recommended you familiarise yourself with Git before attempting to work with the command line. The commands given above are just the bare minimum needed to work with Git.


How to Compile the Documentation from the command line (using sphinx)
=====================================================================
In order to compile the rST files it is necessary to use sphinx which will interpret the files and convert them into html. Note this stage is automatically carried on Github, this information is only given so you can check html pages locally.

In order to use sphinx though it will be necessary to install onto your local machine. This can be done using the following commands::

    pip install -U sphinx
    pip install -U sphinx-tabs
    pip install -U sphinx-rtdtheme

After sphinx is properly installed you will be able to compile the documentation with the command::

    make html

Which should be run in the root directory of your locally checkout out clone of the documentation repository. The html files which are created by this command are located in the build directory. Please ensure that only source (rST) files, and resources used by the rST files (images etc) are checked back into the Gitlab repository - the build process is automated on Gitlab.

How to use Command Line Git to remerge the changes into Gitlab
==============================================================

After you have finished editing and checking you can commit changes you have made, back to the Gitlab repository, with the following commands (adjust these as needed, please do not cut and paste them verbatim)::

    git remote add origin https://gitlab.com/ccp4i2/rstdocs.git
    git add .
    git commit -m "Initial Commit"
    git push -u origin master

Please read the manual properly so you are aware of how Git functions, also ONLY commit rST (and image files) to gitlab - please be careful not to commit build files !
