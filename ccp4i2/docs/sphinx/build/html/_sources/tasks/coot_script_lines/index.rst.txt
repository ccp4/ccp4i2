##########################
Coot Scripting from CCP4i2
##########################

The user can provide their own script to be run non-graphically by Coot.
This is a useful way to access the model building Coot tools quickly.
The task interface has input options for the user to enter models, maps
and difference maps and (in the *Options* tab) a text edit field to
enter a coot script. When you click *Run* a Coot script is generated
that causes Coot to import and display the selected files and then has
the script you have entered and finally has a command to close down
Coot.

Here is a simple example of a script that could be entered to perform a
protein fit:

::

   fit_protein(MolHandle_0)
   write_pdb_file(MolHandle_0,os.path.join(dropDir,"refined.pdb"))

If one each of model, map and difference map have also been selected
then the final Coot script will look like this:

::

   MolHandle_0=read_pdb('/Users/lizp/Desktop/test_projects/test1/CCP4_JOBS/job_52/XYZOUT.pdb')
   MapHandle_0=make_and_draw_map('/Users/lizp/Desktop/test_projects/test1/CCP4_JOBS/job_4/FPHIOUT.mtz', 'F', 'PHI', 'PHI', 0, 0)
   DifmapHandle_0=make_and_draw_map('/Users/lizp/Desktop/test_projects/test1/CCP4_JOBS/job_4/DIFFPHIOUT.mtz', 'F', 'PHI', 'PHI', 0, 1)
   dropDir="/Users/lizp/Desktop/test_projects/test1/CCP4_JOBS/job_79/COOT_FILE_DROP"

   fit_protein(MolHandle_0)
   write_pdb_file(MolHandle_0,os.path.join(dropDir,"refined.pdb"))

   coot_real_exit(0)

| Some point note here:
| The ids for the imported models and maps have been assigned to
  variables called *MolHandle\_*\ i, *MapHandle\_*\ i and
  *DifmapHandle\_*\ i where i is 0 to the number of entered files of
  that type. You will need to substitute these variable names for the
  imported data into your script.
| CCP4i2 automatically creates a 'drop directory' for any Coot output
  files and after Coot has finished CCP4i2 will treat all files in that
  directory as output from the job and save them to the CCP4i2 database
  so they can be easily accessed by subsequent jobs. Within the script
  the drop directory is in variable name *dropDir*. Your script should
  include a command such as *write_pdb_file* to write files to the drop
  directory.
| The choice of output file name is down to you but since it will be
  seen in the GUI try to make it something meaningful.
| If you need to debug the running of your script it can be useful to
  look in the *Project directory* view (one of the tabs on the left hand
  side of the project window) - open the directory for your task (will
  be at the bottom for a recently run task) and look at the files
  *script.py* (the Coot script) and *log.txt* (the Coot output) and the
  sub-directory *COOT_FILE_DROP*.
