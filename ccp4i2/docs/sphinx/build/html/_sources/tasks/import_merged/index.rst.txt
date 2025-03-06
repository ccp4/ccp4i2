#############
Import merged
#############

   You will usually start a project with unmerged data or even images
   but if you have already merged data then this task will import it
   into CCP4I2 - converting to the MTZ format if necessary and creating
   a Reflection data object and also a FreeR data object if FreeR data
   is available. The task will then run some of the tools usually
   associated with Data Reduction to provide a report on the quality of
   the data.

Input
=====

   The input to this task requires spacegroup, cell and wavelength - if
   this information is present in the file it will be read and displayed
   - otherwise you will have to enter the data. MTZ files may contain
   multiple sets of reflection data and you will be required to select
   the set you want to use.
   The task will import files in the usual CCP4 MTZ format, mmCIF format
   (as is usually distributed by the PDB (PDBe, RCSB, PDBj) with a file
   name of the form **xnyz\ sf.ent**) and Scalepack format (usually has
   extension *.sca*). The latter two formats are readable files whose
   first few lines are of the form:

Merged reflection data file formats
-----------------------------------

MMCIF
^^^^^

   ::

      data_r4oogsf
      # 
      _audit.revision_id    1_0
      _audit.creation_date  2014-04-16
      _audit.update_record   'Initial release' 
      # 
      _cell.entry_id      4oog 
      _cell.length_a      157.9720 
      _cell.length_b      183.8040 
      _cell.length_c      61.2890 
      ...
      ...
      loop_
      _refln.crystal_id
      _refln.wavelength_id
      _refln.scale_group_code
      _refln.index_h
      _refln.index_k
      _refln.index_l
      _refln.status
      _refln.F_meas_au
      _refln.F_meas_sigma_au
      1 1 1    0    0    2 o  108.900   1.6000 
      1 1 1    0    0    4 o   73.900   1.1000 
      1 1 1    0    0    6 o   32.700   0.5000 
      1 1 1    0    0   10 o   19.100   0.5000 
      1 1 1    0    0   12 o   56.300   0.9000 

   MMCIF files are very flexible and very variable, beware the same
   format may contain atomic coordinates rather than reflection data.
   The file may or may not contain space group, cell and wavelength
   information.

   .. rubric:: Scalepack
      :name: scalepack

   ::

          1
       -987
         108.742    61.679    71.652    90.000    97.151    90.000 c2        
       -50   0   2  3530.5   733.8
       -50   0   3   210.8   357.2
       -50   0   6  1930.8   547.3

   This file format contains spacegroup and cell data but not
   wavelength.

   .. rubric:: Other formats
      :name: other-formats

   CCP4I2 supports the commonly used reflection file formats but there
   are some others which can be converted to MTZ by the appropriate CCP4
   program. Unmerged data files should be imported via the Data Reduction task

   ================== ==================
   Program:           Format:
   convert2mtz        CNS, PHS and SHELX
   dtrek2mtz          d*trek
   ================== ==================

   Results

   .. container:: tab-content
      :name: tab-content2

      The report from this job is a sub-set of the report from the Data
      reduction task, 
      see `documentation here <../aimless_pipe/aimless_pipe.html>`__
