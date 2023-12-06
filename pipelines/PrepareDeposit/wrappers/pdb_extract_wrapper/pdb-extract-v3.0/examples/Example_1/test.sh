#!/bin/sh

############### testing command line ####################
# use pdb_extract to extract the required statistics and get a mmcif file.
pdb_extract  -e  MAD \
-s HKL -ilog input_data/sclepack1.log  \
-p CNS -iLOG  input_data/mad_sdb.dat input_data/mad_summary.dat input_data/mad_fp.dat \
-d CNS -iLOG input_data/density_modify.dat  \
-r CNS -iCIF input_data/deposit_cns.mmcif \
-iENT input_data/data_template.text \
-o Example_1.cif 


# use pdb_extract_sf to convert the structure factor to mmcif format.
pdb_extract_sf  -rt F -rp CNS -idat input_data/gere-nat.cv  \
-dt I -dp HKL -c 1 -w 1 -idat input_data/w1.sca  \
-c 1 -w 2 -idat input_data/w2.sca  \
-c 1 -w 3 -idat input_data/w3.sca -o Example_1.sf.cif

# use validation-v8 to validate the mmcif file.
#validation-v8 -f Example_1.cif -o 2 -public -exchange -adit


# move the files to some directory and delete some log files. 
mv Example_1.cif deposit
mv Example_1.sf.cif deposit
#mv UNNAMED.lett validation_result/Example_1.lett
#mv UNNAMED.ps   validation_result/Example_1.ps
#mv procheck* SEQUENCE.DAT *ERR validation.alignment UNNAMED* validation_result/
#rm -f *log *err 

