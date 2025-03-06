#!/bin/sh

############### testing the script inp ####################

# use extract to run everything in example_1.inp and get a mmcif file.
extract -ext input_data/example_1.inp

# use validation-v8 to validate the mmcif file.
#validation-v8 -f script_example_1.cif -o 2 -public -exchange -adit

# move the files to some directory and delete some log files. 
mv script_example_1.cif deposit/
mv script_example_1_sf.cif deposit/
rm -f cgi_value *parser.log
#rm -f *log *err procheck* SEQUENCE.DAT *ERR validation.alignment
