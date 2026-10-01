# i2run invocations

Two entry points, both of which parse the same arguments:

    # anywhere a pip-installed ccp4i2 is importable (what the job panel's
    # "i2run command" button renders)
    ccp4-python -m ccp4i2.cli.i2run <task> --project_name <proj> [--PARAM value ...]

    # from a dev checkout's server/ directory
    ccp4-python manage.py i2run <task> --project_name <proj> [--PARAM value ...]

The `i2run` console script declared in pyproject.toml is shadowed on a normal
CCP4 setup by the legacy Qt `$CCP4/bin/i2run`, so prefer `-m`.

## Naming a file

A file parameter can be given a path (`fullPath=...`), a database id
(`dbFileId=...`), or -- usually best -- a reference to where the file came from:

    --XYZIN  "fileOut=[3].XYZOUT"                  XYZOUT of job 3
    --F_SIGF "fileIn=[2].F_SIGF"                   what job 2 used as F_SIGF
    --XYZIN  "fileOut=prosmart_refmac[-1].XYZOUT"  the most recent refmac's XYZOUT
    --DICT_LIST "fileIn=[2].DICT_LIST[0]"          one element of a file list

`fileOut=` is a file a job produced; `fileIn=` is one it consumed. Saying which
matters because four tasks carry the same file parameter name in both
`inputData` and `outputData`. A negative index counts back in creation order and
skips jobs that have not got the file, so `[-1]` is the most recent job that
actually produced one; an unqualified positive index is a job *number*.

Negative indices must be bracketed. `--XYZIN "-1.XYZOUT"` cannot work: argparse
treats a leading `-` as an option.

`fileUse=` is a deprecated alias that names no direction (produced first, then
consumed). It still resolves, with a warning.

A task-qualified non-negative index is a literal position among that task's
jobs, unrun ones included: `aimless_pipe[1]` is the second aimless job, which
may be a clone that has not run. To name one particular job, use its number,
`fileOut=[1].HKLOUT[0]`.

A file inside a list item is named the same way, as a sub-key of the item:

    --ENSEMBLES label=MDMX use=True number=1 \
                pdbItemList/identity_to_target=0.57 \
                "pdbItemList/structure/fileOut=chainsaw[-1].XYZOUT"

It resolves exactly as a top-level reference does, so the job uses that
job's file and records it as used. (`.../fullPath=` would import the file again
as a new one of no known origin.)

Worked parameter examples follow.

ccp4-python manage.py i2run prosmart_refmac --project_name refmac_gamma_test_0
case1 = """aimless_pipe \
 --UNMERGEDFILES \
 crystalName=hg7 \
 dataset=DS1 \
 file=$CCP4I2_TOP/demo_data/mdm2/mdm2_unmerged.mtz \
 --project_name refmac_gamma_test_0"""

case2a = """aimless_pipe \
 --UNMERGEDFILES \
 crystalName=hg7 \
 dataset=DS1 \
 file=$CCP4I2_TOP/demo_data/mdm2/mdm2_unmerged.mtz \
    --XYZIN_REF fullPath=$CCP4I2_TOP/demo_data/mdm2/4hg7.pdb \
 --MODE MATCH \
 --REFERENCE_DATASET XYZ \
 --project_name refmac_gamma_test_0"""

case2b = """prosmart_refmac \
 --F_SIGF fileOut="SubstituteLigand[-1].F_SIGF_OUT" \
 --XYZIN \
 fullPath=$CCP4I2_TOP/demo_data/mdm2/4hg7.pdb \
        selection/text="not (HOH)" \
    --prosmartProtein.REFERENCE_MODELS \
        fullPath=$CCP4I2_TOP/demo_data/mdm2/4qo4.cif \
 --project_name SubstituteLigand_test_0"""

case3 = """phaser_simple \
 --F_SIGF \
 fullPath=$CCP4I2_TOP/demo_data/beta_blip/beta_blip_P3221.mtz \
 columnLabels="/_/_/[Fobs,Sigma]" \
 --F_OR_I F \
 --XYZIN \
 $CCP4I2_TOP/demo_data/beta_blip/beta.pdb \
 --project_name refmac_gamma_test_0"""
