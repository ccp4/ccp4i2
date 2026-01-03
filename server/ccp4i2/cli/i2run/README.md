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
 --F_SIGF fileUse="SubstituteLigand[-1].F_SIGF_OUT" \
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
