#
#     Copyright (C) 1999-2004  Liz Potterton, Peter Briggs
#
#     This code is distributed under the terms and conditions of the
#     CCP4 Program Suite Licence Agreement as a CCP4 Library.
#     A copy of the CCP4 licence can be obtained by writing to the
#     CCP4 Secretary, Daresbury Laboratory, Warrington WA4 4AD, UK.
#
#CCP4i_cvs_Id $Id: refmac5.com,v 1.21.2.5 2009/04/23 18:28:23 mdw Exp $
########################################################################
# REFMAC5 - KEYWORD TEMPLATE
########################################################################

$EXCLUDE_RESOLUTION  refi resolution $EXCLUDE_RESOLUTION_RANGE.start $EXCLUDE_RESOLUTION_RANGE.end
1 make check $MAKE_CHECK
1 make 
 - 1 hydrogen $MAKE_HYDROGEN
 - $MAKE_HOUT  hout YES | hout NO
 - $MAKE_PEPTIDE peptide YES | peptide NO
 - $MAKE_CISPEPTIDE cispeptide YES | cispeptide NO
 - $MAKE_SSBRIDGE ssbridge YES | ssbridge NO
 - $MAKE_SYMMETRY symmetry YES | symmetry NO
 - 1 sugar $MAKE_SUGAR
 - 1 connectivity $MAKE_CONNECTIVITY
 - 1 link $MAKE_LINK
 - $MAKE_NEWLIGAND newligand noexit
 - $MAKE_REVIEW format f exit Y

IF { $REFINE_TYPE != 'REVIEW' }

######################################
# NCS RESTRAINTS
######################################
# The syntax for the NONX/NCSR keyword is
# NCSR NCHAin <nchains> \
#      CHAIns <id_1> <id_2> ... <id_nchains> \
#      NSPAin <nspans> <resno_first_1> <resno_last_1> <ncode_1> ... \
#                      <resno_first_nspans> <resno_last_nspans> <ncode_nspans>
IF { $N_NONX > 0 }
LOOP N 1 $N_NONX
1 ncsr nchains $N_NONX_CHN_TOTAL($N) chains $NONX_CHN_SRC($N) $NONX_CHN_TARG($N)
  LOOP I 1 $N_NONX_CHN($N)
     - 1 $NONX_CHN($N,$I)
  ENDLOOP
  - 1 nspans $N_NONX_SPANS_TOTAL($N)
  - 1 $NONX_INIT_SPANS_RES1($N) $NONX_INIT_SPANS_RES2($N) $NONX_INIT_SPANS_CODE($N)
  LOOP I 1 $N_NONX_SPANS($N)
    - 1 $NONX_SPANS_RES1($N,$I) $NONX_SPANS_RES2($N,$I) $NONX_SPANS_CODE($N,$I)
  ENDLOOP
ENDLOOP
ENDIF

1 refi
 - 1 type $REFINE_TYPE
 - {['PHASE','HL'].count($INPUT_PHASE)} PHASE SCBL $PHASE_SCBLUR BBLU $PHASE_BBLUR
 - 1 resi MLKF
 - 1 meth CGRAD
IF { ['IDEA','RIGID'].count($REFINE_TYPE) } 
 - 1 bref over
ELSE
 - 1 bref $B_REFINEMENT_MODE
ENDIF
 
$IFTLS refi tlsc $TLS_NCYCLES

IF { $INPUT_PHASE == 'NO'  and $TWINREF_TYPE != 'NO' }
 1  twin
ENDIF

IF {  $INPUT_PHASE == 'SAD' }
  IF { not $REF_SUBOCC  }
 1     refi oref no
  ENDIF
  { $WAVELENGTH.isSet() } anom wave $WAVELENGTH
  LOOP N 0 $ANOMALOUS_FORM_FACTORS.__len__()
 1    { $ANOMALOUS_FORM_FACTORS[$N].isSet() } anom form $ANOMALOUS_FORM_FACTORS[$N].atomName $ANOMALOUS_FORM_FACTORS[$N].Fp $ANOMALOUS_FORM_FACTORS[$N].Fpp
  ENDLOOP
ENDIF

IF { $REFINE_TYPE != "RIGID"  } 
  1 ncyc $NCYCLES
ELSE
  1 rigid ncycle $RIGID_NCYCLES
  LOOP N 1 len($RIGID_DOMAINS)
    1 rigid group $N
    LOOP  J 1 {len($RIGID_DOMAINS[$N].selection)}
      - 1 from $RIGID_DOMAINS[$N].selection[$J].firstRes $RIGID_DOMAINS[$N].selection[$J].chainId to $RIGID_DOMAINS[$N].selection[$J].lastRes $RIGID_DOMAINS[$N].selection[$J].chainId
    ENDLOOP
    - { $RIGID_DOMAINS[$N].excludeSideChains } excl schain
    - { $RIGID_DOMAINS[$N].applyRotation and $RIGID_DOMAINS[$N].transformation.rotation.isSet() } euler $RIGID_DOMAINS[$N].transformation.alpha $RIGID_DOMAINS[$N].transformation.beta  $RIGID_DOMAINS[$N].transformation.gamma
    - { $RIGID_DOMAINS[$N].applyTranslation and  $RIGID_DOMAINS[$N].transformation.translation.isSet() } trans $RIGID_DOMAINS[$N].transformation.x $RIGID_DOMAINS[$N].transformation.y $RIGID_DOMAINS[$N].transformation.z
  ENDLOOP
ENDIF

$BLIM blim $BLIM_MIN $BLIM_MAX

1 scal 
 - { [StringSame $BULK_SOLVENT_SCALING BULK ] } type BULK | type SIMP
IF { [StringSame $BULK_SOLVENT_SCALING BULK ] }
 - { [IfSet $BULK_SCALING_RESOLUTION_MIN] && [IfSet $BULK_SCALING_RESOLUTION_MAX] } reso  $BULK_SCALING_RESOLUTION_MAX $BULK_SCALING_RESOLUTION_MIN 
ELSE
 - { [IfSet $SIMPLE_SCALING_RESOLUTION_MIN] && [IfSet $SIMPLE_SCALING_RESOLUTION_MAX] } reso $SIMPLE_SCALING_RESOLUTION_MAX $SIMPLE_SCALING_RESOLUTION_MIN
ENDIF
 - 1 LSSC
 - $SCALING_ANISO ANISO
 - $SCALING_IF_FIXB FIXB bvalue $SCALING_FIXB_BBULK
 - $SCALING_EXPE_SIGMA EXPE
 - { [StringSame $SCALING_REF_SET FREE ] } FREE
 - { [StringSame $SCALING_APPLY NO ] } appl $SCALING_APPLY

IF { [ StringSame $MLSC_REF_SET WORK ] }
1 scale mlsc WORK
ENDIF

{ $EXCLUDE_FREER && ![StringSame $EXCLUDE_FREER_VALUE 0 ] } free $EXCLUDE_FREER_VALUE

######################################
# SOLVENT MASK
######################################
$IF_SOLVENT solvent YES | solvent NO
- { [IfSet $SOLVENT_VDWPROB] } VDWProb $SOLVENT_VDWPROB
- { [IfSet $SOLVENT_IONPROB] } IONProb $SOLVENT_IONPROB
- { [IfSet $SOLVENT_RSHRINK] } RSHRink $SOLVENT_RSHRINK

# Always use MATRIX weighting
1 weight
 - { !$EXPERIMENTAL_WEIGHTING } NOEX
IF { $AUTO_WEIGHTING }
 - 1 AUTO
ELSE
 - 1 MATRIX $MATRIX_WEIGHT
ENDIF

1 monitor $MONI_LEVEL
  - 1 torsion $MONI_TORSION
  - 1 distance $MONI_DISTANCE
  - 1 angle $MONI_ANGLE
  - 1 plane $MONI_PLANE
  - 1 vanderwaals $MONI_VANDERWAALS
  - 1 chiral $MONI_CHIRAL
  - 1 bfactor $MONI_BFACTOR
  - 1 bsphere $MONI_BSPHERE
  - 1 rbond $MONI_RBOND
  - 1 ncsr $MONI_NCSR

IF { ![StringSame $REFINE_TYPE IDEA ] }  

LABELLINE labin $LABIN
 - { [StringSame $INPUT_PHASE PHASE ] } FOM= $FOMB
IF $EXCLUDE_FREER
  LABELLINE - FREE
ENDIF

LABELLINE labout $LABOUT

ENDIF

# This combines the B factor weighting use of TEMP with the TLS initialisation
# use of TEMP
IF { $IFTMP && $IFTLS && $IFBFAC_SET }
  $IFTMP temp $WBSKAL $SIGB1 $SIGB2 $SIGB3 $SIGB4 set $BFAC_SET
ELSE
  $IFTMP temp $WBSKAL $SIGB1 $SIGB2 $SIGB3 $SIGB4 
  IF {$IFTLS}
    $IFBFAC_SET temp set $BFAC_SET
  ENDIF
ENDIF
$IFADDU tlso addu

$IFDIST dist $WDSKAL $SIGD1 $SIGD2 $SIGD3 $SIGD4 $SIGD5
$IFANGL angle $ANGLE_SCALE
$IFPLAN plan $WPSKAL $SIGPP $SIGPA
$IFCHIR chir $WCSKAL $SIGC
$IFNCSR ncsr $WSSKAL $SIGSP1 $SIGSP2 $SIGSP3 $SIGSB1 $SIGSB2 $SIGSB3
$IFTORS torsion $WTSKAL $SIGT1 $SIGT2 $SIGT3 $SIGT4
$IFVAND vand 
- { [IfSet $WVSKAL ] } overall $WVSKAL
- { [IfSet $WAND_SIGMA_VDW ] }  sigma vdw $WAND_SIGMA_VDW
- { [IfSet $WAND_SIGMA_HBOND ] } sigma hbond $WAND_SIGMA_HBOND
- { [IfSet $WAND_SIGMA_METAL ] }  sigma metal $WAND_SIGMA_METAL
- { [IfSet $WAND_SIGMA_TORS ] } sigma tors $WAND_SIGMA_TORS
- { [IfSet $WAND_INCR_TORS ] } incr tors $WAND_INCR_TORS
- { [IfSet $WAND_INCR_ADHB ] } incr adhb $WAND_INCR_ADHB
- { [IfSet $WAND_INCR_AHHB ] } incr ahhb $WAND_INCR_AHHB
IF { $IFISO && [ StringSame $B_REFINEMENT_MODE anis ] }
{ [IfSet $SPHERICITY ] } sphericity $SPHERICITY
{ [IfSet $RBOND ] } rbond $RBOND
ENDIF 

AT { [FileJoin [GetEnvPath CCP4I_TOP] templates harvest.com ] }

ENDIF

# Append an external script file specified by developer
$USE_INCLUDEFILE # External script file:
$USE_INCLUDEFILE $INCLUDEFILE_CMD

1 END
