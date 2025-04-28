###########################################
CCP4I2 Developers Program Command Templates
###########################################

Introduction
~~~~~~~~~~~~

CCP4i2 supports a template mechanism for generating program command
files similar to that in CCP4i. The underlying code has been re-written
as the Python module *CCP4ComTemplate* and the format of the template
file is somewhat changed to be compatible with CCP4i2 *CData* classes
and to be more Python-like. Nevertheless it should be relatively
straight-forward to convert old Tcl-based templates to the new template
style.

The template mechanism can be used from the plugins based on
*CPluginScript*. The template can either be coded within the
*CPluginScript* sub-class as the attribute COMTEMPLATE (see
`here <./pipelines.html#CPluginScript>`__) or coded in a separate file
that should be placed in the $CCP4I2/templates directory. The benefit of
the separate file is that different plugin wrappers performing slightly
different functions but using the same program could use the same
template. If the template is in a separate file then the file should be
specified in the COMTEMLATEFILE attribute of the plugin. The command
file generation is usually handled by
*CPluginScript.makeCommandAndScript* but can be run:

::

       from core.CCP4ComTemplate import CComTemplate
       comTemplate = CComTemplate(parent=self)
       # Load the template from a file
       errorReport = comTemplate.loadTemplateFromFile(fileName)
       # make a command script - pass makeComScript a CContainer containing the data
       comScript,errorReport = comTemplate.makeComScript(container) 

Alternatively if the template is in a text string:

::

       from core.CCP4ComTemplate import CComTemplate
       comTemplate = CComTemplate(parent=self,template=templateString)
       # make a command script - pass makeComScript a CContainer containing the data
       comScript,errorReport = comTemplate.makeComScript(container)

Note that methods return a
`CErrorReport <./data_classes.html#error%20handling>`__ that is
particularly useful in the case of *makeComScript()* for listing any
issues in interpreting the template.

Template Format
~~~~~~~~~~~~~~~

A one-line fragment of REFMAC template:

::

   $EXCLUDE_RESOLUTION  refi resolution $EXCLUDE_RESOLUTION_RANGE.start $EXCLUDE_RESOLUTION_RANGE.end

By default one line of the template is converted into one line in the
command file. The template line is only written to the command file if
the first word of the line, the 'write flag', evaluates to *True*. The
write flag might be a '1' (the line is always written), a reference to a
parameter or a Python expression that may include references to
parameters. In this example the first word begins with a $ which means
that it is a reference to the 'EXCLUDE_RESOLUTION' parameter. A
parameter called 'EXCLUDE_RESOLUTION' is sought in the container passed
to the *CComTemplate*. The search mechanism will search any
sub-containers in the container. If there is no 'EXCLUDE_RESOLUTION'
parameter then a report is appended to the *errorReport* and the line is
not written. If 'EXCLUDE_RESOLUTION' is found and is True then the rest
of the line is written to the command file. The words 'refi resolution'
are written automatically and then the container is searched for a
'EXCLUDE_RESOLUTION_RANGE' parameter and the 'start' and 'end'
attributes are written to the command file. If no
'EXCLUDE_RESOLUTION_RANGE' parameter is found then the whole line is not
written to the command file. Note that the names of attributes such as
'start' and 'end' must correspond to the contents on the *CData* class
that holds the data (in this case it is a *CFloatRange*). The use of
classes with attributes is a major difference of the Python CCP4i2 and
is the major issue if updating from old templates.

A more complex fragment of the REFMAC template:

::

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

In this example the first character of several lines is a hyphon (-)
which means that this line is a continuation of the previous line and
will be reassembled into one line in the command file. The second word
on continuation lines is a write flag for whether this fragment of a
line will be written to the command file. The write flag in the third
line in this fragment is a fragment of Python code (
*{['PHASE','HL'].count($INPUT_PHASE)}*) that is enclosed by curly
braces. This code fragment is testing whether the INPUT_PHASE has one of
two given values; the list *count()* will return a value of 1 if the
parameter has a value given in the list.

This example also has an IF-ELSE-ENDIF statement. The template
mechasnism supports several types of compound statement. These must be
written in upper case and on separate lines as shown in the examples
below and compound statements can be nested.

::

   LOOP N 1 $N_NONX
   1 ncsr nchains $N_NONX_CHN_TOTAL($N) chains $NONX_CHN_SRC($N) $NONX_CHN_TARG($N)
     LOOP I 1 $N_NONX_CHN($N)
        - 1 $NONX_CHN($N,$I)
     ENDLOOP
     - 1 nspans $N_NONX_SPANS_TOTAL($N)
     - 1 $NONX_INIT_SPANS_RES1($N) $NONX_INIT_SPANS_RES2($N) $NONX_INIT_SPANS_CODE($N)
     LOOP J 1 $N_NONX_SPANS($N)
       - 1 $NONX_SPANS_RES1($N,$J) $NONX_SPANS_RES2($N,$J) $NONX_SPANS_CODE($N,$J)
     ENDLOOP
   ENDLOOP

The example writes the non-crstallographic symmetry using a outer
LOOP-ENDLOOP with the variable N and two inner LOOP-ENDLOOP with the
vaiables I and J. Note that the LOOP statement does allow an *increment*
parameter after the *start* and *end* parameters. Note also how arrays
are references.

The other supported compound statement is CASE-CASEMATCH-ENDCASE. For
example:

::

   CASE $WHATEVER
   CASEMATCH 12
     1 matched 12
   CASEMATCH  24
     1 matched 24
   ENDCASE

It is possible to incorporate another template file into the present one
with the AT command.

::

   AT $CCP4I2_TOP/test/data/test_com_template_1.com

Notes
~~~~~

It is sometimes necessary to put quotes around parameters in the program
script files - for space groups is the specific example. This can be
done by appending *.quote* to the data item to be put in quotes:

::

   1 SYMM $SPACEGROUP.quote

Last modified: Thu Jun 16 17:14:24 BST 2011
