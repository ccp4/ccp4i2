###################################
CCP4I2 Generating Reports using XRT
###################################

**This document describes using the XRT mechanism to generate reports -
the**\ `Python api <./reports.html>`__\ **is now the preferred method.**

The report generator is a work-in-progress - please let Liz know about
any features that you need or any holes in the (acknowledged to be
sketchy) documentation.

Introduction
~~~~~~~~~~~~

The report generator tool (CCP4ReportGenerator.py) expects two input
files:

#. Either an **xrt** report template file or a Python script which
   defines the layout of the final report
#. One or more xml data files output by either a script or a program
   which provide the data that is substituted into the report. These
   files must be valid xml with unique paths to the required data but
   otherwise there is no restriction on the format.

From these the report generator will create an html report file (in
principle it could create reports in alternative formats). But note that
the report may contain Qt widgets that can only be displayed in the
CCP4i2 browser. The report generator is run automatically within the
CCP4i2 interface when a job completes. Some extra features such as
input/output data and job status are added to the report automatically
so need not be specified explicitly in the report template. The same
template mechanism can be used to produce 'real time' reports while the
job is running.

| A task developer needs to write the template file as:
| $CCP4I2/tasks/*mytask*/*mytask*.xrt or
| $CCP4I2/tasks/*mytask*/*mytask*.running.xrt (for a 'real-time' updated
  report)

Presently the best example xrt is in the buccaneer_build_refine task. To
understand the xrt file you need to be familiar with XPATH (see the
`w3schools tutorial <http://www.w3schools.com/xpath/default.asp>`__)
which is used to reference the data in the XML files. A common attribute
for elements of the xrt is 'select' whose value should be an XPATH to an
element in the program xml. The text of this element is then substituted
into the report. An alternative source of data is the CCP4i2 database
which contains the jobs input/output files and other data that is listed
below in `Appendix A <appendix%20a>`__. This data is accessed using a
'database' attribute whose value should be an xpath-style reference to
the data outlined in the appendix. As with the select attribute the
value from the database is substituted into the report.

XRT files are constructed from two types of tags: container tags which
may contain any other container or content tags, and content tags which
may only contain specified elements. The top level tag of every XRT file
is the <report> container.

Container tags: <report>, <fold>, <if>, <ifnot> and <loop>
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

<report>
^^^^^^^^

This is the top level tag of a report document. It may contain other
container tags and content tags.

<fold>
^^^^^^

This tag allows a detail section to be added whose content is hidden
until a 'more detail' button is clicked. It may contain other container
tags (including other folds) and content tags.

.. table:: Attributes

   +-----------+----------------+---------------------------------------------------+
   | **label** | *Show details* | Label to appear on fold describing fold contents. |
   +-----------+----------------+---------------------------------------------------+

<if>
^^^^

Only act on the contents of the 'if' container if the 'select' element
in the program xml contains the text 'True'.

.. table:: Attributes

   ========== ========================================================
   **select**  Path of xml data element which can be 'True' or 'False'.
   ========== ========================================================

See the `w3schools
tutorial <http://www.w3schools.com/xpath/default.asp>`__ to explain path
syntax.

<ifnot>
^^^^^^^

Only act on the contents of the 'if' block if the 'select' element in
the program xml contains the text 'False'.

.. table:: Attributes

   ========== ========================================================
   **select**  Path of xml data element which can be 'True' or 'False'.
   ========== ========================================================

See the `w3schools
tutorial <http://www.w3schools.com/xpath/default.asp>`__ to explain path
syntax.

<loop>
^^^^^^

Apply the same processing (i.e. the xrt within the loop block) to all
elements found by applying the 'select' criteria to the data file.

.. table:: Attributes

   ========== ==========================
   **select**  Path of xml data elements.
   ========== ==========================

The files $CCP4I2/test/data/test_report_loops.xrt and
$CCP4I2/test/data/test_report_loops.xml demonstrate this and can be run
in the *ccp4py* interpreter by:

::

       prompt> from report.CCP4ReportParser import test
       prompt> test('test/data/test_report_loops.xrt','test/data/test_report_loops.xml')
     

The resulting file, *report.html* is placed in your current working
directory.

Content tags: <text>, <copy>, <table>, <graph>, <status>, <picture>, <launcher>, <help>
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

<text>
^^^^^^

The text tag is used to add text to a report. The text may come from one
of two sources - it may be included in the body of the tag (e.g.
<text>Job title:<text>), or it may be extracted from the task XML file
by specifying an XPATH in the select attribute (e.g. <text
select="//job"/>).

.. table:: Attributes

   +--------------+---+-------------------------------------------------+
   | **if**       |   | Name of parameter in the XML which must be True |
   |              |   | for text to be output                           |
   +--------------+---+-------------------------------------------------+
   | **select**   |   | Path of xml data element from which to source   |
   |              |   | the text                                        |
   +--------------+---+-------------------------------------------------+
   | **tag.text** |   | Text explicitly specified in the XRT            |
   +--------------+---+-------------------------------------------------+

<copy>
^^^^^^

Copy block of XML directly from the data file. The block should be valid
html.

.. table:: Attributes

   +------------+---+---------------------------------------------------+
   | **select** |   | Path of xml data element from which to source the |
   |            |   | text to be copied                                 |
   +------------+---+---------------------------------------------------+

See the `w3schools
tutorial <http://www.w3schools.com/xpath/default.asp>`__ to explain path
syntax.

<table>
^^^^^^^

The table tag inserts a table into the report. The table layout and
headers are specified in the XRT report file, the values are extracted
from the task XML.The table element has <data> sub-elements which
specify the column(/row) data.The table element can take two attributes,
transpose changes the orientation of the table, and select specifies an
XPATH expression which specifies the elements containing row data to be
selected. Both are optional. Normally the select element will be
specified and returns a list of XML elements from the task XML, each of
which contains all the information for a single table row. Occasionally
the data in the task XML may not be in a suitable form for this, in
which case the attribute is omitted and the < <data> select attributes
used on their own. The risk of using this approach is that there is now
no guarantee that the columns are the same length.<table> <data> Table
data elements are used to specify column(/row) data. The data element
can take four attributes: title, subtitle, expr and select. The title
attribute gives the column header. The subtitle attribute allows tables
with nested headers. The select attribute specifies an xpath with
selects the data to be displayed. If a row element has been selected in
the <table> tag, then this should specify a unique element within that
row element. If no row element has been selected in the <table> tag, the
select should return a list of elements which will become the column
values. In the latter form there is no guarantee columns will be of the
same length.The expr tag allow a python expression to be applied to the
value before display.See example below. <data> with no select but
containing a list of <value> elements should be implemented to allow
lists of row labels, for example.

.. table:: Attributes

   ============= =================================================
   **select**     Path to XML file data
   **transpose**  If this attribute is set then table is horizontal
   ============= =================================================

.. table:: **Data Element Attributes**

   +--------------+---+-------------------------------------------------+
   | **title**    |   | Text to appear in table                         |
   +--------------+---+-------------------------------------------------+
   | **subtitle** |   | Text to appear in table                         |
   +--------------+---+-------------------------------------------------+
   | **select**   |   | Path to data item releative to the table select |
   |              |   | attribute value                                 |
   +--------------+---+-------------------------------------------------+
   | **expr**     |   | ??                                              |
   +--------------+---+-------------------------------------------------+

This example selects data from a refmac task XML file and produces a
table as shown below.

::

   <table select="//Overall_stats/stats_vs_cycle/new_cycle">
   <data title="Cycle"    select="cycle"    expr="int(x)-1" />
   <data title="R-factor" select="r_factor" />
   <data title="R-free"   select="r_free"   />
   <data title="RMS deviations" subtitle="Bond"   select="rmsBOND" />
   <data   subtitle="Angle"  select="rmsANGLE" />
   <data> subtitle="Chiral" select="rmsCHIRAL" /></table>

<graph>
^^^^^^^

This is an interface to the **pimple** graph display program which will
be shown inline within the report. *pimple* can display multiple plots
in one window with two combo-boxes to allow selection of the current
plot. The left-hand combo-box will select a data table and the
right-hand box will select between views on that table. The *graph*
block described here can only specify one data table with, optionally,
multiple views specified by the *plot* block. Multiple *graph* blocks
can be grouped in a *graphgroup* block to be displayed in the same
*pimple* instance with the *title* attributes from each graph apperaring
in the left-hand combo-box.

The *graph* element can have a *select* attribute that gives the xpath
for a base element in the xml from which subsequent xpath selections are
made.

The *graph* element has four possible sub-elements:

<title>
'''''''

The graph title which will appear in the left-hand combo-box in *pimple*

<data>
''''''

These elements work the same as the *data* sub-elements of *table* to
extract data that is dispersed throughout the xml file.

.. _table-1:

<table>
'''''''

An alternative means to extract a table from the xml file if the data is
present in a format acceptable for the *CCP4Table* element in *pimple*
input (i.e. as rows of space separated data values). This element should
have a *select* attribute that gives the xpath of an xml element
containing a table.

<plot>
''''''

This block contains the plotting directives for *pimple*. The *plot*
block of the xrt is copied straight to the output for processing by
*pimple* except for any plot elements with *select* attributes: for
these the text from the element in the data xml specified by the search
xpath is copied into the text of the *plot* element and the select
attribute is deleted.

It is also possible for the top-level *plot* element to have a *select*
attribute which gives the xpath to a block of plot directives in the
xml. In this case the plot directives will be taken from the xml; but,
note, this abuses the principle that scientific data comes from the xml
and presentation directives are in the xrt.

There can be multiple *plot* blocks within one *graph* element - this
implies that multiple views of the same data table are possible and are
listed in the right-hand combo-box in *pimple*.

<picture>
^^^^^^^^^

Will insert a picture in the report with buttons to launch either CCP4mg
or Coot with the same view. The picture is generated by CCP4mg which
runs in background without appearing on the computer screen. CCP4mg
needs a *scene.xml* file to specify what it is to display and save as a
picture. The scene file specifies which data files (either coordinates
or maps) are to be loaded and how they should be drawn. The drawing
specification can be either a *simple scene* or generated using the
CCP4mg Picture Wizard. Many other aspects of the display can also be
specified using CCP4mg's usual XML. The template for the scene.xml file
can be incorporated directly iton the xrt file or (probably better)
provided in a separate file that is incorporated using the *filename*
attribute of the `<insertXrt> <insertXrt>`__ tag. Typically some
parameters such as filenames need to be substituted into the scene.xml
template and this should be done using *database* or *select* attributes
to specify the source of the substituted data.

A simple example is in the *buccaneer.xrt* where a picture element uses
an *xrtInsert* element to specify a separate xrt file which contains the
template for a CCP4mg scene file:

::

   <xrt:picture label="Final structure" >
   <xrtInsert filename=$CCP4I2/tasks/buccaneer/buccaneer.scene.xrt/>
   </xrt:picture>

The contents of the CCP4mg scene file are explained in `Appendix
B <appendix%20b>`__.

.. table:: Attributes

   +----------+---+-----------------------------------------------------+
   | template |   | Name of CCP4mg picture definition file template     |
   |          |   | (see notes above)                                   |
   +----------+---+-----------------------------------------------------+
   | label    |   | label to appear by picture                          |
   +----------+---+-----------------------------------------------------+

<launch>
^^^^^^^^

Insert a button in the report to launch another viewer. Currently
supported CCP4mg and Coot. The initial view in CCP4mg can be set by
specifying a *pictureDefinitionFile* attribute (see notes on <picture>
tag).

.. table:: Attributes

   +---------------------------+---+-------------------------------+
   | **exe**                   |   | Name of program to launch     |
   |                           |   | (ccp4mg/coot)                 |
   +---------------------------+---+-------------------------------+
   | **label**                 |   | Label to appear on launch     |
   |                           |   | button                        |
   +---------------------------+---+-------------------------------+
   | **pictureDefinitionFile** |   | Only relevant for ccp4mg -    |
   |                           |   | specify picture definition    |
   |                           |   | file name                     |
   +---------------------------+---+-------------------------------+
   | **ccp4_data_id**          |   |                               |
   +---------------------------+---+-------------------------------+

<help>
^^^^^^

This puts a *Help* link by tables, graphs etc. and is highly
recommended! The *ref* attribute gives the link reference which can be a
file path starting with an environment variable (eg $CCP4I2) and the
appropriate full path will be substituted into the html.

.. table:: Attributes

   ===== =======================
   label  Text for help link
   ref    Link reference for help
   ===== =======================

<inputData>, <outputData> and <jobDetails
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

These are usually set automatically by the report generator and need not
be in the XRT file.

Other tags
~~~~~~~~~~

<insertXrt>
^^^^^^^^^^^

This provides a mechanism to copy fragments of xrt code generated by the
program or script into the xrt template. Or to insert standard fragments
of from any specified file. The copying is done in the report generator
as part of a pre-processor step. The *select* attribute should be used
to specify an element in the program xml file that contains one or more
xrt elements that are copied into the xrt replacing the *insertXrt*
element. Note that the 'top' level element specified by *select* is not
copied.

Alternatively the *filename* attribute should give the full path name of
a file to be inserted. This path may begin with the 'environment
variable' $CCP4I2 to specify the ccp4i2 directory.

.. table:: Attributes

   +--------------+---+-------------------------------------------------+
   | **select**   |   | Path of element in program output to be         |
   |              |   | inserted in the report template specification   |
   +--------------+---+-------------------------------------------------+
   | **filename** |   | Full path name for an XML file whose contents   |
   |              |   | will be incorporated into the xrt script.       |
   +--------------+---+-------------------------------------------------+

See the `w3schools
tutorial <http://www.w3schools.com/xpath/default.asp>`__ to explain path
syntax. Note that the xrt in the program output xml must specify
namespaces so, for example, to insert a plot specification the program
xml will be

::

   <plotDefinition>
   <plot xmlns:xrt="http://www.ccp4.ac.uk/xrt">
   <title>Completeness by cycle
   ....
   </plot>
   </plotDefinition>

and the xrt file will contain..

::

   <xrt:insertXrt select="plotDefinition"/>

<comment>
^^^^^^^^^

The enclosed comment text will not be interpreted into the report.

<drillDown>
^^^^^^^^^^^

Handling for drill down to reports from sub-jobs is pending. Probably
will be handled automatically.

Testing A Report Template
~~~~~~~~~~~~~~~~~~~~~~~~~

There is a simple script $CCP4I2/bin/makeReport that can be used to test
a template. The console command is:

::

   prompt> makeReport xrtFileName xmlFileName [-j jobId] [-o outputFileName] [-p]

The jobId allows the report generator to access the database for extra
info without which some of the features will be missing. The jobId is
not the same as the job number seen on the gui. A test value can be
found as described below. *-p* list the resultant file to console.

Alternatively the report generator can be run from the `i2 Python
console <./introduction.html#console>`__ with the following:

::

   >>> from report import CCP4ReportParser
   >>> CCP4ReportParser.test(xrtfile,xmldatafile,jobId)

This will list the report to the terminal and write it to *report.html*
in the running directory. The *jobId* parameter is used to retrieve
information from the database. The best way to get this is to run a job
of the appropriate task from the gui and then look in the meta data at
the top of the params.xml file.

Appendix A - The Database data accessible from an XRT file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

   filenames / XYZIN
             / HKLIN
             / XYZOUT
             / HKLOUT
             etc
   runtime
   status
   taskname
   jobid
   jobnumber
   projectid
   projectname
   tasktitle
   fileroot
   descendentjobs

Appendix B - The CCP4mg scene file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

NB this will only work with recent nightly build versions of CCP4mg
(post 16th Oct 2012).

A simple example of scene template:

::

   <scene>
     <data>
       <MolData id='id1'>
          <filename database="filenames/XYZOUT"/>
       </MolData>
     </data>
     <wizard><template>ribbons:colour chains</template>
       <parameters>
         <MolData1>id1</MolData1>
       </parameters>
     </wizard>
   </scene>

This template requires only one peice of substituted data - the *XYZOUT*
file from the *database* is substituted as the filename for the
*<MolData><filename>* element to give the scene file:

::


     <ccp4i2_header>
       <function>MGSCENE</function>
       <projectName>s12</projectName>
       <userId>lizp</userId>
       <ccp4iVersion>0.0.1</ccp4iVersion>
       <jobId>66b40dd7146111e2875a3c0754185dfb</jobId>
       <projectId>969cdfae138111e2acd83c0754185dfb</projectId>
       <creationTime>14:39 12/Oct/12</creationTime>
       <jobNumber>4</jobNumber>
     </ccp4i2_header>
     <ccp4i2_body>
       <scene>
         <data>
          <MolData id="id1">
            <filename>/Users/lizp/Desktop/test_projects/s12/CCP4_JOBS/job_4/XYZOUT.pdb</filename>
          </MolData>
         </data>
         <wizard>
           <template>ribbons:colour chains</template>
           <parameters>
             <MolData1>id1</MolData1>
           </parameters>
         </wizard>
      </scene>
     </ccp4i2_body>
   </ccp4:ccp4i2>

This file has the standard CCP4i2 xml structure with a header containing
metadata describing the provenance of the file. The body of the file is
one or more *<scene>* elements. Each *<scene>* element contains a <data>
block that specifies the files to be loaded and a *<wizard>* block that
specifies how the data should be drawn using a CCP4mg wizard template.

The <data> block should list the files to be loaded as *<MolData>* or
*<MapData>* elements. These elements must have unique id attributes and
a *<filename>* sub-element. Usually the actually filename will need to
be substituted in from the filenames in the CCP4i2 database so the
*database* attribute should be used to specify the source of the
filename.

The *<wizard>* element should contain a *<template>* sub-element that
specifies the CCP4mg wizard template in the form *folder:template name*.
A template can be selected from the existing wizard templates which can
be viewed in CCP4mg or something appropriate can be created (see
ccp4mg/help/picture_wizard.html and talk to Liz or Stuart). Picture
wizard templates are partially a Python script which can apply some
intelligence in setting up a display but they can also have some options
which can be set. The options are in the *CHOICES* section of the wizard
template file. The *<wizard>* sub-element *<parameters>* should contain
sub-elements specifying values for the wizard 'choices' and particularly
should provide values for *MolData* and *MapData* parameters which
cross-refernce to the *id* attributes of the sub-elements of the
*<data>* section. Note that although the wizard is mostly used in CCP4mg
to setup the display of a single coordinate file it can also be used to
specify the display of map objects and of multiple data objects.

There is an alternative mechanism to specify a simple scene that does
not use the wizard, for example:

::

   <scene>
     <data>
       <MolData id='id1'>
          <filename database="filenames/XYZOUT"/>
          <MolDisp>
           <select>peptide</select>
           <colour>bychain</colour>
           <style>CYLINDERS</style>
          </MolDisp>
       </MolData>
       <MapData id='id2'>
         <filename database="filenames/HKLOUT"/>
         <columns> <F>2FOFCWT </F> <PHI>PH2FOFCWT </PHI> </columns>
         <difColumns> <F>FOFCWT </F> <PHI>PHFOFCWT </PHI> </difColumns>
         <model>id1 </model>
         <gridSize>0.5 </gridSize>
         <contourUnits>sigma </contourUnits>
         <MapDisp>
           <contourLevel>3 </contourLevel>
           <difference>1 </difference>
         </MapDisp>
        </MapData>
     </data>
   </scene>

The *MolData* element can have one or more *MolDisp* sub-elements. Each
*MolDisp* element should contain *select*, *colour* and *style* elements
and can contain much more such as control of labels, visibility,
flashing - see the CCP4mg documentation. (????). Note that a *select*
attribute could be added to any of the elements so their values can be
taken from the program output.

Any of the aspects of CCP4mg presentation that can be controlled using
the CCP4mg XML format can also be incorporated in the scene file.
Probably the most useful feature is control of the view using:

::

      <View>
         <centre_MolData>id2 </centre_MolData>
         <centre_selection>A/1170(DUP) </centre_selection>
         <orientation_auto> <molData>id2 </molData> <selection>A/1170(DUP) </selection> </orientation_auto>
         <scale_auto>1 </scale_auto>
      </View>

Here the *MolData* object and (optionally) an atom selection that will
be the centre of the view are specified. The *orientation_auto* element
then tells CCP4mg to try to find the best orientation to show a specied
atom selection. *scale_auto* is a flag to CCP4mg to scale to zoom of the
specified atom selelction rather than show all visible atoms.

Alternatively the scale,orientation and centering can be specified
explicitly. In this example variously other display parameters are set -
these could also be used the the 'auto' view parameters shown above.

::

         
      <View>
        <scale>58.0125
        <orientation>
          <q0>0.0134139428983 </q0>
          <q1>0.336496398576 </q1>
          <q2>-0.507403789037 </q2>
          <q3>-0.793178186004/q3>
       </orientation>
       <centre_xyz>
         <x>-23.676</x>
         <y>8.538</y>
         <z>5.262</z>
       </centre_xyz>
       <fog_enabled>1</fog_enabled>
       <fog_near>-1.1</fog_near>
       <fog_far>26.01</fog_far>
       <slab_offset>27.5</slab_offset>
       <slab_enabled>false</slab_enabled>
       <slab_width>55.0</slab_width>
     </View>

Note that it is possible to have more than one *scene* specified in the
file. It is intended that CCP4mg will use these either to create
multiple pictures in one program run or, if it is opened in interactive
mode, to create a Presentation which gives the user a choice of views.

Last modified: Fri Feb 15 10:20:15 GMT 2013
