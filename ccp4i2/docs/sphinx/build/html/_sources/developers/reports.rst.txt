CCP4I2 Generating Reports
=========================

Introduction
~~~~~~~~~~~~

The layout and content of job reports for CCP4i2 containing text,
graphs, tables, pictures and more can be defined in Python. This
supersedes an earlier version of report generation using an XML report
template (*xrt*) that is described `here <./reports_xrt.html>`__. The
*CCP4ReportParser* module provides classes to represent all the features
of the report and simple functions to specify the report in minimal
code. For each task there should be a sub-class of
*CCP4ReportParser.Report* which is the top-level of the report
definition that I will refer to as the Python report definition - *prd*.

Usually the program or script will output an xml file containing the
job-specific data for the report. When the job completes (or otherwise -
see below) ccp4i2 instantiates the *prd* object and passes it a
representation of the XML data file. The *Report* class should construct
a tree of *Text*, *Graph*, *Table* etc. objects that represent the
required report and substitute data from the XML file (and other sources
if necessary) into the definition. The *Report.as_html_file()* method is
then called by the ccp4i2 framework to save the report to an HTML file.
Note that the report may contain Qt widgets that can only be displayed
in the CCP4i2 browser. Some extra features such as input/output data and
job status are added to the report automatically so need not be
specified explicitly in the report code. In future the functionality
could be extended to create alternative file formats. Note that Appendix
C contains a style guide to current best practice.

The *prd* file can be placed in either the *wrapper*, *pipeline* or
*task* directory for the task. At startup ccp4i2 will check for files
containing sub-classes of *Report* and make a lookup table based on the
*TASKNAME* parameter in the *prd*. If there is an error in the Python
file which prevents it being imported when CCP4i2 starts then this is
reported to the terminal and the *Program print log* (under *Utilities*
pull-down menu). Application developers should not normally need to use
it but the function *TASKMANAGER().getReportClass()* will return the
*Report* sub-class for a given taskname. Within ccp4i2 the *Report*
sub-class is instantiated with two arguments that may determine the sort
of output required: *jobStatus* which currently may have values
'Running', 'Finished', 'Unsatisfactory', or 'Failed' and *reportMode*
which currently may have values None (the default) or 'subjob' which
implies the request for a report is for a pipeline sub-job. Currently
the 'subjob' request will come from a user clicking a sub-job in the
Project Viewer and picking *View->Job report* from the context menu. The
*prd* should be customised to support these options.

| The report class is should create the output appropriate for the
  *jobStatus* which may be:
| **Finished** - the most common case - should create the usual report.
| **Running** - the job is still running and probably this requires
  generating a short report showing progress (e.g. a graph showing
  RFactor v. cycle). If this case is supported then a class variable
  *RUNNING* should be set *True*
| **Unsatisfactory** - the job has completed but without a good
  scientific result. The report may be similar to the 'Finished' report
  but should say what has happened and suggest what to do about it.
| **Failed** - the job has failed; a generic failed job report is
  usually created by the core system but if the *CReport* sub-class
  variable *FAILED* is set True then the a *CReport* sub-class is also
  created and put at the top of the generic report. In this case a short
  report on what has failed and what to do about it is appropriate.

Handling XML
~~~~~~~~~~~~

The CCP4i2 libraries include `lxml <lxml.de>`__ which handles XML within
Python and provides many tools. For users convenience lxml provides the
same api as an earlier *etree* package so you are more likely to see the
word 'etree' rather than 'lxml' in code. The internal representation of
XML is as an `etree.Element <lxml.de/tutorial.html>`__ for which there
is brief (but likely enough for ccp4i2 developers) documentation within
the Python documentation
`here <http://docs.python.org/2/library/xml.etree.elementtree.html#module-xml.etree.ElementTree>`__
(scroll down for Element). The most useful tool is
`xpath <http://www.w3schools.com/xpath>`__ which queries the XML.

In order to tidy up the report definition code there is a small CCP4i2
extension to *etree.Element* called *CCP4ReportParser.SearchableElement*
with the following methods (where input *path* is an xpath expression -
see the `xpath <http://www.w3schools.com/xpath>`__ docs). :

+----------+--------------------+------------------------------------+
| select   | path               | Return the content (element.text   |
|          |                    | or attribute) of the first item    |
|          |                    | satisfying the expression. Return  |
|          |                    | empty string if nothing satisfies  |
+----------+--------------------+------------------------------------+
| ifSelect | path,default=False | Return a Boolean from first item   |
|          |                    | satisfying the expresion or return |
|          |                    | the default                        |
+----------+--------------------+------------------------------------+
| haspath  | path               | Return True/False if there is an   |
|          |                    | item satisfying the expression     |
+----------+--------------------+------------------------------------+
| xpath0   | path               | Return the first element           |
|          |                    | satisfying the expression or       |
|          |                    | return None                        |
+----------+--------------------+------------------------------------+
| xpath1   | path               | Return the first element           |
|          |                    | satisfying the expression or       |
|          |                    | return an empty element (tag       |
|          |                    | ='dummy')                          |
+----------+--------------------+------------------------------------+

To use these extensions the XML file must be loaded using a customised
parser. The 'standard' job XML output is passed to the Report object
already appropriately parsed but if you wish to load any other XML files
use the method *Report.loadXmlFile(fileName)*

CCP4ReportParser classes
~~~~~~~~~~~~~~~~~~~~~~~~

The report definition will have a tree structure with a *Report* object
as the root node. The child objects may be either container objects that
have children or content objects that represent features such as text or
tables.\ *Report* is a sub-class of *Container* (which should not be
used directly). The other useful sub-classes are *Fold* which creates a
'fold' in the report whose contents can be displayed or hidden by
clicking on a title line and *GraphGroup* or *PictureGroup* which groups
together multiple graphs or pictures to be displayed in the same widget.
There are several other *Container* types that are used to support the
logic in the XRT report definitions but these are not necessary to
support Python definitions.

The *\__init_\_* method of all classes have three arguments that pass
along some basic data:

======= ================= =============================================
xrtnode etree.Element     Supports XRT file method of report definition
xmlnode etree.Element     The job XML output
jobInfo Python dictionary Job information from the database
======= ================= =============================================

There are also optional arguments to specify the *id* or *class*
attribute for the html element that may be useful for customising the
appearance of the report. Note that *class* is a reserved word in Python
so the Python argument is called *class\_* and also some report features
have a fixed *class* attribute needed to support the standard
functionality.

The Python report definition should usually begin:

::

   from report.CCP4ReportParser import *
   from core.CCP4ErrorHandling import *

   class my_task_report(Report):
     # Specify which gui task and/or pluginscript this applies to
     TASKNAME = 'my_task'
     def __init__(self,xmlnode=None,jobInfo={},**kw):
       Report. __init__(self,xmlnode=xmlnode,jobInfo=jobInfo,**kw)
       if self.errorReport().maxSeverity()>SEVERITY_WARNING:
         print 'FAILED instantiating my_task report generator'
         self.errorReport().report()
         return
       graph = self.addGraph( title="The best results" , select="//results/best" )

The task-specific report object will be instantiated with appropriate
*xmlnode* and *jobInfo* data by the ccp4i2 framework and must pass these
arguments onto the *Report* sub-class which will hand them on to its
child objects. The *xrtnode* argument is not used in the Python api and
can be ignored.

After calling *Report.__init__()* it is a good idea to check for errors
- *Container.errorReport()* will return an instance of
`CCP4ErrorHandling.CErrorReport <./error_handling.html>`__ which is a
list of errors recorded in the *Report* object and any child objects. If
there are any errors with severity greater than warning (likely to be
failing to read a file) then it is best to stop and, to aid debugging,
you can print out the error report with the *CErrorReport.report()*
method.

It is a good idea to split the report definition into different methods
for each element of the report rather than put the full definition in
the \__init__() method. This will make it easier to change the report
structure or provide customised reports in different context (e.g. a
'Running' report may contain just a graph showing progress per cycle).
See the
`buccaneer_build_refine <../../tasks/buccaneer_build_refine/buccaneer_build_refine_report.py>`__
for an example of this.

The application developer should usually create a new object using
methods such as *addGraph()* which creates a new object passing on the
standard arguments, adds the new object to its list of children and
returns the new object.

The table below lists the *addWhatever* methods and their arguments in
addition to the three standard and two optional (*id* and *class\_*)
arguments. The methods of the individual content classes are also
included.

+----------------+----------------+----------------+----------------+
| **Container**  | -              | -              | Base container |
|                |                |                | class - should |
|                |                |                | not be used    |
|                |                |                | itself         |
+----------------+----------------+----------------+----------------+
| \__init_\_     | -              | -              | -              |
+----------------+----------------+----------------+----------------+
| -              | xrtnode        | etree.Element  | Supports XRT   |
|                |                |                | file method of |
|                |                |                | report         |
|                |                |                | definition     |
+----------------+----------------+----------------+----------------+
| -              | xmlnode        | etree.Element  | The job XML    |
|                |                |                | output         |
+----------------+----------------+----------------+----------------+
| -              | jobInfo        | Python         | Job            |
|                |                | dictionary     | information    |
|                |                |                | from the       |
|                |                |                | database       |
+----------------+----------------+----------------+----------------+
| loadXmlFile    | -              | -              | Read and parse |
|                |                |                | an XML file to |
|                |                |                | return an      |
|                |                |                | etree.Element  |
+----------------+----------------+----------------+----------------+
| -              | fileName       | string         | File path      |
+----------------+----------------+----------------+----------------+
| addDiv         | -              | -              | Return a new   |
|                |                |                | *Div* object   |
+----------------+----------------+----------------+----------------+
| addFold        | -              | -              | Return a new   |
|                |                |                | *Fold* object  |
+----------------+----------------+----------------+----------------+
| -              | label          | string         | Label to       |
|                |                |                | appear on fold |
|                |                |                | bar            |
+----------------+----------------+----------------+----------------+
| -              | initiallyOpen  | Boolean        | Determines     |
|                |                |                | whether the    |
|                |                |                | *Fold* object  |
|                |                |                | is initially   |
|                |                |                | open           |
+----------------+----------------+----------------+----------------+
| -              | brief          | string         | A short label  |
|                |                |                | to appear on   |
|                |                |                | top bar of     |
|                |                |                | report display |
+----------------+----------------+----------------+----------------+
| addGraphGroup  | -              | -              | Return a new   |
|                |                |                | *GraphGroup*   |
|                |                |                | object         |
+----------------+----------------+----------------+----------------+
| -              | launcher       | string         | Do not display |
|                |                |                | the graphs in  |
|                |                |                | the report but |
|                |                |                | display a      |
|                |                |                | button with    |
|                |                |                | this label     |
|                |                |                | that will      |
|                |                |                | launch the     |
|                |                |                | graph viewer   |
|                |                |                | displaying the |
|                |                |                | graphs         |
+----------------+----------------+----------------+----------------+
| a              | -              | -              | Return a new   |
| ddPictureGroup |                |                | *PictureGroup* |
|                |                |                | object (not    |
|                |                |                | fully          |
|                |                |                | implemented)   |
+----------------+----------------+----------------+----------------+
| addText        | -              | -              | Add text to    |
|                |                |                | report         |
+----------------+----------------+----------------+----------------+
| -              | text           | string         | text to add to |
|                |                |                | report - <p>   |
|                |                |                | tags will be   |
|                |                |                | added if text  |
|                |                |                | is not         |
|                |                |                | well-formed    |
|                |                |                | xml            |
+----------------+----------------+----------------+----------------+
| -              | select         | string         | xpath string   |
|                |                |                | to select text |
|                |                |                | from the xml   |
|                |                |                | input          |
+----------------+----------------+----------------+----------------+
| addTable       | -              | -              | Add a table to |
|                |                |                | the report     |
+----------------+----------------+----------------+----------------+
| -              | select         | string         | xpath          |
|                |                |                | expression to  |
|                |                |                | select list on |
|                |                |                | items in xml   |
|                |                |                | input          |
+----------------+----------------+----------------+----------------+
| -              | exclud         | Boolean        | If true        |
|                | eIfDataMissing |                | exclude any    |
|                |                |                | column from    |
|                |                |                | table for      |
|                |                |                | which data not |
|                |                |                | found in xml   |
+----------------+----------------+----------------+----------------+
| -              | transpose      | Boolean        | If true print  |
|                |                |                | out table      |
|                |                |                | transposed     |
|                |                |                | with row       |
|                |                |                | headers        |
+----------------+----------------+----------------+----------------+
| -              | help           | string         | -              |
+----------------+----------------+----------------+----------------+
| addGraph       | -              | -              | Add a graph to |
|                |                |                | the report or  |
|                |                |                | graphgroup     |
+----------------+----------------+----------------+----------------+
| -              | select         | string         | xpath          |
|                |                |                | expression to  |
|                |                |                | select list on |
|                |                |                | items in xml   |
|                |                |                | input          |
+----------------+----------------+----------------+----------------+
| -              | title          | string         | Title for the  |
|                |                |                | graph          |
+----------------+----------------+----------------+----------------+
| -              | help           | string         | -              |
+----------------+----------------+----------------+----------------+
| -              | launcher       | string         | Do not display |
|                |                |                | the graph in   |
|                |                |                | the report but |
|                |                |                | display a      |
|                |                |                | button with    |
|                |                |                | this label     |
|                |                |                | that will      |
|                |                |                | launch the     |
|                |                |                | graph viewer   |
|                |                |                | displaying the |
|                |                |                | graph          |
+----------------+----------------+----------------+----------------+
| addPicture     | -              | -              | Add a picture  |
|                |                |                | to the report  |
+----------------+----------------+----------------+----------------+
| -              | scene          | string         | text ccp4mg    |
|                |                |                | scene          |
|                |                |                | description    |
+----------------+----------------+----------------+----------------+
| -              | sceneFile      | string         | Full path      |
|                |                |                | filename for   |
|                |                |                | file           |
|                |                |                | containing     |
|                |                |                | ccp4mg scene   |
|                |                |                | description    |
+----------------+----------------+----------------+----------------+
| -              | label          | string         | Label for the  |
|                |                |                | picture        |
+----------------+----------------+----------------+----------------+
| addCopy        | -              | -              | Return a Copy  |
|                |                |                | object that    |
|                |                |                | copies xml     |
|                |                |                | 'as-is' into   |
|                |                |                | the report     |
+----------------+----------------+----------------+----------------+
| -              | select         | string         | xpath          |
|                |                |                | expression to  |
|                |                |                | specify        |
|                |                |                | element in     |
|                |                |                | task XML       |
+----------------+----------------+----------------+----------------+
| addLaunch      | -              | -              | Return a       |
|                |                |                | Launch object  |
|                |                |                | that puts a    |
|                |                |                | Launch button  |
|                |                |                | in the report  |
+----------------+----------------+----------------+----------------+
| -              | label          | string         | Label for      |
|                |                |                | launch button  |
+----------------+----------------+----------------+----------------+
| -              | taskName       | string         | Name of task   |
|                |                |                | to be launched |
+----------------+----------------+----------------+----------------+
| addHelp        | -              | -              | Add help link  |
|                |                |                | to report      |
+----------------+----------------+----------------+----------------+
| -              | label          | string         | Label for help |
|                |                |                | link           |
+----------------+----------------+----------------+----------------+
| -              | ref            | string         | path of file   |
|                |                |                | to link to     |
+----------------+----------------+----------------+----------------+
| addReference   | -              | -              | Add            |
|                |                |                | bibliographic  |
|                |                |                | reference to   |
|                |                |                | appear in      |
|                |                |                | 'References'   |
|                |                |                | section of     |
|                |                |                | report         |
+----------------+----------------+----------------+----------------+
| -              | authors        | string         | List of        |
|                |                |                | authors        |
+----------------+----------------+----------------+----------------+
| -              | articleTitle   | string         | Title of       |
|                |                |                | article        |
+----------------+----------------+----------------+----------------+
| -              | source         | string         | Source of      |
|                |                |                | article        |
+----------------+----------------+----------------+----------------+
| -              | articleLink    | string         | link to        |
|                |                |                | article on the |
|                |                |                | web            |
+----------------+----------------+----------------+----------------+
| **Repor        |                | -              | Obligatory top |
| t(Container)** |                |                | container for  |
|                |                |                | a report       |
+----------------+----------------+----------------+----------------+
| \__init_\_     | -              | -              | -              |
+----------------+----------------+----------------+----------------+
| -              | standardise    | Boolean        | If true add    |
|                |                |                | the standard   |
|                |                |                | features       |
|                |                |                | (input/output  |
|                |                |                | files etc) to  |
|                |                |                | the report     |
+----------------+----------------+----------------+----------------+
| -              | xmlFile        | string         | File path for  |
|                |                |                | XML file to be |
|                |                |                | read and used  |
|                |                |                | as xmlnode     |
+----------------+----------------+----------------+----------------+
| -              | jobId          | string         | jobId for a    |
|                |                |                | previously run |
|                |                |                | job whose info |
|                |                |                | can be         |
|                |                |                | extracted from |
|                |                |                | database and   |
|                |                |                | used as        |
|                |                |                | jobInfo        |
+----------------+----------------+----------------+----------------+
| -              | htmlBase       | string         | Optional URL   |
|                |                |                | for            |
|                |                |                | *report_files* |
|                |                |                | directory      |
|                |                |                | containing css |
|                |                |                | and js files - |
|                |                |                | default is for |
|                |                |                | CCP4i2 HTTP    |
|                |                |                | server         |
+----------------+----------------+----------------+----------------+
| -              | cssFile        | string         | Optional       |
|                |                |                | basename of a  |
|                |                |                | css file       |
|                |                |                | specific for   |
|                |                |                | this report    |
+----------------+----------------+----------------+----------------+
| -              | jsFile         | string         | Optional       |
|                |                |                | basename of a  |
|                |                |                | js file        |
|                |                |                | specific for   |
|                |                |                | this report    |
+----------------+----------------+----------------+----------------+
| -              | cssVersion     | string         | The version    |
|                |                |                | sub-directory  |
|                |                |                | of             |
|                |                |                | ccp4i          |
|                |                |                | 2/report_files |
|                |                |                | used to        |
|                |                |                | provide css    |
|                |                |                | and js files - |
|                |                |                | defaults to    |
|                |                |                | the latest     |
|                |                |                | version        |
+----------------+----------------+----------------+----------------+
| as_etree       | -              | -              | Return an      |
|                |                |                | etree          |
|                |                |                | representation |
|                |                |                | of the html    |
|                |                |                | report         |
+----------------+----------------+----------------+----------------+
| as_html        | -              | -              | Return html    |
|                |                |                | text           |
|                |                |                | representation |
|                |                |                | of the report  |
+----------------+----------------+----------------+----------------+
| as_html_file   | -              | -              | Write html     |
|                |                |                | representation |
|                |                |                | to a file      |
+----------------+----------------+----------------+----------------+
| -              | fileName       | string         | Full path file |
|                |                |                | name           |
+----------------+----------------+----------------+----------------+
| -              | pretty_print   | Boolean        | If true output |
|                |                |                | in 'pretty     |
|                |                |                | print' more    |
|                |                |                | readable       |
|                |                |                | format         |
+----------------+----------------+----------------+----------------+
| **Di           | -              | -              | Create an html |
| v(Container)** |                |                | div element -  |
|                |                |                | used to group  |
|                |                |                | other elements |
+----------------+----------------+----------------+----------------+
| \__init_\_     | -              | -              | -              |
+----------------+----------------+----------------+----------------+
| **Fol          | -              | -              | Folder in      |
| d(Container)** |                |                | report that    |
|                |                |                | can be         |
|                |                |                | opened/closed  |
|                |                |                | by clicking    |
|                |                |                | title line     |
+----------------+----------------+----------------+----------------+
| \__init_\_     | -              | -              | -              |
+----------------+----------------+----------------+----------------+
| -              | label          | string         | Label for      |
|                |                |                | title line     |
+----------------+----------------+----------------+----------------+
| **GraphGrou    | -              | -              | Group graphs   |
| p(Container)** |                |                | that appear in |
|                |                |                | one graph      |
|                |                |                | widget         |
+----------------+----------------+----------------+----------------+
| **PictureGrou  | -              | -              | Group of       |
| p(Container)** |                |                | pictures to be |
|                |                |                | displayed      |
|                |                |                | together and   |
|                |                |                | put in one     |
|                |                |                | CCP4mg         |
|                |                |                | presentation.  |
|                |                |                | **Not fully    |
|                |                |                | implemented.** |
+----------------+----------------+----------------+----------------+
| **Table**      | -              | -              | Display a      |
|                |                |                | table          |
+----------------+----------------+----------------+----------------+
| addData        | -              | -              | Specify one    |
|                |                |                | column of data |
|                |                |                | for table      |
+----------------+----------------+----------------+----------------+
| -              | title          | string         | column title   |
+----------------+----------------+----------------+----------------+
| -              | subTitle       | string         | sub-title for  |
|                |                |                | tables with    |
|                |                |                | nested headers |
+----------------+----------------+----------------+----------------+
| -              | select         | string         | xpath          |
|                |                |                | expression to  |
|                |                |                | select data    |
|                |                |                | elements from  |
|                |                |                | xml            |
+----------------+----------------+----------------+----------------+
| -              | data           | list of        | List of data   |
|                |                | floats/ints    | elements       |
+----------------+----------------+----------------+----------------+
| **Graph**      | -              | -              | Display a      |
|                |                |                | graph          |
+----------------+----------------+----------------+----------------+
| addData        | -              | -              | Add the data   |
|                |                |                | for one plot   |
|                |                |                | on the graph   |
+----------------+----------------+----------------+----------------+
| -              | xmldata        | List of        | List of        |
|                |                | etree.Elements | elements       |
|                |                |                | (probably in   |
|                |                |                | the task XML)  |
|                |                |                | from which the |
|                |                |                | plot data is   |
|                |                |                | taken          |
+----------------+----------------+----------------+----------------+
| -              | title          | string         | A legend for   |
|                |                |                | the plot       |
+----------------+----------------+----------------+----------------+
| -              | select         | string         | An xpath       |
|                |                |                | expression to  |
|                |                |                | extract data   |
|                |                |                | from           |
|                |                |                | sub-elements   |
|                |                |                | of xmldata     |
+----------------+----------------+----------------+----------------+
| -              | expr           | string         | A Python       |
|                |                |                | expression to  |
|                |                |                | be applied to  |
|                |                |                | each data      |
|                |                |                | element to     |
|                |                |                | transform it   |
+----------------+----------------+----------------+----------------+
| -              | data           | List of        | The data to be |
|                |                | int/float      | plotted        |
+----------------+----------------+----------------+----------------+
| addTable       | -              | -              | Add a table of |
|                |                |                | data taken     |
|                |                |                | from the XML   |
|                |                |                | as a           |
|                |                |                | s              |
|                |                |                | pace-separated |
|                |                |                | table          |
+----------------+----------------+----------------+----------------+
| -              | xmlnode        | etree.Element  | An element     |
|                |                |                | (probably in   |
|                |                |                | the task XML)  |
|                |                |                | containing     |
|                |                |                | table of data  |
+----------------+----------------+----------------+----------------+
| -              | select         | string         | xpath          |
|                |                |                | expression to  |
|                |                |                | sub-element of |
|                |                |                | xmlnode which  |
|                |                |                | is the actual  |
|                |                |                | table          |
+----------------+----------------+----------------+----------------+
| -              | headers        | string         | xpath          |
|                |                |                | expression to  |
|                |                |                | a list of      |
|                |                |                | sub-elements   |
|                |                |                | containing the |
|                |                |                | column header  |
+----------------+----------------+----------------+----------------+
| addPlot        | -              | -              | Add a plot     |
|                |                |                | definition to  |
|                |                |                | the graph      |
+----------------+----------------+----------------+----------------+
|                | xmlnode        | etree.Element  | An element     |
|                |                |                | (probably in   |
|                |                |                | the task XML)  |
|                |                |                | containing     |
|                |                |                | instructions   |
|                |                |                | for one plot   |
+----------------+----------------+----------------+----------------+
| -              | select         | string         | xpath          |
|                |                |                | expression to  |
|                |                |                | sub-element of |
|                |                |                | xmlnode which  |
|                |                |                | is the actual  |
|                |                |                | plot           |
+----------------+----------------+----------------+----------------+
| -              | plot           | string         | An             |
|                |                |                | etree.Element  |
|                |                |                | which is the   |
|                |                |                | plot           |
|                |                |                | definition     |
+----------------+----------------+----------------+----------------+
| -              | plotFile       | string         | The full path  |
|                |                |                | of file        |
|                |                |                | containing a   |
|                |                |                | plot           |
|                |                |                | definition     |
+----------------+----------------+----------------+----------------+
| addPlotObject  | -              | -              | Return a new   |
|                |                |                | *Plot* object. |
+----------------+----------------+----------------+----------------+
| **Plot**       | -              | -              | Define one     |
|                |                |                | plot to appear |
|                |                |                | on a graph     |
+----------------+----------------+----------------+----------------+
| append         | -              | -              | Add an element |
|                |                |                | to the XML     |
|                |                |                | plot           |
|                |                |                | definition     |
+----------------+----------------+----------------+----------------+
| -              | tag            | string         | The XML        |
|                |                |                | element tag    |
+----------------+----------------+----------------+----------------+
| -              | text           | string         | The XML        |
|                |                |                | element text   |
+----------------+----------------+----------------+----------------+
| -              | \*             | \*             | Additional     |
|                |                |                | arguments      |
|                |                |                | interpreted as |
|                |                |                | XML element    |
|                |                |                | attributes     |
+----------------+----------------+----------------+----------------+
| validate       |                |                | Validate the   |
|                |                |                | plot           |
|                |                |                | definition     |
|                |                |                | against the    |
|                |                |                | schema         |
+----------------+----------------+----------------+----------------+
| **Picture**    | -              | -              | Display a      |
|                |                |                | picture        |
+----------------+----------------+----------------+----------------+

File paths
^^^^^^^^^^

Specification of file paths can begin '$CCP4I2' for files in the CCP4i2
distribution directory.

Tables
^^^^^^

The table layout and headers are specified by the *Container.addTable()*
and *Table.addData()* methods, the values are extracted from the task
XML. The table object has two optional arguments, *transpose* changes
the orientation of the table, and *select*, an XPATH expression which
extracts a list of elements containing row data. Each extracted element
is expected to contain all the information for a single table row.
Occasionally the data in the task XML may not be in a suitable form for
this, in which case the argument is omitted. The *addData()* method can
take four arguments: *title*, *subtitle*, *expr* and *select*. The
*title* argument gives the column header. The *subtitle* argument allows
tables with nested headers. The *select* argument specifies an xpath
expression to select the data to be displayed. If the *addTable()*
*select* argument has selected a list of row elements then this should
specify a unique element within that row element. If *addTable()* has no
*select* argument, the select should return a list of elements which
will become the column values. In the latter form there is no guarantee
columns will be of the same length. The *expr* argument enables a python
expression to be applied to the value before display.

Alternatively the Python script can extract and process the data as it
needs and then load it using the *data* argument of the *addData()*
method; this should be a list of values for a column of data.

Graphs
^^^^^^

This is an interface to the **pimple** graph display program which will
be shown inline within the report. *pimple* can display multiple plots
in one window with two combo-boxes to allow selection of the current
plot. The left-hand combo-box will select a data table and the
right-hand box will select between views on that table. The *Graph*
object can only specify one data table with, optionally, multiple views
specified by the *addPlot* block. Multiple *Graph* objects can be
grouped in a *GraphGroup* object to be displayed in the same *pimple*
instance with the *title* arguments from each *addGraph()* appearing in
the left-hand combo-box.

The plotting directives for *pimple* can be entered as XML using the
*addPlot()* method or they can be specified in Python using the
*addPlotObject()* method. The latter returns a *Plot* object which can
then be populated using the *append()* method; this is a totally generic
method for defining the XML definition from Python with arguments that
specify the tag, text and an arbitrary number of attributes. The XML
data structure for a table and associated graphs is described
`here <./pimple.html>`__ . The *Plot.validate()* method will validate
the generated XML against the schema
`ccp4i2/pimple/CCP4ApplicationsOutput.xsd <../../pimple/CCP4ApplicationOutput.xsd>`__
and returns a list of errors as a
`CErrorReport <./error_handling.html>`__. There is an example of this in
`buccaneer_build_refine <../../tasks/buccaneer_build_refine/buccaneer_build_refine_report.py>`__.

There can be multiple *Plot*\ s within one *Graph* - this implies that
multiple views of the same data table are possible and are listed in the
right-hand combo-box in *pimple*.

If the graph is instantiated with *launcher* parameter set then the
graph will not be shown in the report but there will be a button to
launch the graph viewer. The label on the button will be the value of
*launcher*.

Pictures
^^^^^^^^

This makes a picture in the report with buttons to launch either CCP4mg
or Coot with the same view. The picture is generated by CCP4mg which
runs in background without appearing on the computer screen. CCP4mg
needs a *scene.xml* file to specify what it is to display and save as a
picture. The scene file specifies which data files (either coordinates
or maps) are to be loaded and how they should be drawn. The drawing
specification can be either a *simple scene* or generated using the
CCP4mg Picture Wizard. Many other aspects of the display can also be
specified using CCP4mg's usual XML. The template for the scene.xml file
should be given as the *scene* argument or provided in a separate file
that is specified by the *sceneFile* argument. Typically some parameters
such as filenames need to be substituted into the scene.xml template and
this should be done using *database* or *select* attributes to specify
the source of the substituted data.

The contents of the CCP4mg scene file are explained in `Appendix
B <#appendix%20b>`__.

Launch
^^^^^^

Insert a button in the report to launch another viewer. Currently
supports CCP4mg and Coot. The initial view in CCP4mg can be set by
specifying a *pictureDefinitionFile* argument.

Help
^^^^

This puts a *Help* link by tables, graphs etc. and is highly
recommended! The *ref* argument gives the link reference which can be a
file path starting with an environment variable (eg $CCP4I2) and the
appropriate full path will be substituted into the html.

Bibliographic references
^^^^^^^^^^^^^^^^^^^^^^^^

A section of references is automatically appended to a report provided
that the references information has been provided. It is necessary for
task developers to provide MedLine and BibTeX format references for each
task as described `here <./pipelines.html#bibliographic_references>`__.
The *Report* sub-class should call the *addTaskReferences()* method. If
it is called without an argument this method will look for the
bibliographic reference files for the current task. A pipeline report
generator will probably need to call this method to get the references
for sub-tasks, for example for the buccaneer_build_refine pipeline that
uses Refmac the report generator includes:

::

         self.addTaskReferences('buccaneer')
         self.addTaskReferences('prosmart_refmac')

Alternatively the information can be provided via the
*Report.addReference()* method but since this implies there is no
MedLine or BibTeX formated files available the tools for users to
download references will not be available.

Accessing other program output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Within i2, by default, the Report.__init__() is passed the etree
representation of one data XML file, *program.xml* from the top level
job directory. If you are creating a report for a pipeline you may need
to access data files from sub-jobs or data from other sources. The
report generator can access other files using the usual Python
functionality but it must know where to find the files. Within i2 the
Report.__init__() is passed a Python dictionary *jobInfo* which contains
some data extracted from the database (see `Appendix
A <#appendix%20a>`__). The most useful thing is *jobInfo['fileroot']*
that is the job directory for the pipeline job.

If you are running in test mode outside i2 *jobInfo* is empty. See
`Testing <#testing>`__ on how you should run one job and note the
*jobId* that you then pass into the report and it gets the *jobInfo*
from the database.

A potential problem is that the report class does not know exactly which
subjobs were run - an alternative approach is for the pipeline script to
merge the subjob xml files - there is an example in
`pipelines/bucref/script/buccaneer_build_refine <../../tasks/buccaneer_build_refine/buccaneer_build_refine_report.py>`__.
This example uses xpath to pick bits out of the Refmac subjobs xml
output but you could just aim to put everything in one file with a
structure like:

::

   <thewholelot>
   <pointless>
      copy everything from subjob output
   </pointless>
   <aimless>
     copy everything from subjob output
   </aimless>
      and more
   </thewholelot>

It is the responsibility of the developer to deal with the output from
different versions of CCP4i2 tasks and the different programs. Ideally
the report code in any release must be able to deal with the output from
any earlier version of the task or programs - so I'd suggest:

-  program output has some meta-data header that includes the program
   version
-  put each element definition (for graph,table,etc) in a separate
   method - replicate the method if the given element changes
   significantly
-  use try,except where the required data might not be available
-  python tricks such as 'import whatever-3.1 as whatever' could be used
   to handle big changes

Customising Report Appearance
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The appearance of reports is determined by css and Javascript files that
are referenced in the header of the generated html file. Individual
report html elements or groups of html elements for selection by
Javascript or css can be given *id* and/or *class* attributes by using
*id* or *class\_* arguments to the element's *\__init_\_* method or to
the *Container.addWhatever* methods. Note the underscore in the
*class\_* argument that is necessary to differentiate from a Python
reserved word. For example in the *buccaneer_build_refine_report.py*
there is the line:

::

   graph = parent.addGraph( title="Completeness by build-refine cycle", select="//BuccaneerBuildRefineResult/BuildRefineCycle", id = "completenessGraph")

This will create a graph element in the html with the id
*completenessGraph*.

By default the css *xreport.css* and Javascript *xreport.js* files are
referenced in generated html files. These files have a pathname in the
CCP4i2 release such as *ccp4i2/docs/report_files/x.y.z/xreport.css*
where *x.y.z* is a version directory.Developers may add alternative css
or js files but please consider and discuss how your requirements might
be met by updating the default files. If you wish to use alternative
files then theg the filename (e.g. *whatever.css*) must be passsed as
argument *cssFile* or *jsFile* to *Report.__init__()*.

Future releases of CCP4i2 will retain old *report_files* version
directories in order to support any old report files. By default the
report parser will generate html to use the latest version available
unless *Report.__init__()* class has a *cssVersion* argument which is
set to an earlier version number. Developers are recommended to set this
parameter to the version that they have used in development of the *prd*
and to only update the value after they have tested against any more
recent version of *report_files*. If the *cssVersion* is set but there
is a higher minor version number available then that will be used on the
assumption that this minor version number increments are bug-fixes.

CCP4i2 uses an internal HTTP server to serve the files from
*report_files* and it usually expects these to have a URL beginning
*http://127.0.0.1:43434/report_files*. The port number (default 43434)
may change - reports are edited at display time to correct the port
number. Old report files which do not use the HTTP server are also be
updated to do so.

Real-time updating reports
~~~~~~~~~~~~~~~~~~~~~~~~~~

A report can be generated and updated as a job is running - a drawback
is that at the moment the entire report is redone so speed is a possible
issue and the 'running' report should be not too long - usually an
updating graph is ideal. The report class should have the the class
parameter *RUNNING* set true so the system knows to call the report
generator as the job is running. The *Report.__init__()* method is
passed a *jobStatus* argument which will be *Running* and the
appropriate report should be returned.

By default the system watches a running job's *program.xml* file and
when this file changes there is a redraw of the job report. If you
require to watch an alternative file then set the class attribute
*MyReport.WATCHED_FILE* to the name of the watched file or a relative
path to the file from the job directory.

Creating Reports for Pipelines
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In principle a report from a pipeline is no different from a report from
a wrapper but a pipeline report may want to incorporate some elements
from the sub-job reports. For this to be possible the sub-job report
should be written in modular fashion with separate methods for each
element in the report. The *Report.__init__()* method calls the
appropriate methods to generate the report dependent on the *jobStatus*
argument which is usually either *Finished*, *Unsatisfactory* or
*Running*. The *Report.__init__()* method should also handle a
*jobStatus* value of *noOutput* for which it returns without creating
any report elements. For example:

::

   class MyWrapperReport(Report):
     def __init__(self,xmlnode=None,jobInfo={},jobStatus=None,**kw):
       Report. __init__(self,xmlnode=xmlnode,jobInfo=jobInfo,**kw)
       if self.errorReport().maxSeverity()>SEVERITY_WARNING:
         print 'FAILED instantiating MyWrapper report generator'
         self.errorReport().report()
         return

       # 'nooutput' mode would be used by another report class that wanted
       # to use some method(s) from this class for its own report
       if jobStatus is not None and jobStatus.lower() == 'nooutput':
         return
       elif jobStatus == 'Running':
         self.completenessGraph(self)
       else:
         self.completenessGraph(self)
         self.details(self)
         self.picture(self)

     def completenessGraph(self,parent):
       graph = parent.addGraph( title="Completeness by cycle", select="//MyWrapper/Cycle" )
       graph.addData (title="Cycle",  select="Number" )
       ...
       ...

Note that here the method *completenessGraph* is passed a reference to
*self* so it will (in the first line) add a graph to *self* and then go
on to add details to that graph.

A pipeline that uses the *MyWrapper* wrapper and wanted to incorporate
elements from the *MyWrapperReport* could then do:

::

      import MyWrapperReport
      myWrapperXml = self.xmlnode.xpath0('foo/bar')
      myWapperReport = MyWrapperReport.MyWrapperReport(xmlnode=myWrapperXml,jobStatus='noOutput')
      myWapperReport.picture(parent=self)

Here we are assuming that the program xml for *MyWrapper* is embedded in
the pipeline xml and can be selected by the xpath 'foo/bar' and this is
then passed to *MyWrapperReport.__init__()* as the argument *xmlnode*.
Alternatively data may come from a separate file that should be passed
as the *xmlFile* argument. *MyWrapperReport.__init__()* is called with
*jobStatus='noOutput'* so that it creates no output. The pipeline report
generator can then call the appropriate methods for the report elements
that it requires. It passes *self* as the *parent* argument so the
elements are appended to the pipeline report.

One possible issue with calling Report methods from different places is
that the top level of *Report.xmlNode*, the etree representation of the
XML data file, may not be at the same level in the hierarchy and this
will make writing xpath selections tricky. The solution may be to either
add or remove a top level etree element. Here are two examples assuming
we want the top level element to be 'MYWRAPPER'. To remove some
redundant parent elements:

::

     if not self.xmlNode.tag == 'MYWRAPPER':
       myWrapperNode = self.xmlNode.xpath0('MYWRAPPER')
       if myWrapperNode is not None:
         self.xmlNode = myWrapperNode
       else:
         print "ERROR no MYWRAPPER node"

Alternatively to add an extra parent element:

::

     from lxml import etree
     newNode = etree.Element('NEWNODE')
     newNode.append(self.xmlNode)
     self.xmlNode = newNode

Examples of Python Report Definitions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Currently best example:

| 
| `buccaneer_build_refine <../../tasks/buccaneer_build_refine/buccaneer_build_refine_report.py>`__
  which creates 'Running' reports.
| a fragment for
  `pointless <../../test/report_test/pointless_report.py>`__

Please ask Liz if you would like examples of other specific features.

Testing A Python Report Definitions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If the following script fragment (with the appropriate name of the
report class replacing 'my_task_report') is appended to the python file
the script can be run from a Python prompt.

::

   def test(xmlFile=None,jobId=None,reportFile=None):
     import sys,os
     if reportFile is None:
       if xmlFile is not None:
         reportFile = os.path.join(os.path.split(xmlFile)[0],'report.html')
       else:
         reportFile = os.path.join(os.getcwd(),'report.html')
     r = my_task_report(xmlFile=xmlFile,jobId=jobId)
     r.as_html_file(reportFile)
     if len(r.errorReport())>0: print 'ERRORS:',r.errorReport()

It is easiest to test if you have previously run the task from the gui
so that there is program output xml and a jobId that can be input to the
test. The command line is:

::

   > CCCP4I2_TOP/bin/pyi2

   >>>  import my_task_module
   >>>  my_task_module.test(my_xmlfile,my_jobid)

The xml file should be the output xml from a job run previously and the
jobId (not the same as the job number which is displayed in the project
viewer) can be found in the meta-data header of files such as
*params.xml* and *input_params.xml*). This will create report called
*report.html* in the working directory.

Alternatively to test run the report parser in the GUI where you can see
the result better: after running an appropriate job once, close the GUI
while it is displaying the report for the job, delete the report.html
file for that job and reopen the GUI; it should recreate and display the
report.

A good current working example is in
tasks/buccaneer_build_refine/buccaneer_build_refine_report.py

Note that in the *test()* method above the *Report.errorReport()* method
is used to return a `CErrorReport <./error_handling.html>`__ which is a
list of errors and warnings recorded for the Report and its child
objects. If this has non-zero length then it is printed out.

Appendix A - The Database data accessible in the jobInfo dictionary
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
   fileroot - the job directory
   descendentjobs - a list of jobIds

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

This template requires only one piece of substituted data - the *XYZOUT*
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
cross-reference to the *id* attributes of the sub-elements of the
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
then tells CCP4mg to try to find the best orientation to show a
specified atom selection. *scale_auto* is a flag to CCP4mg to scale to
zoom of the specified atom selelction rather than show all visible
atoms.

Alternatively the scale,orientation and centring can be specified
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

Appendix C - Style Guide
~~~~~~~~~~~~~~~~~~~~~~~~

Firstly the 'running' reports that are updated in real time usually
consist of a graph that shows, for example, progress per cycle and
possibly a table showing initial, previous and current cycle values.

Final reports should begin with a short text summary of what the task
achieved, any issues and warning messages in red. The terminology here
(and throughout the report) should be consistent with the terminology
used in the gui and documentation. Links to documentation are
encouraged.

There is usually an open folder showing a graphical summary of the task
followed by closed folders with more detailed information possibly from
the separate steps of the task.

Visual presentation (i.e. graphs) is preferred over tables - thought
tables presenting the same information may be good. The usual good
practice for graphs - label axes, legends etc. apply. A simple clear
graph or multiple simple clear graphs are preferred over over-loaded
graphs (though beware what users might want to look at simultaneously.)

Last modified: Tue Jun 2 12:32:24 BST 2015
