################################
Example Code for *CData* Classes
################################

Introduction
~~~~~~~~~~~~

Use the `console mode <./introduction.html#console>`__ to try these. I
expect that plugins will mostly initialise data structures from DEF and
DATA files so you will not need to do this basic initialisation stuff
often. Please ask if you need explanation or other examples or feel free
to edit in more examples.

Getting/Setting Data
~~~~~~~~~~~~~~~~~~~~

::

   from core.CCP4Data import *
   i = CInt(5)
   print type(i),i        >>>>>   <class 'CCP4Data.CInt'> 5
   i.set(7)
   print type(i),i        >>>>>   <class 'CCP4Data.CInt'> 7
   i = 7
   print type(i),i        >>>>>   <type 'int'> 7     (Beware!!!)

   # A 'complex' class with two items - several equivalent ways to initialise..
   r = CFloatRange(start=5,end=10)
   r = CFloatRange(value={'start':5,'end':10})
   r = CFloatRange(name='foo',value={'start':5,'end':10})
   print type(r),r        >>>>>   <class 'CCP4Data.CFloatRange'> {'start': 5.0, 'end': 10.0}

   # addressing the individual items
   r = CFloatRange()
   r.start= 5.0
   r.end = 10.0
   print type(r.start),r.start         >>>>>  <class 'CCP4Data.CFloat'> 5.0

   # When you need getattr (or setattr)
   for item in ['start','end']:
     sum = sum  + getattr(r,item)

Qualifiers
~~~~~~~~~~

::


   i = CInt(min=10)
   i.set(5)             >>>>>  'below minimum' exception

   r = CFloatRange(start_min=0.0,end_max=10.0,compare=1)   # Note setting qualifiers for items
   r.start = 5.0
   r.end = 2.0          >>>>> 'End of range less than start' exception
   r.end = 12.0         >>>>> 'above maximum' exception
   r.get()              >>>>>  { 'start' : 5.0, 'end' : None }


   # List all qualifers
   r.qualifiers()       >>>>>  {'compare': 1}
   r.start.qualifiers() >>>>>   {'charWidth': None, 'min': 0.0, 'default': None, 'max': None, 'allowUndefined': True, 'onlyEnumerators': False, 'enumerators': [], 'menuText': []}
   # List qualifiers but not those that have default values
   r.start.qualifiers(default=False)  >>>> {'min': 0.0}
   # Get a specific qualifier
   r.start.qualifier('min')   >>>>>    0.0

Maths
~~~~~

NOTE: The CData classes do a lot of work validating after ever value
update. So probably not wise to use them for serious calculations.

::


   s = r.start + 4
   print r.start,type(s),s    >>>>  5.0  9.0

   r.start.add(2.5)
   print r                  >>>>>  {'start': 7.5, 'end': None}

   r.start>3                >>>>> True

CDataFile
~~~~~~~~~

| Initialising and retrieving a file name. Applies to subclasses such as
  CPdbDataFile and CMtzDataFile
| There is a project 'my_project' with directory
  /Users/lizp/Desktop/myproject.
| N.B. CDataFile.__str__() is best way to return a string with full file
  name. e.g. *str(f)*

::

   from core.CCP4File import *

   f = CDataFile('/Users/lizp/Desktop/whatsit.wht')
   print f                 >>>>>   /Users/lizp/Desktop/whatsit.wht
   str(f)                  >>>>>   /Users/lizp/Desktop/whatsit.wht 
   f.get()                 >>>>>   {'project': None, 'baseName': 'whatsit.wht', 'relPath': '/Users/lizp/Desktop'}

   f = CDataFile(project='my_project',relPath='job_2',baseName='whatsit.wht')
   str(f)                  >>>>>    /Users/lizp/Desktop/myproject/job_2/whatsit.wht

   f = CDataFile('/Users/lizp/Desktop/myproject/job_2/whatsit.wht')
   f.get()                 >>>>>    {'project': 'my_project', 'baseName': 'whatsit.wht', 'relPath': 'job_2'}
   f.relPath = 'job_66'
   f.get()                 >>>>>    {'project': 'my_project', 'baseName': 'whatsit.wht', 'relPath': 'job_66'}

Note that *CDataFile.dbFileId* is an integer id for the file in the
database and is vital for the core i2 code that saves the file
provenance to the database. If you are copying filenames in pipelines
(for example setting the input to one wrapper from the output of
another) please copy the whole of the CDataFile object:

::

     secondPlugin.container.inputData.XYZIN = firstPlugin.container.outputData.XYZOUT

Here, in a pipeline, *firstPlugin* and *secondPlugin* are two
consecutive wrappers and *XYZOUT* and *XYZIN* are instances of
*CPdbDataFile* that are the output and input of the two plugins. This is
explained in the `pipelines <./pipleines.html>`__ section.

CI2XmlFile - XML Files
~~~~~~~~~~~~~~~~~~~~~~

| Reading an XML file. Note that the 'CCP4I2' directory alias can be
  used as 'project'
| The body of the file can be return as an etree which can be addressed
  using the functions from lxml

::

   from core.CCP4File import *
   f = CI2XmlDataFile(project='CCP4I2',relPath='test/data',baseName='pdbset.def.xml')
   f.header.get()           >>>>>  {'function': 'DEF', 'comment': None, 'creationTime': None,
                                    'userId': 'lizp', 'ccp4iVersion': '0.1', 'jobId': None,
                                    'project': None, 'pluginName': 'pdbset',
                                    'pluginVersion': None, 'pluginTitle': 'PDBSet'}
   root = f.getBodyEtree()
   # Use lxml functionality to address the body of the xml file..
   for child in root:  print child.get('id')     >>>>>>>>>   inputData
                                                             outputData
                                                             controlParameters

Creating an XML file

::

   from core.CCP4File import *
   from lxml import etree

   f = CI2XmlDataFile ('/Users/lizp/Desktop/whatsit.xml')
   f.header.setCurrent()    #  Sets the current time and user in header
   f.header.function = 'DEF'

   tree = etree.Element('mystuff')
   tree.append(etree.Element('something'))
   tree.append(etree.Element('something_else'))

   f.saveFile(bodyEtree=tree)

Lists
~~~~~

The list base class *CList* is intended to behave like a Python list but
has the restriction that all items in the list must be of the same type
and that type is a subclass of *CData*. The CList has qualifiers
*listMinLength* and *listMaxLength* which do what they say and
*listCompare* which is an integer - if it is defined with value 1/-1
then each item in the list must be greater/less than preceeding item -
for this to work the list item class must have a \__cmp__() method (`see
Python docs <http://docs.python.org/library/functions.html#cmp>`__). In
the example below a list of integers with minimum length 1 and minimum
integer value 1 is created.

::

   from core.CCP4Data import *
   l = CList(subItem= { 'class' : CInt, 'qualifiers' : { 'min' : 1 , 'default' : 1 }}, listMinLength=1)
   print l >>>>  [1]
   l.append(6)
   print l >>>>  [1,6]
   l[1] = -1  >>>> below minimum error
   l[1] = l[1] + 2
   print l  >>>>> [1,8]
   del l[1]
   del l[0] >>>>> list below minimum length error

Below is an example of a *CList* subclass from the *CCP4ModelData* -
this has list items *CResidueRange* that has components: chainId and
first and last residue.

::

   from core.CCP4ModelData import *
   l = CResidueRangeList()
   l.append({'chainId':'A','firstRes':'23','lastRes':'45'})
   print l,type(l[0])   >>>> [{'firstRes': '23', 'chainId': 'A', 'lastRes': '45'}] , 
   print l[0].chainId  >>>> A

--------------

Last modified: Thu Jun 14 13:03:14 BST 2012
