##########################################################
CCP4I2 Developers Data Classes and Containers
##########################################################

Data Classes
~~~~~~~~~~~~

All data classes are sub-classed from *CCP4Data.CData*. The same classes
are used by the gui and pipelines. The *CCP4Data.CData* (in core
directory) currently subclasses either\ *CCP4QtObject.CObject* or
*CCP4Object.CObject*. The first of these subclasses the Qt class
*QtCore.QObject* which provides `signals and
slots <http://doc.qt.nokia.com/4.6/signalsandslots.html>`__
functionality. The second is intended as an alternative to Qt dependence
but is not yet implemented.

The *CCP4Data.CBaseData* class sub-classes CData and is the base class
for all 'simple' data classes that hold one item of data.

Currently the data classes are in:

-  CCP4Data - general data
-  CCP4XtalData MX experimental data
-  CCP4ModelData MX model data
-  CCP4Files - files and XML handling
-  CCP4Annotation - annotation and meta data
-  CCP4MathsData - basic maths

There is a GUI widget for each data class - there is not necessarilly a
simple ono-to-one relationship - each widget class has a MODEL_CLASS
parameter which lists which data classes it can represent. The gui
classes are in the *qtgui* directory.

To enable easy XML input/output of data all data classes have methods to
get and set an eTree representation of their data content. The eTree
functionality is provided by *lxml* - see http://codespeak.net/lxml/.The
eTree representation is easily imported from or exported as XML.

All CData classes have a few class-wide attributes that are usually
defined at the top of the class. These are:

+-----------------------+----------------------+----------------------+
| **Attribute**         | **Format**           | **Function**         |
+-----------------------+----------------------+----------------------+
| CONTENTS              | Dictionary, keys are | The contents of      |
|                       | sub-object names and | 'complex' classes    |
|                       | values are           |                      |
|                       | dictionaries with at |                      |
|                       | least a *class* key. |                      |
+-----------------------+----------------------+----------------------+
| PYTHONTYPE            | A Python type        | The Python type for  |
|                       |                      | simple data classes  |
+-----------------------+----------------------+----------------------+
| QUALIFIERS            | Dictionary, keys are | Specify the default  |
|                       | qualifier names, and | qualifiers for this  |
|                       | values are           | class. This can be   |
|                       | appropriate type for | over-ridden for any  |
|                       | the qualifier.       | instance of the      |
|                       |                      | class.               |
+-----------------------+----------------------+----------------------+
| QUALIFIERS_ORDER      | A list of qualifier  | The order qualifiers |
|                       | names                | will apppear in      |
|                       |                      | defEd GUI.           |
+-----------------------+----------------------+----------------------+
| QUALIFIERS_DEFINITION | Dictionary, keys are |                      |
|                       | qualifier names,     |                      |
|                       | values are a         |                      |
|                       | dictionary with at   |                      |
|                       | least keys: *type* - |                      |
|                       | a Python type and    |                      |
|                       | *description* - a    |                      |
|                       | brief description    |                      |
|                       | for defEd GUI.       |                      |
+-----------------------+----------------------+----------------------+
| ERROR_CODES           | Dictionary           | Described above      |
+-----------------------+----------------------+----------------------+
| PROPERTIES            | Dictionary with keys | Analogous to the     |
|                       | are property name    | Python *property*    |
|                       | and value is         | command. This must   |
|                       | dictionary with      | be placed after the  |
|                       | *fget* and *fset*    | methods for *fget*   |
|                       | keys.                | and *fset*.          |
+-----------------------+----------------------+----------------------+

These attributes should be accessed by methods such as *contents()* and
*qualifier()* rather than accessed directly. The access methods
implement the principle that sub-classes inherit the properties of base
classes and that the QUALIFIER attributes can be over-ridden in
particular class instances.

The following table explains the main class methods - it is not
comprehensive documentation.

+----------------------+------------+---------------------------+
| CData method         | Argument   | Description               |
+======================+============+===========================+
| \__init_\_           |            |                           |
+----------------------+------------+---------------------------+
|                      | parent     | The Qt system expects     |
|                      |            | each QObject to have a    |
|                      |            | parent QObject with the   |
|                      |            | first ancestor object     |
|                      |            | being the Qt              |
|                      |            | QApplication. The         |
|                      |            | CC                        |
|                      |            | P4Modules.QTAPPLICATION() |
|                      |            | function will return the  |
|                      |            | QApplication if you need  |
|                      |            | it but the usual parent   |
|                      |            | of CData should be the    |
|                      |            | *CContainer* that holds   |
|                      |            | the data.                 |
+----------------------+------------+---------------------------+
|                      | qualifiers | A dictionary of           |
|                      |            | parameters to 'fine-tune' |
|                      |            | the class - particularly  |
|                      |            | specifying validation     |
|                      |            | parameters such as        |
|                      |            | maximum or minimum        |
|                      |            | allowed values. This      |
|                      |            | mechanism is intended to  |
|                      |            | reduce the need for       |
|                      |            | sub-classing. The allowed |
|                      |            | values can be seen in the |
|                      |            | \__init_\_ method. The    |
|                      |            | CData base class expects  |
|                      |            | two parameters: a *name*  |
|                      |            | and a *default* value.    |
|                      |            | The name is the same as   |
|                      |            | the element tag in the    |
|                      |            | XML representation of the |
|                      |            | class. The default        |
|                      |            | *default* value is the    |
|                      |            | Python *None* (i.e. not   |
|                      |            | set).                     |
+----------------------+------------+---------------------------+
| parent               |            | Return the parent object. |
+----------------------+------------+---------------------------+
| objectName           |            | Return the object's name. |
+----------------------+------------+---------------------------+
| objectPath           |            | Return a 'path' name      |
|                      |            | incuding parent object    |
|                      |            | names. An underscore      |
|                      |            | separator is used.        |
+----------------------+------------+---------------------------+
| getEtree             |            | Return an eTree element   |
|                      |            | representing the data     |
|                      |            | object.                   |
+----------------------+------------+---------------------------+
| setEtree             |            | Parse an eTree element    |
|                      |            | and initialise the data   |
|                      |            | object.                   |
+----------------------+------------+---------------------------+
|                      | element    | An etree element          |
+----------------------+------------+---------------------------+
| set                  |            | Set the data. The method  |
|                      |            | calls the *validity*      |
|                      |            | method and will only set  |
|                      |            | the data if the validity  |
|                      |            | is OK. setData returns a  |
|                      |            | CErrorReport from         |
|                      |            | validity().               |
+----------------------+------------+---------------------------+
|                      | data       | The set method attempts   |
|                      |            | to be very tolerant of    |
|                      |            | input e.g. *i = CInt();   |
|                      |            | i.set('12')* is           |
|                      |            | acceptable. Simple, one   |
|                      |            | item, data objects expect |
|                      |            | the appropriate Python    |
|                      |            | type or the appropriate   |
|                      |            | CData class as argument   |
|                      |            | to set. Complex, multi    |
|                      |            | item, data objects,       |
|                      |            | expect Python             |
|                      |            | dictionaries or the       |
|                      |            | different components on   |
|                      |            | the command line. See the |
|                      |            | `example                  |
|                      |            | code                      |
|                      |            |  <./example_code.html>`__ |
+----------------------+------------+---------------------------+
| unSet                |            | Unset the data            |
+----------------------+------------+---------------------------+
| get                  |            | Return the data contents  |
|                      |            | on the object. For a      |
|                      |            | simple class with single  |
|                      |            | item then that is         |
|                      |            | returned; if the data has |
|                      |            | multiple items then they  |
|                      |            | are returned as a Python  |
|                      |            | dictionary.               |
+----------------------+------------+---------------------------+
|                      | name       | Optional name of one item |
|                      |            | in a multi-item data      |
|                      |            | class. If this is set     |
|                      |            | then only the named item  |
|                      |            | is returned.              |
+----------------------+------------+---------------------------+
| validity             |            | Return a *CErrorReport*   |
|                      |            | indicating validity of    |
|                      |            | input data.               |
+----------------------+------------+---------------------------+
| fix                  |            | Return a dictionary       |
|                      |            | representation of the     |
|                      |            | class data that has been  |
|                      |            | fixed (possibly with      |
|                      |            | significant change to     |
|                      |            | content) to be valid.     |
|                      |            | This method is only       |
|                      |            | useful where it has been  |
|                      |            | reimplemented in some     |
|                      |            | sub-classes.              |
+----------------------+------------+---------------------------+
| isSet                |            | Return *False*/ *True*    |
|                      |            | dependent on whether the  |
|                      |            | data is set or is a null  |
|                      |            | value.                    |
+----------------------+------------+---------------------------+
| emitDataChanged      |            | A wrapper for the Qt      |
|                      |            | QObject.emit() method to  |
|                      |            | emit a 'dataChanged'      |
|                      |            | signal.                   |
+----------------------+------------+---------------------------+
| setQualifiersEtree   |            | Parse an etree element    |
|                      |            | containing the qualifiers |
|                      |            | for the class. This is    |
|                      |            | used when the data type   |
|                      |            | definitions are read from |
|                      |            | an XML file. Return a     |
|                      |            | *CErrorReport*            |
+----------------------+------------+---------------------------+
|                      | element    | etree element.            |
+----------------------+------------+---------------------------+
| getQualifiersEtree   |            | Return a tuple of an      |
|                      |            | eTree element containing  |
|                      |            | the data qualifiers for   |
|                      |            | the class and a           |
|                      |            | *CErrorReport*            |
+----------------------+------------+---------------------------+
| contents             |            | Return a dictionary       |
|                      |            | specifying the contents   |
|                      |            | of a complex class        |
+----------------------+------------+---------------------------+
|                      | name       | Name of one item in the   |
|                      |            | class. If this is set     |
|                      |            | then returns definition   |
|                      |            | of just that item.        |
+----------------------+------------+---------------------------+
| qualifiers           |            | Return a dictionary of    |
|                      |            | qualifier values for the  |
|                      |            | class.                    |
+----------------------+------------+---------------------------+
|                      | name       | The name of a qualifier.  |
|                      |            | Return the value of that  |
|                      |            | qualifier only.           |
+----------------------+------------+---------------------------+
|                      | default    | If False do not return    |
|                      |            | the qualifiers that have  |
|                      |            | the default value for the |
|                      |            | class                     |
+----------------------+------------+---------------------------+
|                      | custom     | If False do not return    |
|                      |            | the qualifiers that have  |
|                      |            | the been customised in    |
|                      |            | this instance of the      |
|                      |            | class                     |
+----------------------+------------+---------------------------+
| qualifiersDefinition |            | Return the definition of  |
|                      |            | the qualifiers for the    |
|                      |            | class                     |
+----------------------+------------+---------------------------+
| pythonType           |            | Return the Python type    |
|                      |            | equivalent for a simple   |
|                      |            | class                     |
+----------------------+------------+---------------------------+
|                      | name       | The name of a qualifier.  |
|                      |            | Return the definition of  |
|                      |            | that qualifier only.      |
+----------------------+------------+---------------------------+
| qualifiersOrder      |            | Return the order of the   |
|                      |            | qualifiers used in a GUI. |
+----------------------+------------+---------------------------+
| setDefault           |            | Set the data to the       |
|                      |            | default value.            |
+----------------------+------------+---------------------------+

The Container Class
~~~~~~~~~~~~~~~~~~~

*CCP4Container.CContainer* is a sub-class of *CData* which differs
mostly in that its data contents (for most\ *CData* defined bt the
CONTENTS attribute) are defined at run time. The container holds a set
of data objects in a Python dictionary. Typically one container is
associated with one gui window or one program wrapper. The data
container classes are not sub-classed - their contents are specified
when the class is instantiated. The data objects within the container
apply their *validity()* method to ensure that all loaded data is valid.
A data container may contain sub-containers.

Typically the content of a *CContainer* is defined in a DEF (extension
*def.xml*) file and its data is imported from and exported to PARAMS
(extension *params.xml*) files. When the content definition is read the
data objects to hold the data are created automatically.

Containers are the Python representation of PARAMS files and these files
contain a header holding meta-data such as creation date and project id.
The header can be accessed as *myContainer.header* and is a
CCP4File:CI2XmlHeader class. Tasks developers should not normally need
to access this data as it is handled by the CCP4i2 core.

+-----------------------+-------------+--------------------------+
| Method                | Argument    | Description              |
+=======================+=============+==========================+
| loadContentsFromXml   |             | Load the content         |
|                       |             | definition from a DEF    |
|                       |             | file. Returns a          |
|                       |             | CErrorReport.            |
+-----------------------+-------------+--------------------------+
|                       | fileName    | The full path name of a  |
|                       |             | file                     |
+-----------------------+-------------+--------------------------+
| loadDataFromXml       |             | Load the data values     |
|                       |             | from a PARAMS file.      |
|                       |             | Should only be used      |
|                       |             | after the contents have  |
|                       |             | been defined.Returns a   |
|                       |             | CErrorReport.            |
+-----------------------+-------------+--------------------------+
|                       | fileName    | The full path name of a  |
|                       |             | file                     |
+-----------------------+-------------+--------------------------+
| saveContentsToXml     |             | Save the content         |
|                       |             | definition to a DEF      |
|                       |             | file.Returns a           |
|                       |             | CErrorReport.            |
+-----------------------+-------------+--------------------------+
|                       | fileName    | The full path name of a  |
|                       |             | file                     |
+-----------------------+-------------+--------------------------+
| saveDataToXml         |             | Save the data values to  |
|                       |             | a PARAMS file.Returns a  |
|                       |             | CErrorReport.            |
+-----------------------+-------------+--------------------------+
|                       | fileName    | The full path name of a  |
|                       |             | file                     |
+-----------------------+-------------+--------------------------+
| loadContentsFromEtree |             | Load the content         |
|                       |             | definition from an eTree |
|                       |             | element                  |
+-----------------------+-------------+--------------------------+
|                       | element     | An eTree element         |
+-----------------------+-------------+--------------------------+
| loadDataFromEtree     |             | Load the data values     |
|                       |             | from an eTree element.   |
|                       |             | Should only be used      |
|                       |             | after the contents have  |
|                       |             | been defined.            |
+-----------------------+-------------+--------------------------+
|                       | element     | An eTree element         |
+-----------------------+-------------+--------------------------+
| saveContentsToEtree   |             | Return a Python tuple of |
|                       |             | an eTree element         |
|                       |             | representing the         |
|                       |             | contents and a           |
|                       |             | CErrorReport.            |
+-----------------------+-------------+--------------------------+
| saveDataToEtree       |             | Return a Python tuple of |
|                       |             | an eTree element         |
|                       |             | representing the data    |
|                       |             | values and a             |
|                       |             | CErrorReport.            |
+-----------------------+-------------+--------------------------+
| addContent            |             | The arguments to this    |
|                       |             | method define a data     |
|                       |             | object. The method       |
|                       |             | creates a new object and |
|                       |             | appends it to the        |
|                       |             | container.               |
+-----------------------+-------------+--------------------------+
|                       | name        | A name for the new data  |
|                       |             | object                   |
+-----------------------+-------------+--------------------------+
|                       | cls         | The class of the new     |
|                       |             | data object              |
+-----------------------+-------------+--------------------------+
|                       | qualifiers  | A Python dictionary of   |
|                       |             | qualifiiers for the new  |
|                       |             | object                   |
+-----------------------+-------------+--------------------------+
| addObject             |             | Add an already existing  |
|                       |             | data object ot the       |
|                       |             | container                |
+-----------------------+-------------+--------------------------+
|                       | name        | A name for the new data  |
|                       |             | object                   |
+-----------------------+-------------+--------------------------+
|                       | object      | The CData object         |
+-----------------------+-------------+--------------------------+
|                       | afterObject | Insert the new object in |
|                       |             | the container after the  |
|                       |             | object with this name.   |
+-----------------------+-------------+--------------------------+
| replaceObject         |             | Replace a data object    |
|                       |             | with given name by       |
|                       |             | another existing object  |
+-----------------------+-------------+--------------------------+
|                       | name        | The name of an existing  |
|                       |             | data object in the       |
|                       |             | container                |
+-----------------------+-------------+--------------------------+
|                       | object      | A CData object           |
+-----------------------+-------------+--------------------------+
| deleteObject          |             | Delete an object with    |
|                       |             | given name               |
+-----------------------+-------------+--------------------------+
|                       | name        | The name of an existing  |
|                       |             | data object in the       |
|                       |             | container                |
+-----------------------+-------------+--------------------------+
| renameObject          |             | Rename an object         |
+-----------------------+-------------+--------------------------+
|                       | oldName     | The name of an existing  |
|                       |             | data object in the       |
|                       |             | container                |
+-----------------------+-------------+--------------------------+
|                       | newName     | A new name               |
+-----------------------+-------------+--------------------------+
| clear                 |             | remove all content from  |
|                       |             | the comtainer            |
+-----------------------+-------------+--------------------------+
| dataOrder             |             | Return a list of all the |
|                       |             | names of data objects    |
|                       |             | (and sub-containers) in  |
|                       |             | the container            |
+-----------------------+-------------+--------------------------+
| addHeader             |             | Add a CCP4File.CHeader   |
|                       |             | header to the container. |
|                       |             | This header can be set   |
|                       |             | appropriately prior to   |
|                       |             | exporting an XML file.   |
+-----------------------+-------------+--------------------------+
| parseCommandLine      |             | Load data into the       |
|                       |             | container by parsing a   |
|                       |             | 'command line' formatted |
|                       |             | asa list of words        |
+-----------------------+-------------+--------------------------+
|                       | commandLine | A list of words such as  |
|                       |             | returned by sys.argv()   |
+-----------------------+-------------+--------------------------+
|                       | template    | Optional template for    |
|                       |             | interpreting command     |
|                       |             | line (`more              |
|                       |             | in                       |
|                       |             | fo <#command%20line>`__) |
+-----------------------+-------------+--------------------------+

Creating a New File Type
~~~~~~~~~~~~~~~~~~~~~~~~

This is a brief how-to - mostly to remind Liz.

To register a new file type in CCP4i2 it needs to be added in three
different places (no, that's not good but at least its documented!):

-  A sub-class of *CCP4DataFile.CDataFile* should be created so there is
   an appropriate data object.
-  It should be added to the mimetypes in *CCP4CustomMimeTypes* so the
   Qt file display and drag-n-drop can handle the files.
-  It should be registered with the database by adding to lists in
   *CCP4DbApi*.
-  Also a file icon should be provided in the *qticons* directory.

A simple file class definition might look like:

::

           class CPhaserSolDataFile(CCP4File.CDataFile):
           
           QUALIFIERS = { 'mimeTypeName' : 'application/phaser-sol',
           'mimeTypeDescription' : 'Phaser solution file',
           'fileExtensions' : [ 'phaser_sol.pkl' ],
           'fileContentClassName' : None,
           'fileLabel' : 'phaser_sol',
           'guiLabel' : 'Phaser solution file',
           'toolTip' : "Possible solutions passed between runs of the Phaser program",
           }
       

Note that the *mimeTypeName* must match the mime type in
*CCP4CustomMimeTypes*, *CCP4DbApi*. The *fileContentClassName* is not
set here but if i2 needs to read the file then a file content class that
sub-classes *CCP4File.CDataFileContent* should be implemented and the
class name provided here as a string.

A definition of the mime type looks like this:

::

           mimeType = CMimeType()
           mimeType.name = "application/phaser-sol"
           mimeType.description = "Phaser solution file"
           mimeType.fileExtensions = ['phaser_sol.pkl']
           mimeType.viewers = []
           mimeType.icon = 'PhaserSolDataFile'
           mimeType.className = 'PhaserSolDataFile'
           self.mimeTypes["application/phaser-sol"] = mimeType
       

Note that the *className* should cross reference the class name and by
default the system assumes the file icon is the class name but an
alternative can be provided here to be used in the file browser. Note
that the className and icon have dropped the leading 'C'.

In the database file *CCP4DbApi* it is necessary to add the file type to
three lists:

-  FILETYPES_TEXT - this is a list of the mime types
-  FILETYPES_CLASS - a list of file class names without the initial 'C'
-  FILETYPELIST - a list of items with each element being a list of
   three items: a file type number, the mime type name and a short
   description.

It is imperative to keep these three lists in sync! Also if a file type
is added to the list and then subsequently not required the slot should
not be deleted or reused - it should just be left as a 'dummy' the same
as file type 14.

File content classes and *contentFlag* and *subType* flags
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The classes derived from *CDataFile* such as *CPdbDataFile* and the
various forms of mini MTZ such as *CObsDataFile* basically serve as a
reference to the data file holding the file path and some additional
flags discussed later. A *fileContent* class such as *CPdbData* or
*CMtzData* can hold some data extracted from the file. A example of
using these:

::

           >>> from core import CCP4XtalData
           >>> mtzData = CCP4XtalData.CMtzDataFile('/y/people/lizp/rnase25.mtz')
           >>> print type(f.fileContent)
           <class 'CCP4XtalData.CMtzData'>
           >>> print f.fileContent.cell
           {'a': '64.8970031738', 'c': '38.7919998169', 'b': '78.3229980469', 'beta': '1.57079632679', 'alpha': '1.57079632679', 'gamma': '1.57079632679'}
       

The file content classes are mostly being developed as required - please
contact Liz if you need further functionality.

The file classes also contain two integer flags, *contentFlag* and
*subType* that indicate something about the content of a selected file.
The *contentFlag* is most important for mini MTZs to indicate the
representation of the data (see `miniMtzs <../general/miniMtzs.html>`__
for explanation).

+-------------------------------------+-------+-------------------------------+
| Name                                | Value | Meaning                       |
+=====================================+=======+===============================+
| **CObsDataFile**                                                            |
+-------------------------------------+-------+-------------------------------+
| CObsDa                              | 1     | Freidal's pairs of            |
| CObsDataFile.CONTENT_FLAG_IPAIR     |       | intensities                   |
+-------------------------------------+-------+-------------------------------+
| CObsDataFile.CONTENT_FLAG_FPAIR     | 2     | Freidal's pairs of            |
|                                     |       | structure factors             |
+-------------------------------------+-------+-------------------------------+
| CObsDataFile.CONTENT_FLAG_IMEAN     | 3     | averaged intensities          |
|                                     |       |                               |
+-------------------------------------+-------+-------------------------------+
| CObsDataFile.CONTENT_FLAG_FMEAN     | 4     | averaged structure            |
|                                     |       | factors                       |
+-------------------------------------+-------+-------------------------------+
| **CPhsDataFile**                                                            |
+-------------------------------------+-------+-------------------------------+
| CPhsDataFile.CONTENT_FLAG_HL        | 1     | Hendrickson-Lattmann          |
|                                     |       | coefficients                  |
+-------------------------------------+-------+-------------------------------+
| CPhsDataFile.CONTENT_FLAG_PHIFOM    | 1     | phase and figure of merit     |
+-------------------------------------+-------+-------------------------------+
| **CMapCoeffsDataFile**                                                      |
+-------------------------------------+-------+-------------------------------+
| CMapCoeffsDataFile.CONTENT_FLAG_FPHI| 1     | structure facture and         |
|                                     |       | phase                         |
+-------------------------------------+-------+-------------------------------+

The *contentFlag* is presently only used for mini-MTZs but, in
principle, could be used to distiguish different forms of other data
files. This flag is useful in automating handling of the mini-MTZs and
ensuring programs are given data in a representation that they can
handle.

The *subType* flag is presently used by mini-MTZ classes and
*CPdbDataFile* to indicate the scientific content of the data.

+--------------------------------------+-------+-----------------------------+
| Name                                 | Value | Meaning                     |
+======================================+=======+=============================+
| **CPdbDataFile**                                                           |
+--------------------------------------+-------+-----------------------------+
| CPdbDataFile.SUBTYPE_UNKNOWN         | 0     | unknown                     |
+--------------------------------------+-------+-----------------------------+
| CPdbDataFile.SUBTYPE_MODEL           | 1     | working model               |
+--------------------------------------+-------+-----------------------------+
| C                                    | 2     | homolog                     |
| PdbDataFile.SUBTYPE_HOMOLOG          |       |                             |
+--------------------------------------+-------+-----------------------------+
| CP                                   | 3     | structure fragment (e.g.    |
| dbDataFile.SUBTYPE_FRAGMENT          |       | ligand)                     |
+--------------------------------------+-------+-----------------------------+
| CPdbDataFile.SUBTYPE_HEAVY_ATOMS     | 4     | heavy atoms                 |
|                                      |       |                             |
+--------------------------------------+-------+-----------------------------+
| **CObsDataFile**                                                           |
+--------------------------------------+-------+-----------------------------+
| CObsDataFile.SUBTYPE_OBSERVED        | 1     | observed data               |
|                                      |       |                             |
+--------------------------------------+-------+-----------------------------+
| CObsDataFile.SUBTYPE_DERIVED         | 2     | derived data                |
|                                      |       |                             |
+--------------------------------------+-------+-----------------------------+
| CObsDataFile.SUBTYPE_REFERENCE       | 3     | reference data              |
+--------------------------------------+-------+-----------------------------+
| **CPhsDataFile**                                                           |
+--------------------------------------+-------+-----------------------------+
| CPhsDataFile.SUBTYPE_UNBIASED        | 1     | unbiased phases             |
|                                      |       |                             |
+--------------------------------------+-------+-----------------------------+
| CPhsDataFile.SUBTYPE_BIASED          | 2     | biased phases               |
+--------------------------------------+-------+-----------------------------+
| **CMapCoeffsDataFile**                                                     |
+--------------------------------------+-------+-----------------------------+
| CMapCoeffsDataFile.SUBTYPE_NORMAL    | 1     | normal map                  |
|                                      |       |                             |
+--------------------------------------+-------+-----------------------------+
| CMapCoeffsDataFile.SUBTYPE_DIFFERENCE| 2     | difference map              |
|                                      |       |                             |
+--------------------------------------+-------+-----------------------------+
| **CCootHistoryDataFile**                                                   |
+--------------------------------------+-------+-----------------------------+
| CCootHistoryDataFile.SUBTYPE_INITIAL | 1     | Initialisation for Coot     |
|                                      |       |                             |
+--------------------------------------+-------+-----------------------------+
| CCootHistoryDataFile.SUBTYPE_HISTORY | 2     | Coot history file           |
|                                      |       |                             |
+--------------------------------------+-------+-----------------------------+

This information can be used to ensure that the data is used
appropriately. Obviously the *subType* can be tricky to define so its
use in the code is intended to be flexible.

If these flags are set in a *CDataFile* class that is being saved to the
database then the *contentFlag* and/or *subType* are saved to the
database. It is very helpful therefore if
*CPluginScript.processOutputFiles()* implementations set these flags for
output files. It is also possible to set a default value for
*contentFlag* or *subType* in a wrapper def file if the nature of the
program output file is known in advance (but done here is probably less
obvious and may not be maintained properly).

It is possible to select input files based on the *contentFlag* or
*subType* but note the handling of these two sorts of parameters is
different and potentially could be different again for different data
file classes. The input file selection is specified by the qualifiers:
*requiredContentFlag* and *requiredSubType* both of which are expected
to be a list of integers.

Currently the *requiredContentFlag* is only appropriate for mini-MTZ
classes and should be set if a program can only handle a limited set of
the possible mini-MTZ representations and the GUI will then only allow
selection of files containing the appropriate representation or
containing data that can be converted to the appropriate representation.

The definition of *subType* is less clear-cut than the *contentFlag*
especially for PDB files so when specifying a list of *requiredSubType*
consider giving *CPdbDataFile.SUBTYPE_UNKNOWN* (i.e. 0) as the last
value in the list as this serves as a 'wildcard' which will allow any
PDB file to be selected but the GUI will set the default and the first
choices on a drop-down list to be of the preferred *subType*.

In conclusion: in *CPluginScript.processOutputFiles()* set *contentFlag*
and *subType* for the output file data classes that support these and in
the def file set qualifers *requiredContentFlag* and *requiredSubType*
for input file data classes.

--------------

Last modified: Thu Apr 23 10:21:59 BST 2015
