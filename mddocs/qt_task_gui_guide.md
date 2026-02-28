# CCP4i2 Qt-Based Task GUI — Comprehensive Guide

> **Purpose:** This document describes the architecture, conventions, and patterns used in CCP4i2's
> Python/Qt task GUIs. It is intended as a reference for automated translation of these interfaces
> into React-based equivalents on the Django branch.

---

## Table of Contents

1. [Architecture Overview](#1-architecture-overview)
2. [File Locations and Naming](#2-file-locations-and-naming)
3. [Class Structure and Metadata](#3-class-structure-and-metadata)
4. [The drawContents Method — Layout DSL](#4-the-drawcontents-method--layout-dsl)
5. [Folders and Tabs](#5-folders-and-tabs)
6. [Lines — The Fundamental Unit](#6-lines--the-fundamental-unit)
7. [SubFrames — Grouping Content](#7-subframes--grouping-content)
8. [Stacks — Conditional Display of Alternatives](#8-stacks--conditional-display-of-alternatives)
9. [Toggle System — Conditional Visibility](#9-toggle-system--conditional-visibility)
10. [Data Model and Parameters](#10-data-model-and-parameters)
11. [Widget Types and Data Type Mapping](#11-widget-types-and-data-type-mapping)
12. [Widget Qualifiers](#12-widget-qualifiers)
13. [Signal/Slot Reactive Patterns](#13-signalslot-reactive-patterns)
14. [Validation](#14-validation)
15. [The fix() Method — Pre-submission Cleanup](#15-the-fix-method--pre-submission-cleanup)
16. [Auto-generated GUIs](#16-auto-generated-guis)
17. [Advanced Patterns](#17-advanced-patterns)
18. [Complete Inventory of GUI Files](#18-complete-inventory-of-gui-files)
19. [Translation Notes for React](#19-translation-notes-for-react)

---

## 1. Architecture Overview

CCP4i2 task GUIs are Python classes that use a declarative DSL built on top of Qt (via PyQt/PySide).
Each task (wrapper or pipeline) can have a GUI class that defines how the user interacts with the
task's parameters before running it.

The key architectural layers:

```
┌─────────────────────────────────────────────┐
│              Task GUI Class                 │
│  (e.g., prosmart_refmac_gui.py)            │
│  Defines layout via drawContents()          │
├─────────────────────────────────────────────┤
│         CTaskWidget Base Class              │
│  (qtgui/CCP4TaskWidget.py)                 │
│  Provides openFolder, createLine, etc.      │
├─────────────────────────────────────────────┤
│         Widget Library                      │
│  (CCP4Widgets, CCP4XtalWidgets,            │
│   CCP4ModelWidgets, CCP4RefmacWidgets)     │
│  Maps data types → GUI widgets              │
├─────────────────────────────────────────────┤
│         Data Model (CContainer)             │
│  Defined in def.xml files                   │
│  Holds inputData, controlParameters,        │
│  outputData as typed CData objects          │
├─────────────────────────────────────────────┤
│         CCP4DataManager                     │
│  Registers data type → widget class map     │
│  Factory for creating widgets from models   │
└─────────────────────────────────────────────┘
```

**Data flow:**
1. A `def.xml` file declares all parameters with types and qualifiers
2. A `CContainer` is instantiated from the def.xml
3. The GUI class's `drawContents()` creates the visual layout, referencing parameters by name
4. Widgets are automatically created for each parameter based on its data type
5. Bidirectional binding keeps the GUI and container in sync
6. On "Run", the container is validated and serialized to XML for the wrapper/pipeline script

---

## 2. File Locations and Naming

GUI files live in the `script/` subdirectory of each task:

```
pipelines/<task_name>/script/<task_name>_gui.py
wrappers/<task_name>/script/<task_name>_gui.py
```

Corresponding parameter definition files:

```
pipelines/<task_name>/script/<task_name>.def.xml
wrappers/<task_name>/script/<task_name>.def.xml
```

Some tasks have additional files:
- `<task_name>_report.py` — result report generator
- `<task_name>.py` — the wrapper/pipeline script itself

---

## 3. Class Structure and Metadata

Every task GUI is a subclass of `CCP4TaskWidget.CTaskWidget`:

```python
from qtgui import CCP4TaskWidget

class prosmart_refmac_gui(CCP4TaskWidget.CTaskWidget):

    # Required class attributes
    TASKNAME = 'prosmart_refmac'          # Must match wrapper/pipeline name
    TASKVERSION = 0.1                      # Version number
    TASKMODULE = ['refinement']            # Menu category (list or string)
    TASKTITLE = 'Refinement with ProSMART' # Full title for task menu
    SHORTTASKTITLE = 'Refmac/ProSMART'     # Short title for job list
    DESCRIPTION = '''Refinement with ...'''# Tooltip description

    # Optional class attributes
    MGDISPLAYFILES = ['XYZOUT']           # Files to display in molecular graphics
    WHATNEXT = ['coot_rebuild', 'refmac']  # Suggested follow-up tasks
    PROGRAMHELP = ['refmac5']              # Help file references
    CLONEABLE = True                       # Can the job be cloned (default True)
    EDITOR = False                         # Is this a data editor, not a runner
    AUTOPOPULATEINPUT = True               # Auto-populate from previous job
    GUINAME = None                         # If GUI name differs from TASKNAME

    def __init__(self, parent):
        CCP4TaskWidget.CTaskWidget.__init__(self, parent)

    def drawContents(self):
        # Layout definition goes here
        pass
```

### TASKMODULE Values

Tasks are organized into menu categories. Valid modules (from `CCP4TaskManager.py`):

| Module | Description |
|--------|-------------|
| `data_entry` | Data import/entry |
| `data_reduction` | Data processing |
| `experimental_phasing` | Experimental phasing (SAD/MAD/etc.) |
| `molecular_replacement` | Molecular replacement |
| `density_modification` | Density modification |
| `model_building` | Model building |
| `refinement` | Refinement |
| `validation` | Validation |
| `bioinformatics` | Bioinformatics/sequence analysis |
| `ligand` | Ligand handling |
| `export` | Data export |
| `wrapper` | Developer wrappers (hidden in release) |
| `test` | Test tasks (hidden in release) |

---

## 4. The drawContents Method — Layout DSL

`drawContents()` is the core method where the visual layout is defined using a set of
declarative method calls that form a domain-specific language:

```python
def drawContents(self):
    # First folder — always "Input Data"
    self.openFolder(folderFunction='inputData', title='Input Data')
    self.createLine(['subtitle', 'Coordinates'])
    self.openSubFrame(frame=True)
    self.createLine(['widget', 'XYZIN'])
    self.closeSubFrame()
    self.createLine(['subtitle', 'Reflection data'])
    self.openSubFrame(frame=True)
    self.createLine(['widget', 'F_SIGF'])
    self.createLine(['widget', 'FREERFLAG'])
    self.closeSubFrame()
    self.closeFolder()

    # Second folder — options
    self.openFolder(title='Options')
    self.createLine(['widget', 'NCYCLES', 'label', 'Number of refinement cycles'])
    self.createLine(['widget', 'WEIGHT_OPT', 'label', 'Weighting'])
    self.closeFolder()

    # Third folder — deferred drawing for performance
    self.openFolder(title='Advanced', drawFolder=self.drawAdvanced)
    self.closeFolder()

def drawAdvanced(self):
    self.createLine(['subtitle', 'Advanced options'])
    self.createLine(['widget', 'BREFTYPE', 'label', 'B-factor model'])
```

### Method Summary

| Method | Purpose |
|--------|---------|
| `openFolder(...)` | Start a folder/tab section |
| `closeFolder()` | End a folder/tab section |
| `createLine([...])` | Create a horizontal row of widgets/labels |
| `openSubFrame(...)` | Start a grouped/indented section |
| `closeSubFrame()` | End a grouped section |
| `openStack(controlVar=...)` | Start mutually-exclusive content stack |
| `closeStack()` | End content stack |
| `autoGenerate(...)` | Auto-generate widgets from container |
| `setProgramHelpFile(...)` | Set help file for subsequent widgets |
| `createRadioGroup(...)` | Create a vertical radio button group |

---

## 5. Folders and Tabs

The top-level organizational unit. The GUI is split into folders (collapsible sections) or tabs,
depending on user preference. This is transparent to the developer.

```python
self.openFolder(
    folderFunction='inputData',    # Role: 'inputData', 'controlParameters', 'outputData', or omit
    title='Input Data',            # Display title
    toggle=[...],                  # Optional: conditional visibility
    toggleFunction=[...],          # Optional: function-based visibility
    drawFolder=self.drawMethod,    # Optional: deferred drawing for performance
    followFrom=True                # Optional: enable "follow from" feature (default True)
)
# ... folder contents ...
self.closeFolder()
```

### Folder Functions

| folderFunction | Description |
|----------------|-------------|
| `'inputData'` | First folder — required input files. System adds job title and "follow from" widgets. |
| `'controlParameters'` | Options/settings folder |
| `'outputData'` | Output configuration (rarely used — outputs are usually automatic) |
| `'general'` | Generic folder (default) |

### Deferred Drawing (Performance Optimization)

For complex tasks, only the first folder is drawn immediately. Other folders are drawn on demand:

```python
self.openFolder(title='Advanced Options', drawFolder=self.drawAdvanced)
# No content here — drawAdvanced() is called when user opens this folder
self.closeFolder()

def drawAdvanced(self):
    self.createLine(['subtitle', 'Advanced settings'])
    self.createLine(['widget', 'OBSCURE_PARAM'])
```

---

## 6. Lines — The Fundamental Unit

`createLine()` creates a single horizontal row. Its argument is a list of keywords and values
that are parsed left-to-right:

```python
self.createLine(definition, toggle=[], toggleFunction=[], appendLine=None)
```

### Definition Keywords

| Keyword | Arguments | Description |
|---------|-----------|-------------|
| `'label'` | text | Static text label (supports HTML) |
| `'widget'` | parameter_name | Data-bound widget for the named parameter |
| `'subtitle'` | text, tooltip | Section header with tooltip |
| `'advice'` | text | Italic informational text (should be sole content) |
| `'warning'` | text | Bold warning text |
| `'tip'` / `'message'` | text | Tooltip applied to subsequent widgets |
| `'help'` | target | Help link target for subsequent widgets |
| `'spacing'` | pixels | Fixed-width spacer |
| `'stretch'` | factor | Flexible spacer (usually at end of line) |
| `'launchButton'` | task_name | Button to launch a child task |

### Widget Qualifiers in Line Definitions

Qualifiers are prefixed with `-` and placed between `'widget'` and the parameter name:

```python
# Standard widget
self.createLine(['widget', 'XYZIN'])

# Widget with qualifiers
self.createLine(['widget', '-browseDb', True, 'XYZIN'])
self.createLine(['widget', '-guiMode', 'radio', 'METHOD'])
self.createLine(['widget', '-guiMode', 'multiLine', 'KEYWORDS'])
self.createLine(['widget', '-title', 'Select datasets', 'DATASET_LIST'])
self.createLine(['widget', '-guiLabel', 'Input model', 'XYZIN'])
self.createLine(['widget', '-charWidth', 50, 'TEXT_PARAM'])
self.createLine(['widget', '-jobCombo', True, 'HKLOUT'])
```

### Common Qualifier Summary

| Qualifier | Applies To | Description |
|-----------|-----------|-------------|
| `-browseDb` | CDataFile | Show database browser icon |
| `-guiMode` | CString/CInt/CFloat | Rendering mode: `'radio'`, `'combo'`, `'multiLine'`, `'edit'` |
| `-guiLabel` | CDataFile | Short label on the widget |
| `-title` | CList | Title text above list widget |
| `-charWidth` | CString | Width of text input in characters |
| `-jobCombo` | CDataFile | Show job output selector |
| `-enableEdit` | CDataFile | Enable/disable editing |

### Mixed Content Lines

A single line can contain multiple labels and widgets:

```python
self.createLine([
    'label', 'Number of cycles:',
    'widget', 'NCYCLES',
    'label', 'with weight',
    'widget', 'WEIGHT',
    'stretch', 1
])
```

### Appending to Lines

Multiple "lines" can be combined horizontally:

```python
line1 = self.createLine(['label', 'Method:', 'widget', 'METHOD'])
self.createLine(['widget', 'EXTRA_OPTION'], appendLine=line1)
```

---

## 7. SubFrames — Grouping Content

SubFrames group multiple lines together. They can optionally draw a border and can be toggled
visible/hidden:

```python
# Simple grouping with border
self.openSubFrame(frame=True)
self.createLine(['widget', 'PARAM1'])
self.createLine(['widget', 'PARAM2'])
self.closeSubFrame()

# Conditional subframe — only visible when USE_ADVANCED is True
self.openSubFrame(frame=True, toggle=['USE_ADVANCED', 'open', [True]])
self.createLine(['widget', 'ADVANCED_PARAM'])
self.closeSubFrame()

# Subframe with title and tooltip
self.openSubFrame(frame=True, title='NCS Options', tip='Non-crystallographic symmetry')
self.createLine(['widget', 'NCS_TYPE'])
self.closeSubFrame()
```

---

## 8. Stacks — Conditional Display of Alternatives

Stacks show one of several alternative content blocks based on a control variable (typically a
radio button or combo box). Only the selected alternative is visible at a time.

```python
# Control variable shown as radio buttons
self.createLine(['widget', '-guiMode', 'radio', 'INPUT_MODE'])

# Stack shows different content based on INPUT_MODE value
self.openStack(controlVar='INPUT_MODE')

# First alternative (shown when INPUT_MODE matches first enumerator)
self.createLine(['widget', 'FILE_INPUT'])

# Second alternative
self.createLine(['widget', 'SEQUENCE_INPUT'])

self.closeStack()
```

The control variable's `enumerators` qualifier determines which stack page maps to which value.

---

## 9. Toggle System — Conditional Visibility

The toggle system controls whether lines, subframes, or folders are visible. There are two
mechanisms:

### Simple Toggle (Value-Based)

```python
toggle = ['PARAMETER_NAME', 'open'|'closed', [value1, value2, ...]]
```

- **Parameter**: name of the controlling parameter
- **State**: `'open'` means visible when values match; `'closed'` means hidden when values match
- **Values**: list of values to test against (default `[True]`)

```python
# Show line only when USE_TLS is True
self.createLine(['widget', 'TLS_FILE'], toggle=['USE_TLS', 'open', [True]])

# Show subframe only when METHOD is 'SAD' or 'MAD'
self.openSubFrame(toggle=['METHOD', 'open', ['SAD', 'MAD']])

# Hide folder when SIMPLE_MODE is True
self.openFolder(title='Advanced', toggle=['SIMPLE_MODE', 'closed', [True]])
```

### Function Toggle (Complex Logic)

```python
toggleFunction = [function, [dependency1, dependency2, ...]]
```

- **Function**: a method returning `True` if the element should be visible
- **Dependencies**: parameter names that trigger re-evaluation when changed

```python
self.openSubFrame(
    toggleFunction=[self.showRosettaOptions, ['USE_EXISTING_MODELS', 'MODEL_TYPE']]
)

def showRosettaOptions(self):
    return (str(self.container.controlParameters.USE_EXISTING_MODELS) == 'False' and
            str(self.container.controlParameters.MODEL_TYPE) == 'ROSETTA')
```

### How Toggles Work Internally

1. A `CToggle` object is created for each toggle
2. It connects to the controlling parameter's `dataChanged` signal
3. When the parameter changes, the toggle evaluates and shows/hides the target widget
4. After initial drawing, `applyToggles()` is called to set correct initial visibility

---

## 10. Data Model and Parameters

### def.xml Structure

Each task's parameters are declared in a `def.xml` file with this structure:

```xml
<ccp4i2>
  <ccp4i2_header>
    <function>DEF</function>
    <pluginName>task_name</pluginName>
    <pluginTitle>Human-readable Title</pluginTitle>
  </ccp4i2_header>

  <ccp4i2_body id="task_name">

    <container id="inputData">
      <content id="XYZIN">
        <className>CPdbDataFile</className>
        <qualifiers>
          <mustExist>True</mustExist>
          <allowUndefined>False</allowUndefined>
          <fromPreviousJob>True</fromPreviousJob>
          <ifAtomSelection>True</ifAtomSelection>
        </qualifiers>
      </content>

      <content id="F_SIGF">
        <className>CObsDataFile</className>
        <qualifiers>
          <mustExist>True</mustExist>
          <allowUndefined>False</allowUndefined>
          <fromPreviousJob>True</fromPreviousJob>
        </qualifiers>
      </content>
    </container>

    <container id="controlParameters">
      <content id="NCYCLES">
        <className>CInt</className>
        <qualifiers>
          <default>5</default>
          <min>1</min>
          <max>200</max>
          <guiLabel>Number of cycles</guiLabel>
        </qualifiers>
      </content>

      <content id="BREFTYPE">
        <className>CString</className>
        <qualifiers>
          <default>ISOT</default>
          <enumerators>ISOT,ANIS,OVER,MIXE</enumerators>
          <menuText>Isotropic,Anisotropic,Overall,Mixed</menuText>
          <onlyEnumerators>True</onlyEnumerators>
          <guiLabel>B-factor refinement type</guiLabel>
        </qualifiers>
      </content>

      <content id="USE_TLS">
        <className>CBoolean</className>
        <qualifiers>
          <default>False</default>
          <guiLabel>Use TLS refinement</guiLabel>
        </qualifiers>
      </content>
    </container>

    <container id="outputData">
      <content id="XYZOUT">
        <className>CPdbDataFile</className>
        <qualifiers>
          <default>XYZOUT.pdb</default>
          <saveToDb>True</saveToDb>
        </qualifiers>
      </content>
    </container>

  </ccp4i2_body>
</ccp4i2>
```

### Container Hierarchy

```
CContainer (root)
├── inputData: CContainer
│   ├── XYZIN: CPdbDataFile
│   ├── F_SIGF: CObsDataFile
│   ├── FREERFLAG: CFreeRDataFile
│   └── ...
├── controlParameters: CContainer
│   ├── NCYCLES: CInt
│   ├── BREFTYPE: CString (enum)
│   ├── USE_TLS: CBoolean
│   └── ...
├── outputData: CContainer
│   ├── XYZOUT: CPdbDataFile
│   └── ...
└── guiAdmin: CContainer (system-managed)
    ├── jobTitle: CString
    ├── followFrom: CFollowFromJob
    └── patchSelection: CPatchSelection
```

### Accessing Parameters in Code

```python
# Through the container property
self.container.inputData.XYZIN                    # The CData model object
self.container.inputData.XYZIN.isSet()            # Check if value is set
self.container.inputData.XYZIN.exists()           # Check if file exists
self.container.controlParameters.NCYCLES          # Integer parameter
bool(self.container.controlParameters.USE_TLS)    # Boolean value
str(self.container.controlParameters.BREFTYPE)    # String value

# Setting values
self.container.controlParameters.NCYCLES.set(10)
self.container.controlParameters.BREFTYPE.set('ANIS')

# Getting/setting qualifiers
self.container.inputData.XYZIN.setQualifiers({'allowUndefined': True})
```

---

## 11. Widget Types and Data Type Mapping

The system automatically selects the appropriate widget class for each data type. The
`CCP4DataManager` maintains a registry mapping data type classes to widget classes.

### Core Data Types → Widgets

| Data Type | Python Class | Default Widget | Rendering |
|-----------|-------------|----------------|-----------|
| String | `CString` | `CStringView` | Text input, or combo/radio if enumerators defined |
| Integer | `CInt` | `CIntView` | Numeric input with min/max validation |
| Float | `CFloat` | `CFloatView` | Numeric input with min/max validation |
| Boolean | `CBoolean` | `CBooleanView` | Checkbox |
| Range | `CRange` | `CRangeView` | Min-max pair input |
| List | `CList` | `CListView` | Add/edit/delete list with item widgets |

### File Data Types → Widgets

| Data Type | Python Class | Widget | Features |
|-----------|-------------|--------|----------|
| Generic file | `CDataFile` | `CDataFileView` | File browser + optional DB browser |
| PDB model | `CPdbDataFile` | `CPdbDataFileView` | + atom selection interface |
| Sequence | `CSeqDataFile` | `CSeqDataFileView` | + inline sequence editor |
| Dictionary | `CDictDataFile` | `CDictDataFileView` | Restraint dictionary file |
| ASU content | `CAsuDataFile` | `CAsuDataFileView` | ASU composition |

### Crystallographic Data Types → Widgets

| Data Type | Python Class | Widget | Features |
|-----------|-------------|--------|----------|
| Obs reflections | `CObsDataFile` | `CMtzDataFileView` | MTZ with F/I column selection |
| Phases | `CPhsDataFile` | `CMtzDataFileView` | MTZ with phase column selection |
| Map coefficients | `CMapCoeffsDataFile` | `CMtzDataFileView` | MTZ with map coeff columns |
| Free R | `CFreeRDataFile` | `CMtzDataFileView` | MTZ with Free R column |
| Space group | `CSpaceGroup` | `CSpaceGroupView` | Space group selector |
| Unit cell | `CCell` | `CCellView` | Cell parameter inputs |
| Space group + cell | `CSpaceGroupCell` | `CSpaceGroupCellView` | Combined selector |

### Rendering Modes for Enumerated Types

When a `CString`, `CInt`, or `CFloat` has `enumerators` defined, the widget rendering depends
on the `guiMode` qualifier:

| guiMode | Rendering | When to Use |
|---------|-----------|-------------|
| `'combo'` (default) | Drop-down select | 4+ options, or space-constrained |
| `'radio'` | Radio button group | 2-4 mutually exclusive options |
| `'multiLine'` | Multi-line text editor | Free-form text entry |
| `'edit'` | Single-line text input | Default for strings without enumerators |

```python
# Combo box (default when enumerators exist)
self.createLine(['widget', 'BREFTYPE'])

# Explicit radio buttons
self.createLine(['widget', '-guiMode', 'radio', 'BREFTYPE'])

# Multi-line text editor
self.createLine(['widget', '-guiMode', 'multiLine', 'KEYWORDS'])
```

---

## 12. Widget Qualifiers

Qualifiers are metadata attributes on data objects that control both validation and widget behavior.
They can be set in `def.xml` or dynamically in code.

### Common Qualifiers (All Types)

| Qualifier | Type | Description |
|-----------|------|-------------|
| `default` | varies | Default value |
| `allowUndefined` | bool | Whether None/unset is acceptable |
| `guiLabel` | str | Display label for the parameter |
| `toolTip` | str | Tooltip help text |
| `guiMode` | str | Widget rendering mode |
| `guiDefinition` | dict | Advanced GUI configuration |

### String Qualifiers

| Qualifier | Type | Description |
|-----------|------|-------------|
| `enumerators` | list | Allowed/suggested values |
| `menuText` | list | Display labels for enumerators |
| `onlyEnumerators` | bool | Restrict to enumerated values only |
| `minLength` | int | Minimum string length |
| `maxLength` | int | Maximum string length |
| `charWidth` | int | Widget width in characters |

### Numeric Qualifiers (CInt, CFloat)

| Qualifier | Type | Description |
|-----------|------|-------------|
| `min` | number | Minimum allowed value |
| `max` | number | Maximum allowed value |
| `enumerators` | list | Predefined value options |

### File Qualifiers

| Qualifier | Type | Description |
|-----------|------|-------------|
| `mustExist` | bool | File must exist on disk |
| `fromPreviousJob` | bool | Can select from previous job outputs |
| `browseDb` | bool | Show database browser icon |
| `browseFiles` | bool | Show file browser |
| `ifAtomSelection` | bool | Show atom selection interface |
| `requiredContentFlag` | list | Required data content types |
| `requiredSubType` | list | Required file subtypes |

### List Qualifiers

| Qualifier | Type | Description |
|-----------|------|-------------|
| `listMinLength` | int | Minimum number of items |
| `listMaxLength` | int | Maximum number of items |

### Dynamic Qualifier Changes

Qualifiers can be changed at runtime to modify validation requirements:

```python
@QtCore.Slot()
def handleModeChange(self):
    if bool(self.container.controlParameters.USE_TLS):
        self.container.inputData.TLS_FILE.setQualifiers({
            'allowUndefined': False,
            'mustExist': True
        })
    else:
        self.container.inputData.TLS_FILE.setQualifiers({
            'allowUndefined': True
        })
    self.getWidget('TLS_FILE').validate()
```

---

## 13. Signal/Slot Reactive Patterns

Task GUIs use Qt signals and slots to react to user input changes. The primary signal is
`dataChanged`, emitted by every data object when its value changes.

### Connecting to Data Changes

Two equivalent approaches:

```python
# Approach 1: Direct Qt connection
from PySide2 import QtCore
self.container.inputData.XYZIN.dataChanged.connect(self.handleXyzinChange)

# Approach 2: Convenience method (simpler, dataChanged only)
self.connectDataChanged('XYZIN', self.handleXyzinChange)
```

### Typical Callback Patterns

**1. Update dependent requirements:**

```python
@QtCore.Slot()
def handleModeChange(self):
    mode = str(self.container.controlParameters.MODE)
    if mode == 'SAD':
        self.container.inputData.SAD_DATA.setQualifiers({'allowUndefined': False})
        self.container.inputData.NATIVE_DATA.setQualifiers({'allowUndefined': True})
    else:
        self.container.inputData.SAD_DATA.setQualifiers({'allowUndefined': True})
        self.container.inputData.NATIVE_DATA.setQualifiers({'allowUndefined': False})
    self.getWidget('SAD_DATA').validate()
    self.getWidget('NATIVE_DATA').validate()
```

**2. Set defaults from input data:**

```python
@QtCore.Slot()
def handleModelChange(self):
    if not self.container.inputData.XYZIN.isSet():
        return
    if not self.container.inputData.XYZIN.exists():
        return
    # Analyze the model file and set defaults
    chains = self.getChainList()
    self.container.controlParameters.CHAIN_LIST.set(chains)
```

**3. Cross-validate inputs:**

```python
@QtCore.Slot()
def handleReflectionChange(self):
    if self.container.inputData.F_SIGF.isSet():
        wavelength = self.container.inputData.F_SIGF.getWavelength()
        self.container.controlParameters.WAVELENGTH.set(wavelength)
```

### Important: Always Initialize

Signal connections are typically made at the end of `drawContents()`, and the handler is called
once immediately to set the correct initial state:

```python
def drawContents(self):
    # ... layout code ...

    # Connect and initialize
    self.connectDataChanged('XYZIN', self.handleModelChange)
    self.handleModelChange()  # Set initial state
```

### editingFinished Signal

For text input widgets, `dataChanged` fires on every keystroke. Use `editingFinished` instead
to wait for the user to finish typing:

```python
self.getWidget('TEXT_PARAM').widget.editingFinished.connect(self.handleTextComplete)
```

---

## 14. Validation

### Automatic Validation

Each widget validates its value based on qualifiers:
- File widgets check `mustExist`, `allowUndefined`
- Numeric widgets check `min`, `max`
- String widgets check `minLength`, `maxLength`, `onlyEnumerators`
- List widgets check `listMinLength`, `listMaxLength`

Invalid widgets are highlighted in red.

### Custom Validation: isValid()

Override `isValid()` to add cross-parameter validation:

```python
def isValid(self):
    invalidElements = super().isValid()

    # Cross-check: space groups must match
    if (self.container.inputData.XYZIN.isSet() and
        self.container.inputData.F_SIGF.isSet()):
        sg_model = self.container.inputData.XYZIN.getSpaceGroup()
        sg_data = self.container.inputData.F_SIGF.getSpaceGroup()
        if sg_model != sg_data:
            # Add to invalid list or show warning
            pass

    return invalidElements
```

### Custom Validation: taskValidation()

Override `taskValidation()` for validation that returns an error report:

```python
def taskValidation(self):
    report = CErrorReport()
    if some_condition_fails:
        report.append(self.__class__, ERROR_CODE, details='...')
    return report
```

### Triggering Widget Re-validation

After changing qualifiers, explicitly re-validate the affected widget:

```python
self.container.inputData.PARAM.setQualifiers({'allowUndefined': False})
self.getWidget('PARAM').validate()
```

---

## 15. The fix() Method — Pre-submission Cleanup

`fix()` is called just before job submission. It allows the GUI to clean up, normalize, or
compute derived values:

```python
def fix(self):
    # Unset unused parameters to keep the container clean
    if not bool(self.container.controlParameters.USE_TLS):
        self.container.inputData.TLS_FILE.unSet()

    # Compute derived values
    self.container.controlParameters.COMPUTED_VALUE.set(
        self.computeSomething()
    )

    # Check for errors
    if something_wrong:
        return CErrorReport(self.__class__, ERROR_CODE)

    return CErrorReport()  # Empty = success
```

Common patterns in `fix()`:
- Unset parameters that aren't relevant for the selected mode
- Set content flags on data files (`setContentFlag()`)
- Validate data types/formats that can't be checked by simple qualifiers
- Compute derived control parameters

---

## 16. Auto-generated GUIs

For simple parameters, the GUI can be auto-generated from the container definition:

### Basic Auto-generation

```python
self.openFolder(title='Options')
self.autoGenerate(
    container=self.container.controlParameters,
    selection={'includeParameters': ['PARAM1', 'PARAM2', 'PARAM3']}
)
self.closeFolder()
```

### Selection Options

```python
# Include specific parameters
selection={'includeParameters': ['PARAM1', 'PARAM2']}

# Exclude specific parameters (include everything else)
selection={'excludeParameters': ['INTERNAL_PARAM']}

# Filter by guiDefinition keys
selection={'keyValues': {'expertLevel': '0'}}
```

### Nested Auto-generation

Some tasks implement recursive auto-generation for deeply nested containers:

```python
def nestedAutoGenerate(self, container, expertLevel=['0']):
    dataOrder = container.dataOrder()
    contents = [getattr(container, name) for name in dataOrder]

    # Filter by expertLevel
    filtered = [name for name, c in zip(dataOrder, contents)
                if c.qualifiers('guiDefinition').get('expertLevel', '0') in expertLevel]

    # Separate regular data from sub-containers
    data_names = [n for n in filtered if not isinstance(getattr(container, n), CContainer)]
    sub_containers = [getattr(container, n) for n in filtered
                      if isinstance(getattr(container, n), CContainer)]

    if data_names:
        self.autoGenerate(container, selection={'includeParameters': data_names})

    for sub in sub_containers:
        self.openSubFrame(frame=True)
        self.nestedAutoGenerate(sub, expertLevel=expertLevel)
        self.closeSubFrame()
```

### guiDefinition Qualifiers for Auto-generation

In `def.xml`, parameters can have `guiDefinition` qualifiers that control auto-generated layout:

```xml
<content id="OBSCURE_PARAM">
  <className>CFloat</className>
  <qualifiers>
    <guiLabel>Obscure threshold</guiLabel>
    <toolTip>Threshold for obscure calculation</toolTip>
    <guiDefinition>
      <expertLevel>1</expertLevel>
      <toggleParameter>USE_OBSCURE</toggleParameter>
      <toggleState>open</toggleState>
      <toggleValues>True</toggleValues>
    </guiDefinition>
  </qualifiers>
</content>
```

---

## 17. Advanced Patterns

### 17.1 Multiple GUI Interfaces for One Wrapper

A single wrapper can have multiple GUI interfaces (e.g., different modes of the same program):

```python
class molrep_den_gui(CTaskWidget):
    TASKNAME = 'molrep_mr'       # The wrapper to run
    GUINAME = 'molrep_den'        # Unique GUI identifier
    TASKMODULE = 'model_building'
```

The derived def.xml includes the base and overrides specific parameters:

```xml
<ccp4i2_body id="molrep_den">
  <file>
    <CI2XmlDataFile>
      <project>CCP4I2_TOP</project>
      <relPath>wrappers/molrep_mr/script</relPath>
      <baseName>molrep_mr.def.xml</baseName>
    </CI2XmlDataFile>
  </file>
  <!-- Override specific parameters here -->
</ccp4i2_body>
```

### 17.2 GUI Inheritance Between Tasks

GUI classes can inherit from other GUI classes:

```python
from pipelines.phaser_pipeline.wrappers.phaser_MR_AUTO.script import phaser_MR_AUTO_gui

class phaser_pipeline_gui(phaser_MR_AUTO_gui.phaser_MR_AUTO_gui):
    TASKNAME = 'phaser_pipeline'
    TASKTITLE = 'Phaser Pipeline'

    def drawContents(self):
        self.drawPhaserFrontPage()          # Inherited method
        self.openFolder(title='Extra Steps', drawFolder=self.drawExtraSteps)
        self.closeFolder()

    def drawExtraSteps(self):
        # Additional content specific to the pipeline
        self.createLine(['widget', 'RUNREFMAC', 'label', 'Run REFMAC after'])
```

### 17.3 Child Pop-out Tasks

A task can launch another task in a pop-out window:

```python
self.createLine(['launchButton', 'gesamt', 'label', 'Superpose using Gesamt'])

def handleLaunchedJob(self, jobId=None, status=None, taskWidget=None):
    if status == 1 and taskWidget is not None:
        # Child task just opened — set its inputs from our data
        taskWidget.container.inputData.XYZIN_QUERY.set(
            self.container.inputData.XYZIN
        )
    elif status == 6:
        # Child task finished — use its output
        from core import CCP4Modules
        fileList = CCP4Modules.PROJECTSMANAGER().db().getJobFilesInfo(
            jobId=jobId, jobParamName='XYZOUT'
        )
        if fileList:
            self.getWidget('XYZIN').model.setDbFileId(fileList[0]['fileId'])
```

### 17.4 Showing Atom Selection

For PDB file widgets, atom selection can be enabled:

```python
self.createLine(['widget', 'XYZIN'])
self.getWidget('XYZIN').showAtomSelection()
```

### 17.5 Platform-Specific UI

Some tasks adjust their UI based on the operating system:

```python
import platform

def drawContents(self):
    if platform.system() == 'Windows':
        self.createLine(['widget', 'IMAGE_DIRECTORY'])
    else:
        self.createLine(['widget', 'IMAGE_FILE'])
        self.createLine(['advice', '...or specify a directory'])
        self.createLine(['widget', 'IMAGE_DIRECTORY'])
```

### 17.6 Dynamic Messages

The task widget has a message display area for status/warning messages:

```python
self.setMessage('Warning: space groups do not match', parameter='SG_WARNING')
self.unsetMessage(parameter='SG_WARNING')
```

### 17.7 Running Analysis from the GUI

GUIs can run wrapper scripts for analysis:

```python
from core import CCP4Modules, CCP4PluginScript

workDir = CCP4Modules.PROJECTSMANAGER().jobDirectory(self.jobId(), subDir='TMP')
plugin = CCP4PluginScript.CPluginScript(
    dummy=True, workDirectory=workDir
).makePluginObject(pluginName='analysis_task', reportToDatabase=False)
result = plugin.process(self.container)
```

---

## 18. Complete Inventory of GUI Files

### Pipelines (22 files)

| File | Task | Module |
|------|------|--------|
| `pipelines/boilerplate/script/boilerplate_gui.py` | Template | — |
| `pipelines/LidiaAcedrg/script/lidiaAcedrg_gui.py` | Lidia+AceDRG | ligand |
| `pipelines/LidiaAcedrgNew/script/lidiaAcedrgNew_gui.py` | Lidia+AceDRG (new) | ligand |
| `pipelines/MakeLink/script/MakeLink_gui.py` | Make link | ligand |
| `pipelines/MakeProjectsAndDoLigandPipeline/script/MakeProjectsAndDoLigandPipeline_gui.py` | Ligand pipeline | ligand |
| `pipelines/PrepareDeposit/script/PrepareDeposit_gui.py` | Deposition prep | export |
| `pipelines/SubstituteLigand/script/SubstituteLigand_gui.py` | Substitute ligand | ligand |
| `pipelines/dr_mr_modelbuild_pipeline/script/dr_mr_modelbuild_pipeline_gui.py` | DR+MR+build | molecular_replacement |
| `pipelines/dr_mr_modelbuild_pipeline/wrappers/mrbump_model_prep/script/mrbump_model_prep_gui.py` | MrBUMP prep | molecular_replacement |
| `pipelines/molrep_pipe/script/molrep_gui.py` | Molrep pipeline | molecular_replacement |
| `pipelines/nautilus_build_refine/script/nautilus_build_refine_gui.py` | Nautilus build+refine | model_building |
| `pipelines/phaser_ep/script/phaser_EP_gui.py` | Phaser EP | experimental_phasing |
| `pipelines/phaser_pipeline/script/phaser_pipeline_gui.py` | Phaser pipeline | molecular_replacement |
| `pipelines/phaser_pipeline/wrappers/phaser_EP_AUTO/script/phaser_EP_AUTO_gui.py` | Phaser EP AUTO | experimental_phasing |
| `pipelines/phaser_pipeline/wrappers/phaser_EP_LLG/script/phaser_EP_LLG_gui.py` | Phaser EP LLG | experimental_phasing |
| `pipelines/phaser_pipeline/wrappers/phaser_MR_AUTO/script/phaser_MR_AUTO_gui.py` | Phaser MR AUTO | molecular_replacement |
| `pipelines/phaser_pipeline/wrappers/phaser_MR_RNP/script/phaser_MR_RNP_gui.py` | Phaser MR RNP | molecular_replacement |
| `pipelines/phaser_rnp_pipeline/script/phaser_rnp_pipeline_gui.py` | Phaser RNP pipeline | molecular_replacement |
| `pipelines/phaser_simple/script/phaser_simple_gui.py` | Phaser simple | molecular_replacement |
| `pipelines/pisapipe/script/pisapipe_gui.py` | PISA pipeline | validation |
| `pipelines/prosmart_refmac/script/prosmart_refmac_gui.py` | ProSMART+Refmac | refinement |
| `pipelines/servalcat_pipe/script/servalcat_pipe_gui.py` | Servalcat pipeline | refinement |

### Wrappers (57 files)

| File | Task | Module |
|------|------|--------|
| `wrappers/AMPLE/script/AMPLE_gui.py` | AMPLE | molecular_replacement |
| `wrappers/AcedrgLink/script/AcedrgLink_gui.py` | AceDRG Link | ligand |
| `wrappers/AlternativeImportXIA2/script/AlternativeImportXIA2_gui.py` | Import XIA2 | data_entry |
| `wrappers/MakeMonster/script/MakeMonster_gui.py` | Make Monster | model_building |
| `wrappers/ProvideAlignment/script/ProvideAlignment_gui.py` | Provide alignment | bioinformatics |
| `wrappers/ProvideAsuContents/script/ProvideAsuContents_gui.py` | ASU contents | data_entry |
| `wrappers/ProvideSequence/script/ProvideSequence_gui.py` | Provide sequence | data_entry |
| `wrappers/SIMBAD/script/SIMBAD_gui.py` | SIMBAD | molecular_replacement |
| `wrappers/ShelxCDE/script/ShelxCD_gui.py` | SHELXC/D | experimental_phasing |
| `wrappers/ShelxCDE/script/ShelxCECompareHands_gui.py` | SHELXE compare | experimental_phasing |
| `wrappers/ShelxCDE/script/ShelxCE_gui.py` | SHELXC/E | experimental_phasing |
| `wrappers/SubtractNative/script/SubtractNative_gui.py` | Subtract native | data_reduction |
| `wrappers/TestObsConversions/script/TestObsConversions_gui.py` | Test conversions | test |
| `wrappers/add_fractional_coords/script/add_fractional_coords_gui.py` | Add frac coords | — |
| `wrappers/adding_stats_to_mmcif_i2/script/adding_stats_to_mmcif_i2_gui.py` | Add stats mmCIF | export |
| `wrappers/arcimboldo/script/arcimboldo_gui.py` | ARCIMBOLDO | molecular_replacement |
| `wrappers/boilerplate/script/boilerplate_gui.py` | Boilerplate | — |
| `wrappers/ccp4mg_edit_model/script/ccp4mg_edit_model_gui.py` | CCP4mg edit | model_building |
| `wrappers/ccp4mg_edit_nomrbump/script/ccp4mg_edit_nomrbump_gui.py` | CCP4mg edit | model_building |
| `wrappers/ccp4mg_general/script/ccp4mg_general_gui.py` | CCP4mg general | model_building |
| `wrappers/chltofom/script/chltofom_gui.py` | CHL to FOM | data_reduction |
| `wrappers/clustalw/script/clustalw_gui.py` | ClustalW | bioinformatics |
| `wrappers/cmapcoeff/script/cmapcoeff_gui.py` | Map coefficients | data_reduction |
| `wrappers/comit/script/comit_gui.py` | Composite omit | refinement |
| `wrappers/coordinate_selector/script/coordinate_selector_gui.py` | Coord selector | data_entry |
| `wrappers/coot1/script/coot1_gui.py` | Coot | model_building |
| `wrappers/coot_find_waters/script/coot_find_waters_gui.py` | Coot find waters | model_building |
| `wrappers/coot_rebuild/script/coot_rebuild_gui.py` | Coot rebuild | model_building |
| `wrappers/coot_rsr_morph/script/coot_rsr_morph_gui.py` | Coot RSR morph | model_building |
| `wrappers/coot_script_lines/script/coot_script_lines_gui.py` | Coot scripting | model_building |
| `wrappers/cpatterson/script/cpatterson_gui.py` | Patterson | experimental_phasing |
| `wrappers/csymmatch/script/csymmatch_gui.py` | Csymmatch | validation |
| `wrappers/ctruncate/script/ctruncate_gui.py` | Ctruncate | data_reduction |
| `wrappers/density_calculator/script/density_calculator_gui.py` | Density calc | — |
| `wrappers/editbfac/script/editbfac_gui.py` | Edit B-factors | model_building |
| `wrappers/i2Dimple/script/i2Dimple_gui.py` | DIMPLE | refinement |
| `wrappers/import_mosflm/script/import_mosflm_gui.py` | Import Mosflm | data_entry |
| `wrappers/lorestr_i2/script/lorestr_i2_gui.py` | LORESTR | refinement |
| `wrappers/mergeMtz/script/mergeMtz_gui.py` | Merge MTZ | data_reduction |
| `wrappers/metalCoord/script/metalCoord_gui.py` | Metal coordination | ligand |
| `wrappers/modelASUCheck/script/modelASUCheck_gui.py` | ASU check | validation |
| `wrappers/modelcraft/script/modelcraft_gui.py` | ModelCraft | model_building |
| `wrappers/morda_i2/script/morda_i2_gui.py` | MoRDa | molecular_replacement |
| `wrappers/mosflm/script/mosflm_gui.py` | Mosflm | data_reduction |
| `wrappers/pdb_redo_api/script/pdb_redo_api_gui.py` | PDB-REDO API | refinement |
| `wrappers/pdbview_edit/script/pdbview_edit_gui.py` | PDB view/edit | model_building |
| `wrappers/phaser_ensembler/script/phaser_ensembler_gui.py` | Phaser ensembler | molecular_replacement |
| `wrappers/phaser_phil/script/phaser_phil_gui.py` | Phaser PHIL | molecular_replacement |
| `wrappers/pointless_reindexToMatch/script/pointless_reindexToMatch_gui.py` | Pointless reindex | data_reduction |
| `wrappers/privateer/script/privateer_gui.py` | Privateer | validation |
| `wrappers/qtpisa/script/qtpisa_gui.py` | PISA | validation |
| `wrappers/sheetbend/script/sheetbend_gui.py` | Sheetbend | refinement |
| `wrappers/validate_protein/script/validate_protein_gui.py` | Validate protein | validation |
| `wrappers/xia2_dials/script/xia2_dials_gui.py` | XIA2/DIALS | data_reduction |
| `wrappers/xia2_multiplex/script/xia2_multiplex_gui.py` | XIA2 multiplex | data_reduction |
| `wrappers/xia2_ssx_reduce/script/xia2_ssx_reduce_gui.py` | XIA2 SSX reduce | data_reduction |
| `wrappers/xia2_xds/script/xia2_xds_gui.py` | XIA2/XDS | data_reduction |
| `wrappers/zanuda/script/zanuda_gui.py` | Zanuda | validation |

---

## 19. Translation Notes for React

This section provides guidance for mapping Qt GUI patterns to React components.

### 19.1 Structural Mapping

| Qt Concept | React Equivalent |
|------------|------------------|
| `CTaskWidget` class | React component (functional or class) |
| `drawContents()` | Component render/return JSX |
| `openFolder()` / `closeFolder()` | `<Folder>` or `<Tabs.Panel>` component |
| `createLine([...])` | `<FormRow>` component with child elements |
| `openSubFrame()` / `closeSubFrame()` | `<FieldGroup>` or `<Card>` component |
| `openStack()` / `closeStack()` | Conditional rendering or `<Switch>`/`<Tabs>` |
| `'subtitle'` keyword | `<SectionHeader>` component |
| `'advice'` keyword | `<HelpText>` or `<Alert variant="info">` |
| `'warning'` keyword | `<Alert variant="warning">` |
| `'label'` keyword | `<label>` or `<span>` |

### 19.2 Widget Mapping

| Qt Widget | React Component |
|-----------|----------------|
| `CStringView` (edit) | `<TextInput>` |
| `CStringView` (combo) | `<Select>` |
| `CStringView` (radio) | `<RadioGroup>` |
| `CStringView` (multiLine) | `<TextArea>` |
| `CIntView` | `<NumberInput type="integer">` |
| `CFloatView` | `<NumberInput type="float">` |
| `CBooleanView` | `<Checkbox>` |
| `CDataFileView` | `<FileSelector>` with optional DB browser |
| `CPdbDataFileView` | `<FileSelector>` with atom selection |
| `CMtzDataFileView` | `<FileSelector>` with column selection |
| `CSeqDataFileView` | `<FileSelector>` with inline editor |
| `CListView` | `<DynamicList>` with item add/edit/delete |
| `CSpaceGroupView` | `<SpaceGroupSelector>` |
| `CRangeView` | `<RangeInput>` (min/max pair) |

### 19.3 Data Binding

| Qt Pattern | React Equivalent |
|------------|------------------|
| `CContainer` model | Form state (React Hook Form, Formik, or custom context) |
| `dataChanged` signal | `onChange` callbacks or state subscriptions |
| `connectDataChanged()` | `useEffect` with dependency array |
| `model.set(value)` | State setter / form setValue |
| `model.isSet()` | Check for non-null/non-undefined |
| `model.get()` | Read from form state |
| `setQualifiers()` | Dynamic validation rule updates |
| `validate()` | Form validation (Yup, Zod, or custom) |

### 19.4 Toggle System

| Qt Pattern | React Equivalent |
|------------|------------------|
| `toggle=['PARAM', 'open', [values]]` | Conditional render: `{paramValue === value && <Component>}` |
| `toggleFunction=[fn, deps]` | `useMemo` or derived state for visibility |
| `openSubFrame(toggle=...)` | Wrapper component with visibility logic |

### 19.5 Key Behavioral Patterns to Preserve

1. **Initial State**: Signal handlers in Qt are called once after `drawContents()` to set initial state.
   In React, use `useEffect` with empty deps or initial state computation.

2. **Deferred Rendering**: Folders with `drawFolder` callbacks only render when opened.
   In React, use lazy rendering or `React.lazy` for tab/accordion panels.

3. **Cascading Validation**: When one parameter changes, it may alter qualifiers of other parameters.
   In React, use `useEffect` to watch dependencies and update validation schemas.

4. **Guard Clauses in Handlers**: Qt handlers always check `isSet()` and `exists()` before processing.
   React handlers should similarly check for null/undefined values.

5. **Pre-submission Cleanup**: The `fix()` method runs before submission.
   In React, implement as a form `onSubmit` pre-processor or middleware.

6. **Bidirectional Binding**: Qt widgets auto-sync with the model.
   React controlled components with state management achieve the same effect.

### 19.6 createLine Definition Parsing

The `createLine` definition list is a mini-DSL that should be parsed left-to-right into React JSX.
Each keyword consumes a fixed number of arguments:

```
keyword     args    produces
-------     ----    --------
'label'     1       <span>{text}</span>
'widget'    1+      <WidgetForParam name={paramName} qualifiers={...} />
'subtitle'  1-2     <SectionHeader title={text} tooltip={tip} />
'advice'    1       <HelpText>{text}</HelpText>
'warning'   1       <WarningText>{text}</WarningText>
'tip'       1       sets tooltip for next widget
'message'   1       same as 'tip'
'help'      1       sets help target for next widget
'spacing'   1       <Spacer width={px} />
'stretch'   0-1     <FlexSpacer />
'-qualifier' 1      adds qualifier to next 'widget'
'launchButton' 1    <LaunchTaskButton task={name} />
```

A parser for this list would iterate through the definition array, accumulating qualifiers and
tooltips until a `'widget'` keyword is encountered, at which point a widget component is emitted
with all accumulated context.

### 19.7 Complexity Spectrum

GUI files range from trivial to very complex. When planning automated translation, note:

- **Trivial** (< 50 lines): `mergeMtz_gui.py`, `coordinate_selector_gui.py`, `coot1_gui.py` —
  Just folder + widget declarations, no callbacks.

- **Simple** (50-100 lines): `comit_gui.py`, `cmapcoeff_gui.py`, `chltofom_gui.py` —
  Basic folders with some toggles.

- **Medium** (100-300 lines): `modelcraft_gui.py`, `ctruncate_gui.py`, `PrepareDeposit_gui.py` —
  Multiple folders, signal connections, dynamic qualifiers.

- **Complex** (300+ lines): `prosmart_refmac_gui.py` (882 lines), `AMPLE_gui.py`,
  `arcimboldo_gui.py`, `xia2_dials_gui.py` — Multiple deferred folders, many toggle functions,
  extensive signal handling, custom validation.

---

*This guide was generated from analysis of 79 task GUI files, the CCP4TaskWidget base class
(1858 lines), the task_guis.html developer documentation, and the CCP4 widget/data type libraries.*
