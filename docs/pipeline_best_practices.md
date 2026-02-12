# CCP4i2 Pipeline Development Best Practices

*Reference implementations: `SubstituteLigand` pipeline (basic pattern),
`servalcat_pipe` pipeline (real-time progress monitoring)*

## 1. Architecture Principles

### 1.1 Synchronous Execution Model

All pipeline phases execute **synchronously** via `plugin.process()`. This eliminates
race conditions, simplifies error propagation, and makes the control flow readable
top-to-bottom.

```python
# CORRECT: Synchronous, blocking
status = plugin.process()
if status != CPluginScript.SUCCEEDED:
    return error

# WRONG: Asynchronous with signal handlers
plugin.doAsync = True
plugin.connectSignal(plugin, 'finished', self.onFinished)
plugin.process()
# Control returns here before job finishes!
```

When a pipeline needs to call another pipeline that declares `ASYNCHRONOUS = True`,
force it synchronous:

```python
plugin.doAsync = False
status = plugin.process()
```

### 1.2 Override `startProcess()`, Not `process()`

The base class `process()` method orchestrates the full lifecycle:

```
process()
  1. validity()
  2. checkOutputData()
  3. processInputFiles()
  4. makeCommandAndScript()     # No-op for pipelines
  5. startProcess()             # <-- Override this
  6. processOutputFiles()
```

Override `startProcess()` for pipeline orchestration logic. Override
`processInputFiles()` for input validation and mode detection.
Override `processOutputFiles()` for output harvesting.

### 1.3 Clear Phase Separation

Each pipeline phase should be a private method with a consistent pattern:

```python
def startProcess(self):
    error = CErrorReport()

    # Phase 1: Optional step
    if self._needsPhase1:
        phase1_error = self._runPhase1()
        if phase1_error and phase1_error.maxSeverity() >= 4:
            return phase1_error

    # Phase 2: Required step
    phase2_error = self._runPhase2()
    if phase2_error and phase2_error.maxSeverity() >= 4:
        return phase2_error

    return error
```

## 2. Class Structure

### 2.1 Class-Level Constants

```python
class MyPipeline(CPluginScript):
    TASKMODULE = 'refinement'              # Module category
    TASKTITLE = 'My Pipeline'              # Display title
    TASKNAME = 'my_pipeline'               # Unique identifier
    MAINTAINER = 'dev@example.com'         # Contact
    WHATNEXT = ['coot_rebuild', 'my_pipeline']  # Suggested next tasks

    ERROR_CODES = {
        201: {'description': 'General pipeline failure'},
        202: {'description': 'Failed in harvest operation'},
        203: {'description': 'Failed in phase 1 (task_name)'},
        204: {'description': 'Failed in phase 2 (task_name)'},
        209: {'description': 'Failed to create sub-plugin'},
        211: {'description': 'Missing required output from sub-plugin'},
        212: {'description': 'Invalid input configuration'},
    }
```

### 2.2 Instance Initialization

```python
def __init__(self, *args, **kws):
    super(MyPipeline, self).__init__(*args, **kws)
    self.xmlroot = etree.Element('MyPipeline')

    # Intermediate state - populated during startProcess()
    self.intermediateFile = None

    # Sub-plugin references - for processOutputFiles() to access
    self.phase1Plugin = None
    self.phase2Plugin = None

    # Pipeline mode flags - determined in processInputFiles()
    self._mode = None
```

**Rules:**
- Initialize XML root element
- Declare all intermediate state variables as `None`
- Store sub-plugin references for later harvesting
- Declare mode flags for conditional execution

## 3. Validation

### 3.1 The `validity()` Override

Override `validity()` to filter validation errors for conditionally-required inputs:

```python
def validity(self):
    error = super(MyPipeline, self).validity()

    mode = str(self.container.controlParameters.MODE)
    if mode != 'ADVANCED':
        # Filter out errors for ADVANCED-only inputs
        filtered = CErrorReport()
        for err in error.getErrors():
            if err.get('name', '') == 'ADVANCED_INPUT':
                continue
            filtered.append(
                err.get('class', ''),
                err.get('code', 0),
                err.get('details', ''),
                err.get('name', ''),
                err.get('severity', 0)
            )
        return filtered
    return error
```

**Rules:**
- Always call `super().validity()` first
- Only filter; never add side effects or modify qualifiers
- Use the DEF XML `allowUndefined` qualifier for truly optional inputs

### 3.2 Input Validation in `processInputFiles()`

```python
def processInputFiles(self):
    error = CErrorReport()

    # Check required files exist
    invalidFiles = self.checkInputData()
    # Filter out conditionally-required files
    for f in list(invalidFiles):
        if str(self.container.controlParameters.MODE) == 'SIMPLE' and str(f) == 'ADVANCED_INPUT':
            invalidFiles.remove(f)

    if len(invalidFiles) > 0:
        for f in invalidFiles:
            self.appendErrorReport(212, f'Missing required input: {f}')
        error.append(self.__class__.__name__, 212,
                    f'Missing required inputs: {invalidFiles}',
                    'processInputFiles', 4)
        return error

    # Determine pipeline modes from parameters
    self._mode = str(self.container.controlParameters.MODE)

    return error
```

## 4. Sub-Plugin Execution Pattern

Every sub-plugin phase follows the same template:

```python
def _runSubTask(self):
    """Run sub-task with full error handling."""
    error = CErrorReport()

    # 1. Create plugin
    try:
        self.subTaskPlugin = self.makePluginObject('task_name')
    except Exception as e:
        self.appendErrorReport(209,
            f'Failed to create task_name plugin: {e}\n{traceback.format_exc()}')
        error.append(self.__class__.__name__, 209,
                    f'Failed to create task_name: {e}', 'task_name', 4)
        return error

    try:
        plugin = self.subTaskPlugin

        # 2. Configure inputs and parameters
        plugin.container.inputData.XYZIN = self.container.inputData.XYZIN
        plugin.container.controlParameters.NCYCLES.set(5)

        # 3. Force synchronous if needed
        plugin.doAsync = False

        # 4. Execute
        print(f"[MyPipeline] Running task_name...")
        status = plugin.process()

        # 5. Check status
        if status != CPluginScript.SUCCEEDED:
            self.appendErrorReport(203, 'task_name failed')
            error.append(self.__class__.__name__, 203,
                        'task_name failed', 'task_name', 4)
            return error

        # 6. Verify critical outputs exist
        if not os.path.isfile(str(plugin.container.outputData.XYZOUT.fullPath)):
            self.appendErrorReport(211, 'task_name did not produce coordinate output')
            error.append(self.__class__.__name__, 211,
                        'task_name did not produce coordinate output', 'task_name', 4)
            return error

        # 7. Capture intermediate state for subsequent phases
        self.intermediateFile = plugin.container.outputData.XYZOUT

        # 8. Append XML
        self._appendPluginXml(plugin)

        print(f"[MyPipeline] task_name completed")

    except Exception as e:
        self.appendErrorReport(203,
            f'Exception in task_name: {e}\n{traceback.format_exc()}')
        error.append(self.__class__.__name__, 203,
                    f'Exception in task_name: {e}', 'task_name', 4)
        return error

    return error
```

### 4.1 Parameter Passing Techniques

```python
# Scalar values - use .set()
plugin.container.controlParameters.NCYCLES.set(5)
plugin.container.controlParameters.MODE.set('AUTO')

# File/object assignment - direct
plugin.container.inputData.XYZIN = self.container.inputData.XYZIN

# Bulk copy - .copyData()
plugin.container.inputData.copyData(self.container.inputData, ['FIELD1', 'FIELD2'])

# Copy all control parameters with exclusion list
EXCLUDED_PARAMS = {'VALIDATE_IRIS', 'REFMAC_CLEANUP', 'RUN_ADP_ANALYSIS'}
for attr in self.container.controlParameters.dataOrder():
    if attr not in EXCLUDED_PARAMS:
        setattr(result.container.controlParameters, attr,
                getattr(self.container.controlParameters, attr))
```

## 5. Error Handling

### 5.1 Dual Error Reporting

Always report errors through **both** channels:

```python
# Channel 1: Pipeline's error report XML (stored in diagnostic.xml)
self.appendErrorReport(error_code, 'Human-readable description')

# Channel 2: CErrorReport return value (used by process() for flow control)
error.append(self.__class__.__name__, error_code,
            'Description', 'method_name', severity)
```

### 5.2 Severity Levels

| Level | Constant | Meaning | Pipeline behavior |
|-------|----------|---------|-------------------|
| 0 | SEVERITY_OK | No error | Continue |
| 2 | SEVERITY_WARNING | Non-blocking | Continue |
| 3 | SEVERITY_UNDEFINED_ERROR | Uncertain | Continue |
| 4 | SEVERITY_ERROR | Blocking | Return early |

### 5.3 Error Checking Pattern

```python
result_error = self._runPhase()
if result_error and result_error.maxSeverity() >= 4:
    return result_error
```

### 5.4 Exception Wrapping

All plugin execution code must be wrapped in try/except with tracebacks:

```python
try:
    status = plugin.process()
except Exception as e:
    self.appendErrorReport(203,
        f'Exception in task_name: {e}\n{traceback.format_exc()}')
    error.append(...)
    return error
```

**Never use bare `except:`** - always catch `Exception` at minimum.

## 6. Output Harvesting

### 6.1 processOutputFiles() Pattern

```python
def processOutputFiles(self):
    error = CErrorReport()

    try:
        # Harvest from each phase that ran
        if self.phase1Plugin is not None:
            out = self.phase1Plugin.container.outputData
            if os.path.isfile(str(out.RESULT.fullPath)):
                self._harvestFile(out.RESULT, self.container.outputData.RESULT)

        # Flush final XML
        self._flushXML()

    except Exception as e:
        self.appendErrorReport(202,
            f'Exception in processOutputFiles: {e}\n{traceback.format_exc()}')
        error.append(self.__class__.__name__, 202,
                    f'Exception in processOutputFiles: {e}', 'harvest', 3)

    return error
```

### 6.2 File Harvesting

```python
def _harvestFile(self, sourceFile, destFile):
    """Copy file and metadata from source to destination."""
    try:
        shutil.copyfile(str(sourceFile.fullPath), str(destFile.fullPath))
        destFile.annotation = str(sourceFile.annotation)
        destFile.contentFlag = int(sourceFile.contentFlag)
        destFile.subType = int(sourceFile.subType)
    except Exception as e:
        self.appendErrorReport(202,
            f'Failed to harvest {sourceFile.fullPath} -> {destFile.fullPath}: {e}')
```

**Rules:**
- Always check `os.path.isfile()` before harvesting
- Copy metadata (annotation, contentFlag, subType), not just the file
- Only harvest from phases that actually ran (check plugin is not None)

## 7. XML Management

### 7.1 Accumulating Sub-Plugin XML

```python
def _appendPluginXml(self, plugin, tag=None):
    """Safely append sub-plugin XML to our xmlroot, replacing existing tag if present."""
    try:
        pluginRoot = CCP4Utils.openFileToEtree(plugin.makeFileName('PROGRAMXML'))
        servalcatXML = pluginRoot.xpath("//SERVALCAT")
        if len(servalcatXML) == 1:
            node = servalcatXML[0]
            if tag:
                node.tag = tag
                # Remove existing element with same tag to avoid duplicates
                # (important when real-time monitoring has already added interim data)
                existing = self.xmlroot.find(tag)
                if existing is not None:
                    self.xmlroot.remove(existing)
            self.xmlroot.append(node)
        else:
            self.xmlroot.append(pluginRoot)
        self._flushXML()
    except Exception as e:
        self.appendErrorReport(201, f'Failed to append plugin XML: {e}')
```

The optional `tag` parameter renames the sub-plugin's root element. This is essential
when a pipeline runs the same wrapper multiple times (e.g. servalcat before and after
water addition → `SERVALCAT_FIRST` and `SERVALCAT_WATERS`). The remove-before-append
pattern ensures that interim updates from real-time monitoring don't create duplicates.

### 7.2 Flushing XML

```python
def _flushXML(self):
    """Write current XML to program.xml via atomic tmp+move."""
    try:
        tmpFileName = self.makeFileName('PROGRAMXML') + '_tmp'
        with open(tmpFileName, 'w') as f:
            CCP4Utils.writeXML(f, etree.tostring(self.xmlroot, pretty_print=True))
        shutil.move(tmpFileName, self.makeFileName('PROGRAMXML'))
    except Exception as e:
        self.appendErrorReport(201, f'Failed to write program.xml: {e}')
```

Use atomic tmp+move to prevent readers from seeing partially-written XML.

## 8. Conditional Execution

### 8.1 Toggle-Controlled Phases

When a phase is controlled by a user toggle, the toggle check should be **inside**
the private method:

```python
def _runOptionalPhase(self):
    """Run optional phase if enabled by user toggle."""
    error = CErrorReport()

    if not self.container.controlParameters.ENABLE_PHASE:
        return error  # Empty error = success, phase skipped

    # ... execute phase ...
    return error
```

This keeps `startProcess()` clean and uniform:

```python
def startProcess(self):
    error = CErrorReport()

    phase1_error = self._runOptionalPhase()  # Internally checks toggle
    if phase1_error and phase1_error.maxSeverity() >= 4:
        return phase1_error

    phase2_error = self._runRequiredPhase()
    if phase2_error and phase2_error.maxSeverity() >= 4:
        return phase2_error

    return error
```

### 8.2 Mode-Dependent Phases

```python
def startProcess(self):
    error = CErrorReport()

    if self._mode == 'ADVANCED':
        phase_error = self._runAdvancedStep()
    else:
        phase_error = self._runSimpleStep()

    if phase_error and phase_error.maxSeverity() >= 4:
        return phase_error

    return error
```

## 9. Anti-Patterns to Avoid

### 9.1 Asynchronous Signal Chains

```python
# WRONG: Async signal chain
plugin.connectSignal(plugin, 'finished', self.onFinished)
plugin.doAsync = True
plugin.process()
return CPluginScript.SUCCEEDED  # Returns before job finishes!

@QtCore.Slot(dict)
def onFinished(self, statusDict):
    if statusDict['finishStatus'] == CPluginScript.FAILED:
        self.reportStatus(CPluginScript.FAILED)
    else:
        self._nextPhase()  # Triggers another async chain
```

### 9.2 File Watching for Pipeline Flow Control

```python
# WRONG: Using file watching to chain pipeline phases
self.watchFile(xmlPath, handler=self.handleXmlChanged,
               minDeltaSize=34, unwatchWhileHandling=True)

@QtCore.Slot(str)
def handleXmlChanged(self, xmlFilename):
    # Pipeline phases should NOT be triggered by file system events
    self._nextPhase()
```

**Note:** File watching IS appropriate for **real-time progress reporting** at the
wrapper level (see Section 10). The anti-pattern is using it for **control flow**
between pipeline phases.

### 9.3 Bare Except Clauses

```python
# WRONG
try:
    something()
except:
    print('Failed')

# CORRECT
try:
    something()
except Exception as e:
    self.appendErrorReport(201, f'Failed: {e}\n{traceback.format_exc()}')
```

### 9.4 Debug Prints in Production

```python
# WRONG
print("AAA1")
print("AAA15.1")

# CORRECT
print(f"[MyPipeline] Phase 1 completed: {output_path}")
```

### 9.5 Validation with Side Effects

```python
# WRONG: validity() modifying qualifiers
def validity(self):
    self.container.someWrapper.inputData.XYZIN.set_qualifier('allowUndefined', True)
    return super().validity()

# CORRECT: Set allowUndefined in the DEF XML
# <qualifiers><allowUndefined>True</allowUndefined></qualifiers>
```

### 9.6 Monolithic finishUp() Methods

```python
# WRONG: One method doing 6 things
def finishUp(self):
    # copy files, set annotations, cleanup, validate, analyze, report

# CORRECT: Separate into processOutputFiles() + dedicated methods
def processOutputFiles(self):   # Harvest files
def _runValidation(self):       # Optional validation phase
def _runAnalysis(self):         # Optional analysis phase
```

## 10. Real-Time Progress Monitoring

When a pipeline runs a long-running wrapper (e.g. iterative refinement), the UI
needs cycle-by-cycle progress updates. This is achieved through a two-level
monitoring chain using pure Python threading (no Qt).

### 10.1 Architecture

```
External binary → output file (e.g. refined_stats.json)
  → Wrapper watchFile thread → handleJsonChanged → wrapper program.xml
    → Pipeline monitor thread → _syncSubPluginXml → pipeline program.xml
      → UI polls pipeline program.xml for display
```

**Level 1 (Wrapper):** The base class `watchFile()` method starts a daemon thread
that polls a file for size changes. When the external binary writes new output, the
handler converts it and writes to the wrapper's `program.xml`.

**Level 2 (Pipeline):** The pipeline starts its own daemon thread that polls the
wrapper's `program.xml`. When it grows, the pipeline reads it, renames the root tag
(e.g. `SERVALCAT` → `SERVALCAT_FIRST`), and writes to the pipeline's `program.xml`.

### 10.2 Wrapper: startProcess() with File Watching

Override `startProcess()` to set up file watching **before** starting the process.
Call `super().startProcess()` to actually launch the subprocess:

```python
def startProcess(self):
    """Set up real-time progress monitoring, then start the external process."""
    jsonFilePath = os.path.join(self.getWorkDirectory(), "refined_stats.json")
    self.watchFile(jsonFilePath, handler=self.handleJsonChanged,
                   unwatchWhileHandling=True)
    return super().startProcess()
```

The `unwatchWhileHandling=True` flag prevents re-entrant handler calls while
the previous invocation is still processing.

**Important:** Always call `super().startProcess()` at the end. In the Django
baselayer, `startProcess()` is the method that actually launches the subprocess
(delegating to `_startProcessSync()` or `_startProcessAsync()`). Returning
`SUCCEEDED` without calling super means the external command never runs.

### 10.3 Wrapper: Progress Handler

```python
def handleJsonChanged(self, jsonFilePath):
    """Parse output and update program.xml in real-time."""
    self.xmlroot.clear()
    if os.path.isfile(jsonFilePath):
        try:
            with open(jsonFilePath, "r") as f:
                data = json.loads(f.read())
            xmlText = convert_to_xml(data)
            self.xmlroot = ET.fromstring(xmlText)
            self.flushXml()
        except Exception:
            traceback.print_exc()
```

The handler runs in a daemon thread, so it must be thread-safe. The `xmlroot.clear()`
+ full rebuild pattern is safe because each JSON snapshot is cumulative (contains all
cycles completed so far).

### 10.4 Pipeline: Monitoring Sub-Plugin XML

```python
import threading

def _monitorSubPluginXml(self, plugin, tag):
    """Start a background thread monitoring a sub-plugin's program.xml.

    Returns a stop function to call when monitoring should end.
    """
    stop_event = threading.Event()
    xml_path = plugin.makeFileName('PROGRAMXML')
    last_size = [0]

    def monitor():
        while not stop_event.wait(2.0):  # Poll every 2 seconds
            try:
                if os.path.isfile(xml_path):
                    current_size = os.path.getsize(xml_path)
                    if current_size > last_size[0]:
                        last_size[0] = current_size
                        self._syncSubPluginXml(xml_path, tag)
            except Exception:
                pass

    thread = threading.Thread(target=monitor, daemon=True,
                              name=f"XmlMonitor-{tag}")
    thread.start()

    def stop():
        stop_event.set()
        thread.join(timeout=5)

    return stop

def _syncSubPluginXml(self, xml_path, tag):
    """Read sub-plugin XML, rename root tag, and update pipeline XML."""
    try:
        pluginRoot = CCP4Utils.openFileToEtree(xml_path)
        servalcatXML = pluginRoot.xpath("//SERVALCAT")
        if len(servalcatXML) == 1:
            node = servalcatXML[0]
            node.tag = tag
            existing = self.xmlroot.find(tag)
            if existing is not None:
                self.xmlroot.remove(existing)
            self.xmlroot.append(node)
            self._flushXML()
    except Exception:
        pass
```

### 10.5 Pipeline: Using the Monitor

```python
def _runRefinement(self):
    error = CErrorReport()
    plugin = self.makePluginObject('servalcat')
    # ... configure plugin ...
    plugin.doAsync = False

    # Start monitoring BEFORE process() blocks
    stopMonitor = self._monitorSubPluginXml(plugin, 'SERVALCAT_FIRST')

    status = plugin.process()  # Blocks while monitor thread updates XML

    # Stop monitoring BEFORE final sync
    stopMonitor()

    if status == CPluginScript.FAILED:
        # ... error handling ...
        return error

    # Final sync gets complete XML including processOutputFiles() data
    self._appendPluginXml(plugin, tag='SERVALCAT_FIRST')
    return error
```

**Key points:**
- Start the monitor **before** `plugin.process()` blocks
- Stop the monitor **before** the final `_appendPluginXml()` call
- The final `_appendPluginXml()` replaces interim data with the complete result
- Use `daemon=True` threads so they don't block process exit
- `threading.Event.wait(2.0)` provides both the polling interval and clean shutdown

## 11. Complete Pipeline Template

*(Note: For brevity, real-time monitoring code is omitted — see Section 10 for
the `_monitorSubPluginXml` and `_syncSubPluginXml` methods to add when wrapping
long-running sub-plugins.)*

```python
import os
import shutil
import traceback

from lxml import etree

from ccp4i2.core import CCP4Utils
from ccp4i2.core.CCP4PluginScript import CPluginScript
from ccp4i2.core.CCP4ErrorHandling import CErrorReport


class MyPipeline(CPluginScript):

    TASKMODULE = 'category'
    TASKTITLE = 'My Pipeline Title'
    TASKNAME = 'my_pipeline'
    MAINTAINER = 'dev@example.com'
    WHATNEXT = ['next_task']

    ERROR_CODES = {
        201: {'description': 'General pipeline failure'},
        202: {'description': 'Failed in harvest operation'},
        203: {'description': 'Failed in phase 1'},
        204: {'description': 'Failed in phase 2'},
        209: {'description': 'Failed to create sub-plugin'},
        211: {'description': 'Missing required output'},
        212: {'description': 'Invalid input configuration'},
    }

    def __init__(self, *args, **kws):
        super(MyPipeline, self).__init__(*args, **kws)
        self.xmlroot = etree.Element('MyPipeline')

        # Intermediate state
        self.intermediateCoords = None

        # Sub-plugin references
        self.phase1Plugin = None
        self.phase2Plugin = None

        # Pipeline mode flags
        self._mode = None

    def processInputFiles(self):
        error = CErrorReport()
        invalidFiles = self.checkInputData()
        if len(invalidFiles) > 0:
            for f in invalidFiles:
                self.appendErrorReport(212, f'Missing required input: {f}')
            error.append(self.__class__.__name__, 212,
                        f'Missing inputs: {invalidFiles}', 'processInputFiles', 4)
            return error
        self._mode = str(self.container.controlParameters.MODE)
        return error

    def startProcess(self):
        error = CErrorReport()

        # Phase 1: Optional step
        phase1_error = self._runPhase1()
        if phase1_error and phase1_error.maxSeverity() >= 4:
            return phase1_error

        # Phase 2: Required step
        phase2_error = self._runPhase2()
        if phase2_error and phase2_error.maxSeverity() >= 4:
            return phase2_error

        return error

    def _runPhase1(self):
        error = CErrorReport()
        if not self.container.controlParameters.ENABLE_PHASE1:
            return error

        try:
            self.phase1Plugin = self.makePluginObject('task1')
        except Exception as e:
            self.appendErrorReport(209,
                f'Failed to create task1: {e}\n{traceback.format_exc()}')
            error.append(self.__class__.__name__, 209,
                        f'Failed to create task1: {e}', 'task1', 4)
            return error

        try:
            plugin = self.phase1Plugin
            plugin.container.inputData.XYZIN = self.container.inputData.XYZIN
            plugin.container.controlParameters.OPTION.set('VALUE')

            print(f"[MyPipeline] Running task1...")
            status = plugin.process()

            if status != CPluginScript.SUCCEEDED:
                self.appendErrorReport(203, 'task1 failed')
                error.append(self.__class__.__name__, 203,
                            'task1 failed', 'task1', 4)
                return error

            self.intermediateCoords = plugin.container.outputData.XYZOUT
            self._appendPluginXml(plugin)
            print(f"[MyPipeline] task1 completed")

        except Exception as e:
            self.appendErrorReport(203,
                f'Exception in task1: {e}\n{traceback.format_exc()}')
            error.append(self.__class__.__name__, 203,
                        f'Exception in task1: {e}', 'task1', 4)
            return error

        return error

    def _runPhase2(self):
        error = CErrorReport()
        # ... same pattern ...
        return error

    def processOutputFiles(self):
        error = CErrorReport()
        try:
            if self.phase1Plugin is not None:
                out = self.phase1Plugin.container.outputData
                if os.path.isfile(str(out.RESULT.fullPath)):
                    self._harvestFile(out.RESULT,
                                     self.container.outputData.RESULT)
            self._flushXML()
        except Exception as e:
            self.appendErrorReport(202,
                f'Exception in harvest: {e}\n{traceback.format_exc()}')
            error.append(self.__class__.__name__, 202,
                        f'Exception in harvest: {e}', 'harvest', 3)
        return error

    def _harvestFile(self, sourceFile, destFile):
        try:
            shutil.copyfile(str(sourceFile.fullPath),
                           str(destFile.fullPath))
            destFile.annotation = str(sourceFile.annotation)
            destFile.contentFlag = int(sourceFile.contentFlag)
            destFile.subType = int(sourceFile.subType)
        except Exception as e:
            self.appendErrorReport(202,
                f'Failed to harvest: {e}')

    def _appendPluginXml(self, plugin):
        try:
            pluginRoot = CCP4Utils.openFileToEtree(
                plugin.makeFileName('PROGRAMXML'))
            self.xmlroot.append(pluginRoot)
        except Exception as e:
            self.appendErrorReport(201,
                f'Failed to append plugin XML: {e}')

    def _flushXML(self):
        try:
            with open(self.makeFileName('PROGRAMXML'), 'w') as f:
                CCP4Utils.writeXML(f,
                    etree.tostring(self.xmlroot, pretty_print=True))
        except Exception as e:
            self.appendErrorReport(201,
                f'Failed to write program.xml: {e}')
```

## 12. Wrapper Best Practices

Wrappers (as opposed to pipelines) run a single external program. They are simpler:

- **Do** override `makeCommandAndScript()` and `processOutputFiles()`
- **Do** override `processInputFiles()` if input transformation is needed
- **Do** override `startProcess()` for real-time progress monitoring (Section 10.2),
  but always call `super().startProcess()` to actually launch the process
- **Do not** override `process()` directly
- **Do not** use bare `except:` clauses
- Use specific error codes for different failure modes
- Verify output files exist before processing them
