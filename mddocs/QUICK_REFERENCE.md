# Quick Reference: Async Execution Infrastructure

**For developers working with the new Qt-free CCP4i2 execution system**

---

## Essential Concepts

### 1. Plugin Types

```python
# Synchronous plugin (blocks until complete)
class MySyncPlugin(CPluginScript):
    ASYNCHRONOUS = False  # Default
    TASKCOMMAND = '/path/to/executable'

# Asynchronous plugin (returns immediately, monitors in background)
class MyAsyncPlugin(CPluginScript):
    ASYNCHRONOUS = True
    TASKCOMMAND = '/path/to/executable'
```

### 2. Template System

Two template types define how your plugin runs:

```python
class MyPlugin(CPluginScript):
    # COMLINETEMPLATE → command-line arguments
    COMLINETEMPLATE = '1 HKLIN $HKLIN\n1 HKLOUT $HKLOUT'

    # COMTEMPLATE → stdin (piped to program)
    COMTEMPLATE = '1 LABIN FP=FP SIGFP=SIGFP\n1 END'
```

**Variables**: Use `$VARNAME` or `$VARNAME.attribute`
- `$HKLIN` → looks up `container.find("HKLIN")` → calls `str(object)`
- `$CELL.a` → looks up `container.find("CELL").a` → calls `str(value)`

**Numeric Prefixes**: Lines like `"1 HKLIN"` have the `"1 "` stripped automatically.

### 3. Signal Flow in Pipelines

```
Plugin.process()
    ↓
startProcess() runs executable
    ↓
[Async: monitor in background] or [Sync: wait for completion]
    ↓
processOutputFiles() parses results
    ↓
reportStatus(SUCCEEDED/FAILED)
    ↓
finished signal emitted
    ↓
Parent's callback runs
```

---

## Common Patterns

### Pattern 1: Simple Plugin

```python
class ConvertData(CPluginScript):
    TASKNAME = 'convert_data'
    TASKCOMMAND = 'converter'
    ASYNCHRONOUS = False

    COMLINETEMPLATE = '1 INPUT $INPUT\n1 OUTPUT $OUTPUT'
    COMTEMPLATE = '1 MODE CONVERT\n1 END'

    def processOutputFiles(self):
        """Extract data from log file."""
        log = self.logFileText()
        # Parse log and populate self.container.outputData
        return CErrorReport()

# Usage
plugin = ConvertData(workDirectory='/tmp/work', name='convert')
plugin.container.inputData.INPUT.setFullPath('/data/input.dat')
plugin.container.outputData.OUTPUT.setFullPath('/data/output.dat')
status = plugin.process()  # Blocks until complete
```

### Pattern 2: Two-Step Pipeline

```python
class TwoStepPipeline(CPluginScript):
    TASKNAME = 'two_step'
    ASYNCHRONOUS = True

    def process(self):
        # Validate inputs
        if self.checkInputData():
            return self.reportStatus(self.FAILED)

        # Set default output names
        self.checkOutputData()

        # Create and run first plugin
        self.step1 = self.makePluginObject('plugin_one')
        self.step1.container.inputData.HKLIN.set(self.container.inputData.HKLIN)
        self.connectSignal(self.step1, 'finished', self.on_step1_done)
        self.step1.process()

    def on_step1_done(self, status):
        """Callback when step 1 finishes."""
        if status != self.SUCCEEDED:
            return self.reportStatus(self.FAILED)

        # Create and run second plugin
        self.step2 = self.makePluginObject('plugin_two')
        self.step2.container.inputData.HKLIN.set(self.step1.container.outputData.HKLOUT)
        self.step2.container.outputData.HKLOUT.set(self.container.outputData.HKLOUT)

        # Use postProcessWrapper to auto-forward final status
        self.connectSignal(self.step2, 'finished', self.postProcessWrapper)
        self.step2.process()

# Usage
pipeline = TwoStepPipeline(workDirectory='/tmp/pipeline', name='pipe')
pipeline.container.inputData.HKLIN.setFullPath('/data/input.mtz')
pipeline.finished.connect(lambda d: print(f"Done! Status={d['finishStatus']}"))
pipeline.process()
```

### Pattern 3: Conditional Branching

```python
class ConditionalPipeline(CPluginScript):
    def process(self):
        # Check input type
        input_file = self.container.inputData.HKLIN
        if self._is_unmerged(input_file):
            self.branch = self.makePluginObject('merge_data')
        else:
            self.branch = self.makePluginObject('validate_data')

        self.branch.container.inputData.HKLIN.set(input_file)
        self.connectSignal(self.branch, 'finished', self.on_branch_done)
        self.branch.process()

    def _is_unmerged(self, mtz_file):
        """Check if MTZ file is unmerged."""
        # Read MTZ header, check for BATCH column, etc.
        pass
```

### Pattern 4: Parallel Execution

```python
class ParallelPipeline(CPluginScript):
    def process(self):
        self.completed_count = 0

        # Launch multiple plugins in parallel
        self.plugin_a = self.makePluginObject('task_a')
        self.plugin_b = self.makePluginObject('task_b')
        self.plugin_c = self.makePluginObject('task_c')

        for plugin in [self.plugin_a, self.plugin_b, self.plugin_c]:
            self.connectSignal(plugin, 'finished', self.on_one_complete)
            plugin.process()

    def on_one_complete(self, status):
        """Called each time one plugin finishes."""
        if status != self.SUCCEEDED:
            return self.reportStatus(self.FAILED)

        self.completed_count += 1
        if self.completed_count == 3:
            # All done!
            self.reportStatus(self.SUCCEEDED)
```

---

## Important Methods

### CPluginScript Methods

```python
# Override these:
def process(self):
    """Main entry point. Start your pipeline here."""

def processOutputFiles(self):
    """Parse output files and populate container.outputData."""
    return CErrorReport()

# Call these:
self.makePluginObject('task_name')  # Create sub-plugin
self.connectSignal(obj, 'signal_name', callback)  # Connect signals
self.reportStatus(status)  # Emit finished signal
self.checkInputData()  # Validate input files exist
self.checkOutputData()  # Generate output file names
self.logFileText()  # Read log.txt contents
self.makeFileName('LOG')  # Get path to log.txt
```

### CData Methods

```python
# Get/Set values
obj.value = 42  # For fundamental types (CInt, CFloat, etc.)
obj.set({'field1': val1, 'field2': val2})  # Set multiple fields
obj.set(other_cdata_object)  # Copy from another CData
data_dict = obj.get()  # Extract as dict

# Value state
obj.isSet('field_name')  # Check if field has been set
obj.unSet('field_name')  # Clear field value
obj.setToDefault('field_name')  # Reset to default

# File operations
file_obj.setFullPath('/path/to/file.dat')
path = file_obj.getFullPath()
```

### HierarchicalObject Methods

```python
# Navigation
obj = container.find('HKLIN')  # Simple search
obj = container.find('inputData.protein.XYZIN')  # Path search

# Signals
obj.finished.connect(callback, weak=False)
obj.finished.emit({'finishStatus': 0})
self.connectSignal(sub_plugin, 'finished', self.on_done)

# Hierarchy
obj.set_parent(parent_obj)
parent = obj.get_parent()
children = obj.get_children()
path = obj.object_path()  # e.g., "pipeline.step1.HKLIN"
```

---

## Debugging Tips

### Enable Debug Logging

```python
# In your plugin or test:
import logging
logging.basicConfig(level=logging.DEBUG)

# You'll see:
# [DEBUG] process() calling makeCommandAndScript()
# [DEBUG] Template expanded to: 'HKLIN /path/to/file.mtz'
# [DEBUG] startProcess() running: /bin/prog HKLIN /path/to/file.mtz
```

### Check Job Directories

```bash
# After running pipeline, check subdirectories:
ls -R /tmp/my_pipeline/

/tmp/my_pipeline/:
pipeline.params.xml  log.txt  job_1/  job_2/

/tmp/my_pipeline/job_1/:
mtzdump.params.xml  log.txt  stderr.txt  com.txt

/tmp/my_pipeline/job_2/:
pdbset.params.xml  log.txt  stderr.txt  com.txt
```

### Inspect Signal Connections

```python
# Check what's connected to a signal:
print(f"Connections: {len(plugin.finished._connections)}")

# Check if callback is still alive (weak reference):
for ref in plugin.finished._connections:
    callback = ref()
    print(f"Callback: {callback}")
```

### Common Errors

**"Pipeline times out"**
- Check that `reportStatus()` is called after sync execution
- Verify signal connections aren't using lambdas without `weak=False`
- Ensure `finished` signal is emitted

**"Variable not found in container"**
- Check object names match template variables
- Use `container.find('VAR')` to test lookup
- Check object hierarchy with `obj.object_path()`

**"TypeError: set() expects dict"**
- Ensure CData objects have `get()` method
- Check that you're not passing primitive types
- Use `set({'field': value})` for single fields

**"FileNotFoundError: No such file or directory"**
- Ensure workDirectory is created before saveParams()
- Check that makePluginObject() creates job_N subdirectories
- Verify file paths are absolute, not relative

---

## Testing Your Plugin

### Unit Test Template

```python
import unittest
from pathlib import Path
import tempfile
from your_plugin import YourPlugin

class TestYourPlugin(unittest.TestCase):
    def test_basic_execution(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            # Setup
            plugin = YourPlugin(workDirectory=tmpdir, name='test')
            input_path = Path(tmpdir) / 'input.dat'
            output_path = Path(tmpdir) / 'output.dat'

            # Create test input
            with open(input_path, 'w') as f:
                f.write("test data\n")

            # Configure plugin
            plugin.container.inputData.INPUT.setFullPath(str(input_path))
            plugin.container.outputData.OUTPUT.setFullPath(str(output_path))

            # Run
            status = plugin.process()

            # Verify
            self.assertEqual(status, plugin.SUCCEEDED)
            self.assertTrue(output_path.exists())

            # Check output data was populated
            self.assertTrue(plugin.container.outputData.RESULT.isSet())
```

### Integration Test Template

```python
def test_pipeline_with_real_data(self):
    with tempfile.TemporaryDirectory() as tmpdir:
        # Track completion
        result = {'completed': False, 'status': None}

        def on_finished(status_dict):
            result['completed'] = True
            result['status'] = status_dict['finishStatus']

        # Setup pipeline
        pipeline = YourPipeline(workDirectory=tmpdir, name='test_pipe')
        pipeline.container.inputData.HKLIN.setFullPath('/data/test.mtz')
        pipeline.finished.connect(on_finished, weak=False)

        # Run
        pipeline.process()

        # Wait for completion
        for i in range(600):  # 60 seconds max
            if result['completed']:
                break
            time.sleep(0.1)

        # Verify
        self.assertTrue(result['completed'], "Pipeline did not complete")
        self.assertEqual(result['status'], pipeline.SUCCEEDED)
```

---

## Performance Guidelines

### Do's
- ✅ Use `ASYNCHRONOUS = False` for fast plugins (<1 second)
- ✅ Use `weak=True` (default) for signal connections to prevent memory leaks
- ✅ Create sub-plugins with `makePluginObject()` for proper isolation
- ✅ Call `reportStatus()` exactly once per plugin execution
- ✅ Use `find()` with paths for nested objects: `find("data.HKLIN")`

### Don'ts
- ❌ Don't call `process()` multiple times on same plugin instance
- ❌ Don't use lambdas with `weak=True` (they get garbage collected immediately)
- ❌ Don't forget to call `reportStatus()` in synchronous plugins
- ❌ Don't modify `container` during `processOutputFiles()` execution
- ❌ Don't create circular signal connections (A→B→A)

### Memory Management
```python
# ✅ Good: explicit reference kept
def on_done(status):
    print(f"Done: {status}")
plugin.finished.connect(on_done, weak=False)

# ❌ Bad: lambda collected immediately with weak=True
plugin.finished.connect(lambda s: print(s), weak=True)

# ✅ Good: use weak=False for lambdas
plugin.finished.connect(lambda s: print(s), weak=False)

# ✅ Best: use methods instead of lambdas
plugin.finished.connect(self.on_plugin_done, weak=True)
```

---

## Reference: Signal Types

All signals emit `dict` with these keys:

### `finished` Signal
```python
{
    'finishStatus': int,  # CPluginScript.SUCCEEDED or FAILED
    'jobId': str or None  # Database job ID if applicable
}
```

### Accessing in Callbacks

```python
# Modern style (dict parameter)
def on_finished(self, status_dict):
    status = status_dict['finishStatus']
    job_id = status_dict.get('jobId')

# Legacy style (int parameter) - automatically adapted
@QtCore.Slot(int)
def on_finished(self, status):
    # connectSignal() extracts int from dict automatically
    if status == CPluginScript.SUCCEEDED:
        # ...
```

---

## Quick Troubleshooting

| Symptom | Likely Cause | Solution |
|---------|-------------|----------|
| Pipeline never finishes | Missing `reportStatus()` | Add `reportStatus()` after sync execution |
| "Variable not found" | Object not in hierarchy | Check with `container.find('VAR')` |
| Signal not received | Weak reference collected | Use `weak=False` for lambdas |
| Test fails with FileNotFound | Directory not created | Ensure `workDirectory` exists |
| Process exits with code 127 | Command not found | Check `TASKCOMMAND` path, ensure CCP4 env loaded |
| Process succeeds but no output | Missing `processOutputFiles()` | Implement method to parse log and populate output |
| Segfault or crash | Circular references | Check parent/child relationships, avoid cycles |

---

## Further Reading

### Pipeline Development
- **[pipeline/ERROR_HANDLING_PATTERNS.md](pipeline/ERROR_HANDLING_PATTERNS.md)** - CErrorReport, try/except patterns, ERROR_CODES
- **[pipeline/VALIDITY_PATTERNS.md](pipeline/VALIDITY_PATTERNS.md)** - Content-aware validation with validity() overrides

### Architecture
- **MILESTONE_ASYNC_EXECUTION.md** - Comprehensive architecture documentation
- **CLAUDE.md** - Project overview and metadata system
- **core/CCP4PluginScript.py** - Main plugin implementation (heavily documented)
- **core/CCP4ComTemplate.py** - Template expansion system
- **tests/test_demo_copycell_integration.py** - Real-world pipeline example

---

**Document Version**: 1.0
**Last Updated**: October 30, 2025
**Maintained by**: Development Team
