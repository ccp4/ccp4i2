# Logging Migration Plan - Replace Print Statements with Python Logging

**Date**: 2025-10-31
**Issue**: Excessive print statements (especially `[DEBUG]` output) interfere with JSON output in management commands and make testing difficult
**Solution**: Migrate to Python's standard `logging` module with configurable log levels

---

## Problem Summary

### Current Issues

1. **Debug output overwhelms test results**
   - `[DEBUG]` prints from CData classes flood stdout
   - Test scripts can't cleanly parse command output
   - JSON output from management commands gets mixed with debug messages

2. **Print statement locations** (audit results):
   - **Core CData system**: ~220 print statements across 17 files
   - **Server code**: ~261 print statements across 49 files
   - **Key offenders**:
     - `core/base_object/cdata.py` - 11 print statements (including `[DEBUG]` initialization)
     - `core/base_object/cdata_file.py` - 17 print statements (including `[DEBUG]` in setFullPath)
     - `server/ccp4x/db/async_db_handler.py` - 22 print statements
     - `server/ccp4x/lib/cdata_utils.py` - 16 print statements

3. **Example problematic output**:
   ```
   [DEBUG] Initializing CData instance of type CObsDataFile with name 'F_SIGF'
   [SETATTR] Setting attribute contentFlag on CObsDataFile
   [DEBUG] Branch: Value type smart assign for contentFlag
   Using selector: BestSelector
   File upload initiated for upload_id: xyz
   {"status": "success", "job_uuid": "..."}  # ← Actual output we want!
   ```

---

## Solution: Python Logging Module

### Benefits

1. **Configurable severity levels**: DEBUG, INFO, WARNING, ERROR, CRITICAL
2. **Easy to disable debug output**: Single configuration change
3. **Structured logging**: Consistent format across codebase
4. **Multiple handlers**: Console, file, syslog, etc.
5. **Performance**: Logging checks level before formatting expensive strings

### Logging Levels (Proposed Mapping)

| Current Print Statement | New Logging Level | When Shown |
|------------------------|-------------------|------------|
| `print(f"[DEBUG] ...")` | `logger.debug()` | LOG_LEVEL=DEBUG |
| `print(f"[SETATTR] ...")` | `logger.debug()` | LOG_LEVEL=DEBUG |
| `print(f"[WARNING] ...")` | `logger.warning()` | LOG_LEVEL=WARNING+ |
| `print("Using ...")` | `logger.info()` | LOG_LEVEL=INFO+ |
| `print("Error: ...")` | `logger.error()` | LOG_LEVEL=ERROR+ |
| Management command output | `self.stdout.write()` | Always (not logging) |

---

## Implementation Strategy

### Phase 1: Core Logging Setup (30 min)

Create centralized logging configuration:

**File**: `core/base_object/logging_config.py`

```python
"""
Centralized logging configuration for CCP4i2 codebase.
"""
import logging
import os
from typing import Optional

# Default log level from environment or INFO
DEFAULT_LOG_LEVEL = os.environ.get('CCP4_LOG_LEVEL', 'INFO').upper()

def get_logger(name: str, level: Optional[str] = None) -> logging.Logger:
    """
    Get a configured logger instance.

    Args:
        name: Logger name (typically __name__)
        level: Optional override level (DEBUG, INFO, WARNING, ERROR, CRITICAL)

    Returns:
        Configured logger instance
    """
    logger = logging.getLogger(name)

    if not logger.handlers:
        # Only configure if not already configured
        handler = logging.StreamHandler()
        formatter = logging.Formatter(
            fmt='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )
        handler.setFormatter(formatter)
        logger.addHandler(handler)

    # Set level
    log_level = level or DEFAULT_LOG_LEVEL
    logger.setLevel(getattr(logging, log_level))

    return logger

def configure_quiet_mode():
    """
    Suppress all output except ERROR and CRITICAL.
    Useful for management commands that need clean JSON output.
    """
    logging.basicConfig(level=logging.ERROR)

def configure_debug_mode():
    """Enable verbose debug output for troubleshooting."""
    logging.basicConfig(level=logging.DEBUG)
```

**File**: `server/ccp4x/config/logging_config.py` (Django-specific)

```python
"""
Django-specific logging configuration.
"""
import os

# Django LOGGING configuration
LOGGING = {
    'version': 1,
    'disable_existing_loggers': False,
    'formatters': {
        'standard': {
            'format': '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
        },
        'simple': {
            'format': '%(levelname)s - %(message)s'
        },
    },
    'handlers': {
        'console': {
            'class': 'logging.StreamHandler',
            'formatter': 'standard',
            'level': os.environ.get('CCP4_LOG_LEVEL', 'INFO'),
        },
        'file': {
            'class': 'logging.FileHandler',
            'filename': 'ccp4i2.log',
            'formatter': 'standard',
            'level': 'DEBUG',
        },
    },
    'loggers': {
        'ccp4x': {
            'handlers': ['console', 'file'],
            'level': 'DEBUG',
            'propagate': False,
        },
        'core': {
            'handlers': ['console', 'file'],
            'level': 'DEBUG',
            'propagate': False,
        },
    },
    'root': {
        'handlers': ['console'],
        'level': os.environ.get('CCP4_LOG_LEVEL', 'INFO'),
    },
}
```

---

### Phase 2: Migrate Core CData Classes (1 hour)

**Priority files** (most impactful):

1. **`core/base_object/cdata.py`** (11 prints)
2. **`core/base_object/cdata_file.py`** (17 prints)
3. **`core/CCP4PluginScript.py`** (75 prints)
4. **`core/base_object/hierarchy_system.py`** (15 prints)

**Migration pattern**:

```python
# BEFORE
def __init__(self, value=None, parent=None, name=None):
    super().__init__(parent=parent, name=name)
    print(f"[DEBUG] Initializing CData instance of type {self.__class__.__name__} with name '{name}'")

# AFTER
import logging
logger = logging.getLogger(__name__)

def __init__(self, value=None, parent=None, name=None):
    super().__init__(parent=parent, name=name)
    logger.debug("Initializing CData instance of type %s with name '%s'",
                 self.__class__.__name__, name)
```

**Key changes**:
- Add `logger = logging.getLogger(__name__)` at module top
- Replace `print(f"[DEBUG] ...")` with `logger.debug(...)`
- Replace `print(f"[WARNING] ...")` with `logger.warning(...)`
- Use `%s` placeholders instead of f-strings (lazy evaluation)

---

### Phase 3: Migrate Management Commands (1 hour)

**Files to update**:
- `server/ccp4x/db/management/commands/create_job.py`
- `server/ccp4x/db/management/commands/list_project.py`
- All new commands (clone_job, validate_job, etc.)

**Management command pattern**:

```python
# BEFORE
class Command(BaseCommand):
    def handle(self, *args, **options):
        print("Creating job...")
        # Do work
        print(f"Job created: {job_uuid}")
        return json.dumps({"status": "success"})

# AFTER
import logging
logger = logging.getLogger(__name__)

class Command(BaseCommand):
    def handle(self, *args, **options):
        # User-facing output (always shown)
        if not options.get('json'):
            self.stdout.write(self.style.SUCCESS('Creating job...'))

        # Internal logging (respects log level)
        logger.debug("Starting job creation with options: %s", options)

        # Do work

        # User-facing output
        if options.get('json'):
            self.stdout.write(json.dumps({"status": "success", "job_uuid": str(job_uuid)}))
        else:
            self.stdout.write(self.style.SUCCESS(f'Job created: {job_uuid}'))
```

**Key distinction**:
- `self.stdout.write()` - User-facing output (always shown)
- `logger.debug()` - Internal logging (respects LOG_LEVEL)

---

### Phase 4: Migrate Server Code (1 hour)

**Priority files**:
- `server/ccp4x/db/async_db_handler.py` (22 prints)
- `server/ccp4x/lib/cdata_utils.py` (16 prints)
- `server/ccp4x/api/JobViewSet.py` (7 prints)

**API endpoint pattern**:

```python
# BEFORE
@action(detail=True, methods=['post'])
def clone(self, request, pk=None):
    print(f"Cloning job {pk}")
    # Do work
    return Response(result.to_dict())

# AFTER
import logging
logger = logging.getLogger(__name__)

@action(detail=True, methods=['post'])
def clone(self, request, pk=None):
    logger.info("Cloning job %s for user %s", pk, request.user)
    # Do work
    logger.debug("Clone result: %s", result)
    return Response(result.to_dict())
```

---

## Usage Examples

### For Developers

**Default behavior** (INFO level):
```bash
python manage.py clone_job --jobuuid abc123
# Output:
# 2025-10-31 10:00:00 - ccp4x.db.commands.clone_job - INFO - Cloning job abc123
# ✓ Job cloned successfully
```

**Debug mode** (see everything):
```bash
export CCP4_LOG_LEVEL=DEBUG
python manage.py clone_job --jobuuid abc123
# Output:
# 2025-10-31 10:00:00 - core.base_object.cdata - DEBUG - Initializing CData instance...
# 2025-10-31 10:00:00 - core.base_object.cdata - DEBUG - Branch: CData to CData for contentFlag
# 2025-10-31 10:00:00 - ccp4x.db.commands.clone_job - INFO - Cloning job abc123
# 2025-10-31 10:00:00 - ccp4x.lib.utils.jobs.clone - DEBUG - Loading source job container
# ✓ Job cloned successfully
```

**Quiet mode** (clean JSON for scripts):
```bash
export CCP4_LOG_LEVEL=ERROR
python manage.py clone_job --jobuuid abc123 --json
# Output:
# {"status": "success", "new_job_uuid": "def456", "new_job_number": "2"}
```

---

### For Testing

**Test scripts with clean output**:

```bash
#!/bin/bash
# test_complete_workflow.sh with logging

# Suppress debug output for clean parsing
export CCP4_LOG_LEVEL=ERROR

# Now commands produce clean JSON output
JOB_OUTPUT=$(python manage.py create_job -pn "$PROJECT_NAME" -tn parrot --json)
JOB_UUID=$(echo "$JOB_OUTPUT" | jq -r '.job_uuid')

# No more grepping out [DEBUG] lines!
```

**pytest configuration**:

```python
# tests/conftest.py
import logging
import pytest

@pytest.fixture(autouse=True)
def configure_logging(request):
    """Configure logging for tests."""
    # Suppress debug output unless running with -vv
    if request.config.getoption("verbose") < 2:
        logging.basicConfig(level=logging.ERROR)
    else:
        logging.basicConfig(level=logging.DEBUG)
```

---

## Migration Checklist

### Phase 1: Setup (30 min)
- [ ] Create `core/base_object/logging_config.py`
- [ ] Create `server/ccp4x/config/logging_config.py`
- [ ] Add logging config to Django settings
- [ ] Test basic logger functionality

### Phase 2: Core CData (1 hour)
- [ ] Migrate `core/base_object/cdata.py`
- [ ] Migrate `core/base_object/cdata_file.py`
- [ ] Migrate `core/base_object/hierarchy_system.py`
- [ ] Test CData initialization with different log levels

### Phase 3: Management Commands (1 hour)
- [ ] Update base command template
- [ ] Migrate all new commands (clone_job, validate_job, etc.)
- [ ] Migrate existing commands (create_job, list_project)
- [ ] Test JSON output with LOG_LEVEL=ERROR

### Phase 4: Server Code (1 hour)
- [ ] Migrate `server/ccp4x/db/async_db_handler.py`
- [ ] Migrate `server/ccp4x/lib/cdata_utils.py`
- [ ] Migrate ViewSets in `server/ccp4x/api/`
- [ ] Test API responses

### Phase 5: Testing & Documentation (30 min)
- [ ] Update test scripts to use LOG_LEVEL
- [ ] Update `test_complete_workflow.sh`
- [ ] Add logging examples to README
- [ ] Document LOG_LEVEL environment variable

### Phase 6: Cleanup (30 min)
- [ ] Search for remaining print statements
- [ ] Add linting rule to prevent new prints
- [ ] Update contribution guidelines

---

## Estimated Timeline

| Phase | Duration | Description |
|-------|----------|-------------|
| Phase 1 | 30 min | Setup logging infrastructure |
| Phase 2 | 1 hour | Core CData classes |
| Phase 3 | 1 hour | Management commands |
| Phase 4 | 1 hour | Server code |
| Phase 5 | 30 min | Testing & docs |
| Phase 6 | 30 min | Cleanup |
| **Total** | **4.5 hours** | Complete migration |

---

## Expected Benefits

### Immediate Benefits
1. ✅ Clean JSON output from management commands
2. ✅ Test scripts can parse output reliably
3. ✅ `export CCP4_LOG_LEVEL=ERROR` gives silent mode

### Long-term Benefits
1. ✅ Professional logging infrastructure
2. ✅ Debugging flexibility (change level without code changes)
3. ✅ Log to files for production debugging
4. ✅ Integration with log aggregation systems (future)
5. ✅ Better performance (lazy evaluation of log messages)

---

## Risk Mitigation

### Low Risk Changes
- Adding logging alongside existing prints (transitional)
- Gradually removing prints after validation

### Testing Strategy
1. Keep existing prints initially, add logging
2. Test with LOG_LEVEL=DEBUG (should see both)
3. Test with LOG_LEVEL=ERROR (should see only errors)
4. Remove prints after validation

### Rollback Plan
- All changes in git, easy to revert if issues
- Can keep prints in critical paths temporarily
- Phased migration allows incremental validation

---

## Linting Rules (Future)

Add to `.pylintrc` or `pyproject.toml`:

```toml
[tool.ruff]
# Forbid print statements outside of management commands
select = ["T20"]  # flake8-print

[tool.ruff.per-file-ignores]
"**/management/commands/*.py" = ["T20"]  # Allow prints in commands (for stdout.write)
"tests/**/*.py" = ["T20"]  # Allow prints in tests
```

---

## Priority Recommendation

**Start with Phase 1 + Phase 3** (2 hours total):

1. **Phase 1**: Setup logging infrastructure (30 min)
2. **Phase 3**: Migrate management commands (1 hour)
3. **Test**: Run `test_complete_workflow.sh` with `LOG_LEVEL=ERROR`

**This solves the immediate pain point** (clean test output) while establishing the pattern for the rest of the codebase.

**Then continue with Phase 2** (Core CData) to eliminate the bulk of debug output.

---

**Status**: Ready to implement
**Blocking Issues**: None
**Dependencies**: Python 3.6+ (logging module is standard library)

---

**Next Step**: Shall we start with Phase 1 + Phase 3 to fix the immediate testing issue?
