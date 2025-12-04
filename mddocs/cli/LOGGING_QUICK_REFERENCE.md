# Logging Quick Reference

## Environment Variable

Control log output with `CCP4_LOG_LEVEL`:

```bash
export CCP4_LOG_LEVEL=ERROR    # Quiet mode - only errors (for testing/JSON output)
export CCP4_LOG_LEVEL=WARNING  # Only warnings and errors
export CCP4_LOG_LEVEL=INFO     # Normal output (default)
export CCP4_LOG_LEVEL=DEBUG    # Verbose - see everything
```

## Common Use Cases

### Testing / JSON Output
```bash
export CCP4_LOG_LEVEL=ERROR
python manage.py clone_job --jobuuid abc-123 --json
# Clean JSON output, no debug noise
```

### Normal Development
```bash
# No export needed - defaults to INFO
python manage.py create_job -pn my_project -tn parrot
# Shows normal progress messages
```

### Debugging Issues
```bash
export CCP4_LOG_LEVEL=DEBUG
python manage.py set_job_parameter --jobuuid abc-123 --path "inputData.XYZIN" --value "/path/to/file.pdb"
# See detailed CData attribute assignments, file path handling, etc.
```

### Unit Tests
```bash
export CCP4_LOG_LEVEL=ERROR
pytest tests/
# Suppress debug output during testing
```

## Using Logging in Your Code

### Basic Pattern
```python
import logging
logger = logging.getLogger(__name__)

def my_function():
    logger.debug("Detailed debug info")      # Only with DEBUG
    logger.info("Normal operation")          # With INFO and above
    logger.warning("Something suspicious")   # With WARNING and above
    logger.error("Something failed")         # With ERROR and above
```

### Management Commands
```python
from django.core.management.base import BaseCommand
import logging

logger = logging.getLogger(__name__)

class Command(BaseCommand):
    def handle(self, *args, **options):
        # User output (always shown)
        self.stdout.write("Processing...")

        # Debug logging (respects CCP4_LOG_LEVEL)
        logger.debug("Internal state: %s", data)
```

## What Got Fixed

✅ **Before**: Overwhelming `[DEBUG]` messages polluted output
✅ **After**: Clean output with `CCP4_LOG_LEVEL=ERROR`

✅ **Before**: Couldn't parse JSON from commands
✅ **After**: Perfect JSON with `--json` flags

✅ **Before**: `print()` statements everywhere
✅ **After**: Professional logging with timestamps

## Files to Know

- **`core/base_object/logging_config.py`** - Core logging setup
- **`server/ccp4x/config/logging_config.py`** - Django logging config
- **`LOGGING_MIGRATION_PLAN.md`** - Full migration plan
- **`LOGGING_MIGRATION_COMPLETE.md`** - What was done

## Pro Tips

1. **Set in your shell profile** for persistent quiet mode:
   ```bash
   echo 'export CCP4_LOG_LEVEL=ERROR' >> ~/.zshrc
   ```

2. **Per-command override**:
   ```bash
   CCP4_LOG_LEVEL=DEBUG python manage.py create_job ...
   ```

3. **In test scripts**:
   ```bash
   #!/bin/bash
   export CCP4_LOG_LEVEL=ERROR
   # Rest of script has clean output
   ```
