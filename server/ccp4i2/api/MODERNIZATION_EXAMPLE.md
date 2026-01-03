# JobViewSet Modernization Example

## Step-by-Step: Updating the `container` Endpoint

This example shows how to update the JobViewSet `container` endpoint from legacy `job_utils` to modern `lib.utils` approach.

## Current Implementation (Legacy)

**File**: `server/ccp4i2/api/JobViewSet.py` (lines ~792-800)

```python
# Legacy imports
from ..lib.job_utils.get_job_container import get_job_container
from ..lib.job_utils.json_for_job_container import json_for_job_container

@action(detail=True, methods=["get"], permission_classes=[])
def container(self, request, pk=None):
    """Get job container as JSON."""
    try:
        job = models.Job.objects.get(id=pk)
        result = json_for_job_container(job)  # Returns dict directly
        return JsonResponse({"status": "Success", "result": result})
    except Exception as e:
        logger.exception("Failed to get container for job %s", pk, exc_info=e)
        return JsonResponse({"status": "Failed", "reason": str(e)}, status=500)
```

### Problems with Legacy Approach:
1. **No Result pattern**: Returns dict directly, inconsistent error handling
2. **No database context**: May not sync properly
3. **Implicit failures**: Exception handling is catch-all
4. **Inconsistent**: Management commands use different pattern

## Modern Implementation

### Step 1: Update Imports

Replace:
```python
from ..lib.job_utils.get_job_container import get_job_container
from ..lib.job_utils.json_for_job_container import json_for_job_container
```

With:
```python
from ..lib.utils.plugins.get_plugin import get_job_plugin
from ..lib.utils.containers.json_encoder import CCP4i2JsonEncoder
import json
```

### Step 2: Update Endpoint Implementation

```python
@action(detail=True, methods=["get"], permission_classes=[])
def container(self, request, pk=None):
    """
    Get job container as JSON.

    Returns the job's CContainer serialized to JSON format, including:
    - All parameters (input/output/control)
    - File references
    - Metadata
    - Current values

    Returns:
        Response: JSON with container data
    """
    try:
        job = models.Job.objects.get(id=pk)

        # Get plugin (loads container with params)
        plugin = get_job_plugin(job)

        # Serialize container to JSON
        container_json = json.dumps(
            plugin.container,
            cls=CCP4i2JsonEncoder,
            indent=2
        )

        return Response({
            "status": "Success",
            "result": json.loads(container_json)  # Convert back to dict for Response
        })

    except models.Job.DoesNotExist:
        logger.warning("Job %s not found", pk)
        return Response(
            {"status": "Failed", "reason": "Job not found"},
            status=status.HTTP_404_NOT_FOUND
        )
    except Exception as e:
        logger.exception("Failed to get container for job %s", pk, exc_info=e)
        return Response(
            {"status": "Failed", "reason": str(e)},
            status=status.HTTP_500_INTERNAL_SERVER_ERROR
        )
```

### Step 3: Test

```bash
export CCP4I2_ROOT=$CCP4I2_ROOT
pytest server/ccp4i2/tests/api/test_viewsets_comprehensive.py::JobViewSetTests::test_job_container -xvs
```

## Alternative: Full Result Pattern

For even better consistency with management commands, we could create a wrapper utility:

### Create: `lib/utils/jobs/get_container_json.py`

```python
"""Get job container as JSON using Result pattern."""

import logging
import json
from typing import Dict, Any

from ccp4i2.db import models
from ccp4i2.lib.response import Result
from ccp4i2.lib.utils.plugins.get_plugin import get_job_plugin
from ccp4i2.lib.utils.containers.json_encoder import CCP4i2JsonEncoder

logger = logging.getLogger(__name__)


def get_job_container_json(job: models.Job) -> Result[Dict[str, Any]]:
    """
    Get job container as JSON dictionary.

    Args:
        job: Job model instance

    Returns:
        Result[Dict] containing container JSON or error
    """
    try:
        plugin = get_job_plugin(job)

        container_json = json.dumps(
            plugin.container,
            cls=CCP4i2JsonEncoder
        )

        return Result.ok(json.loads(container_json))

    except FileNotFoundError as e:
        return Result.fail(f"Job parameters not found: {e}")
    except Exception as e:
        logger.exception("Failed to get container JSON for job %s", job.uuid)
        return Result.fail(f"Failed to load container: {e}")
```

### Then in JobViewSet:

```python
from ..lib.utils.jobs.get_container_json import get_job_container_json

@action(detail=True, methods=["get"], permission_classes=[])
def container(self, request, pk=None):
    """Get job container as JSON."""
    try:
        job = models.Job.objects.get(id=pk)
        result = get_job_container_json(job)

        if result.success:
            return Response({
                "status": "Success",
                "result": result.data
            })
        else:
            return Response(
                {"status": "Failed", "reason": result.error},
                status=status.HTTP_400_BAD_REQUEST
            )

    except models.Job.DoesNotExist:
        return Response(
            {"status": "Failed", "reason": "Job not found"},
            status=status.HTTP_404_NOT_FOUND
        )
```

## Comparison: Before vs After

| Aspect | Legacy | Modern |
|--------|--------|--------|
| Error handling | Catch-all exception | Specific error types + Result pattern |
| Return type | Dict (untyped) | Result[Dict] (typed) |
| Database sync | Implicit | Explicit via plugin context |
| Consistency | Different from mgmt commands | Same as mgmt commands |
| Testability | Hard to mock | Easy to test with Result |
| Type safety | None | Full type hints |

## Benefits

1. **Consistency**: API endpoints and management commands use identical logic
2. **Type Safety**: Result[T] provides type information
3. **Better Errors**: Specific error messages, proper HTTP status codes
4. **Maintainability**: Single source of truth for container JSON
5. **Testability**: Can test `get_job_container_json()` independently

## Other Endpoints to Update

Following the same pattern:

### 1. `digest` and `digest_param_file`

Current:
```python
from ..lib.job_utils.digest_param_file import digest_param_file
```

Modern:
```python
from ..lib.utils.files.digest import digest_file
```

### 2. `upload_file_param`

Current:
```python
from ..lib.job_utils.upload_file_param import upload_file_param
```

Modern:
```python
from ..lib.utils.files.upload_param import upload_file_parameter
```

### 3. `object_method`

Current:
```python
from ..lib.job_utils.object_method import object_method
```

Modern:
```python
from ..lib.utils.helpers.object_method import execute_object_method
```

### 4. `what_next`

Current:
```python
from ..lib.job_utils.get_what_next import get_what_next
```

Modern:
```python
from ..lib.utils.navigation.what_next import get_what_next_suggestions
```

## Migration Checklist

For each endpoint:

- [ ] Identify modern utility in `lib/utils/`
- [ ] Check if it returns `Result[T]`
- [ ] Update imports in ViewSet
- [ ] Update endpoint to handle `Result` pattern
- [ ] Add proper error handling (404, 400, 500)
- [ ] Run specific test
- [ ] Verify behavior unchanged
- [ ] Update API documentation

## Notes

- Some utilities in `lib/utils/` are thin wrappers around `job_utils` (temporary)
- Focus on utilities that use `Result[T]` pattern
- Prioritize endpoints that break frequently
- Test after each change
- Can update incrementally (one endpoint at a time)
