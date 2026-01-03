# CCP4i2 REST API Overview

This document provides an overview of the CCP4i2 REST API endpoints and their usage.

## Base URL

All API endpoints are accessible at:

```
http://localhost:8000/api/
```

## Authentication

Currently, the API uses session-based authentication. For programmatic access, ensure cookies are preserved across requests.

---

## Core Resources

### Projects

| Method | Endpoint | Description |
|--------|----------|-------------|
| GET | `/api/projects/` | List all projects |
| POST | `/api/projects/` | Create a new project |
| GET | `/api/projects/{uuid}/` | Get project details |
| DELETE | `/api/projects/{uuid}/` | Delete a project |
| GET | `/api/projects/{uuid}/jobs/` | List jobs in project |
| POST | `/api/projects/{uuid}/export/` | Export project to ZIP |

### Jobs

| Method | Endpoint | Description |
|--------|----------|-------------|
| GET | `/api/jobs/` | List all jobs |
| POST | `/api/jobs/` | Create a new job |
| GET | `/api/jobs/{uuid}/` | Get job details |
| DELETE | `/api/jobs/{uuid}/` | Delete a job |
| POST | `/api/jobs/{uuid}/run/` | Execute a job |
| GET | `/api/jobs/{uuid}/validation/` | Validate job parameters |
| POST | `/api/jobs/{uuid}/set_parameter/` | Set a job parameter |
| GET | `/api/jobs/{uuid}/params_xml/` | Get job parameters as XML |
| GET | `/api/jobs/{uuid}/report_xml/` | Get job report XML |
| GET | `/api/jobs/{uuid}/container/` | Get job container structure |

### Files

| Method | Endpoint | Description |
|--------|----------|-------------|
| GET | `/api/files/` | List all files |
| GET | `/api/files/{uuid}/` | Get file details |
| POST | `/api/files/{uuid}/preview/` | Preview file contents |

### Plugins

| Method | Endpoint | Description |
|--------|----------|-------------|
| GET | `/api/plugins/` | List available plugins |
| GET | `/api/plugins/{name}/` | Get plugin details |

---

## Job Lifecycle Endpoints

### Parameter Management

#### Set Parameter
```http
POST /api/jobs/{uuid}/set_parameter/
Content-Type: application/json

{
    "object_path": "container.inputData.NCYCLES",
    "value": 10
}
```

#### Upload File Parameter
```http
POST /api/jobs/{uuid}/upload_file_param/
Content-Type: multipart/form-data

file: <binary>
job_param_name: XYZIN
```

#### Set Context Job
```http
POST /api/jobs/{uuid}/set_context_job/
Content-Type: application/json

{
    "context_job_uuid": "abc123..."
}
```

### Validation

#### Validate Job
```http
GET /api/jobs/{uuid}/validation/

Response:
{
    "status": "Success",
    "xml": "<errorReportList>...</errorReportList>",
    "maxSeverity": 0,
    "hasErrors": false
}
```

Severity levels:
- 0: OK (no errors)
- 1: UNDEFINED (values not set)
- 2: WARNING (non-blocking issues)
- 3: UNDEFINED_ERROR (required values missing)
- 4: ERROR (fatal errors)

### Execution

#### Run Job
```http
POST /api/jobs/{uuid}/run/
Content-Type: application/json

{
    "detach": true  // Optional: run in background
}
```

#### Get Job Status
```http
GET /api/jobs/{uuid}/

Response:
{
    "uuid": "...",
    "status": 3,  // 1=pending, 2=queued, 3=running, 4=finished, 5=failed...
    "title": "Refinement round 1",
    ...
}
```

### Reports and Output

#### Get Report XML
```http
GET /api/jobs/{uuid}/report_xml/

Response:
{
    "success": true,
    "data": {
        "xml": "<Report>...</Report>"
    }
}
```

#### Get Diagnostic XML
```http
GET /api/jobs/{uuid}/diagnostic_xml/

Response:
{
    "success": true,
    "data": {
        "xml": "<diagnostic>...</diagnostic>"
    }
}
```

---

## Proxy Endpoints

For frontend compatibility, proxy endpoints are available:

```
/api/proxy/projects/{id}/
/api/proxy/jobs/{id}/
/api/proxy/jobs/{id}/report_xml/
/api/proxy/jobs/{id}/validation/
```

These wrap the main API endpoints with consistent response formatting.

---

## Error Responses

All errors follow a consistent format:

```json
{
    "status": "Failed",
    "reason": "Error description",
    "detail": "Additional details (optional)"
}
```

HTTP status codes:
- `400` - Bad Request (invalid parameters)
- `404` - Not Found (resource doesn't exist)
- `500` - Server Error (unexpected failure)

---

## Example Workflows

### Create and Run a Job

```python
import requests

BASE = "http://localhost:8000/api"

# 1. Create project
project = requests.post(f"{BASE}/projects/", json={
    "name": "my_project",
    "description": "Test project"
}).json()

# 2. Create job
job = requests.post(f"{BASE}/jobs/", json={
    "project": project["uuid"],
    "task_name": "refmac5",
    "title": "Refinement"
}).json()

# 3. Set parameters
requests.post(f"{BASE}/jobs/{job['uuid']}/set_parameter/", json={
    "object_path": "container.inputData.XYZIN.fullPath",
    "value": "/path/to/model.pdb"
})

requests.post(f"{BASE}/jobs/{job['uuid']}/set_parameter/", json={
    "object_path": "container.controlParameters.NCYCLES",
    "value": 10
})

# 4. Validate
validation = requests.get(f"{BASE}/jobs/{job['uuid']}/validation/").json()
if validation["hasErrors"]:
    print(f"Validation failed: {validation['xml']}")
    exit(1)

# 5. Run
requests.post(f"{BASE}/jobs/{job['uuid']}/run/", json={"detach": True})

# 6. Poll for completion
import time
while True:
    status = requests.get(f"{BASE}/jobs/{job['uuid']}/").json()
    if status["status"] in [4, 5, 6]:  # finished, failed, or stopped
        break
    time.sleep(5)

# 7. Get report
report = requests.get(f"{BASE}/jobs/{job['uuid']}/report_xml/").json()
```

---

## See Also

- [API Endpoints Analysis](API_ENDPOINTS_ANALYSIS.md) - Detailed endpoint analysis
- [CLI Reference](../cli/CLI.md) - Command-line interface
- [i2run Guide](../cli/I2RUN_GUIDE.md) - Task runner command
