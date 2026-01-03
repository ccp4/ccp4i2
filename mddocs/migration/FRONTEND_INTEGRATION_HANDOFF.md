# Frontend Integration Handoff

This document summarizes the backend migration state for frontend integration work.

## Backend Migration Status: COMPLETE

The CCP4i2 backend has been migrated from Qt-based ccp4-python to a modern Python virtual environment with Django REST API. The API lifecycle is working end-to-end.

## Key Backend Changes Affecting Frontend

### 1. Execution Environment Change

**Before**: Jobs executed via `ccp4-python`
**After**: Jobs executed via project's `.venv/bin/python`

Location: [server/ccp4i2/lib/utils/jobs/context_run.py](../../server/ccp4i2/lib/utils/jobs/context_run.py)

```python
# Now uses project venv instead of ccp4-python
venv_python = project_root / ".venv" / "bin" / "python"
if not venv_python.exists():
    venv_python = project_root / ".venv.py311" / "bin" / "python"
```

**Frontend Impact**: Any Electron/desktop configuration that references `ccp4-python` must be updated.

### 2. API Response Format Standardization

All API endpoints now return consistent response structures:

```json
// Success response
{
  "success": true,
  "data": { ... },
  "status": 200
}

// Error response
{
  "success": false,
  "error": "Error message",
  "status": 400
}
```

**Frontend Impact**: Response handling may need updates if it expects different field names or structures.

### 3. Job Status Codes

| Status | Code | Description |
|--------|------|-------------|
| PENDING | 1 | Job created, not started |
| RUNNING | 3 | Job executing |
| FAILED | 5 | Job failed |
| FINISHED | 6 | Job completed successfully |

### 4. API Endpoints (Working)

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/api/projects/` | GET | List all projects |
| `/api/projects/<uuid>/` | GET | Get project details |
| `/api/jobs/create_task/` | POST | Create new job |
| `/api/jobs/<uuid>/upload_file_param/` | POST | Upload file and set parameter |
| `/api/jobs/<uuid>/set_parameter/` | POST | Set job parameter |
| `/api/jobs/<uuid>/run/` | POST | Start job execution |
| `/api/jobs/<uuid>/` | GET | Get job status/details |

### 5. File Upload Pattern

Files are uploaded via `upload_file_param` endpoint with:
- `file`: The file data (multipart)
- `param_path`: Dot-notation path to parameter (e.g., `inputData.XYZIN`)
- `subtype`: Optional subtype for file classification

Example:
```bash
curl -X POST "http://localhost:8000/api/jobs/${JOB_UUID}/upload_file_param/" \
  -F "file=@model.pdb" \
  -F "param_path=inputData.XYZIN"
```

### 6. Parameter Setting Pattern

Parameters set via `set_parameter` endpoint with:
- `path`: Dot-notation parameter path
- `value`: Parameter value (string, will be type-coerced)

Example:
```bash
curl -X POST "http://localhost:8000/api/jobs/${JOB_UUID}/set_parameter/" \
  -H "Content-Type: application/json" \
  -d '{"path": "controlParameters.NCYCLES", "value": "10"}'
```

## Project Structure (Reorganized)

```
ccp4i2/
├── server/ccp4i2/          # Django backend
│   ├── api/               # REST API endpoints
│   ├── db/                # Database models
│   └── lib/utils/jobs/    # Job execution (context_run.py)
├── mddocs/                # Documentation (reorganized)
│   ├── api/               # API docs
│   ├── architecture/      # System design docs
│   ├── cli/               # CLI reference
│   ├── migration/         # Migration docs (this file)
│   └── setup/             # Setup guides
├── tests/                 # Test files (reorganized)
│   ├── shell/             # Shell script tests
│   └── standalone/        # Python test files
├── scripts/               # Utility scripts
│   └── run_refmac_via_api.sh  # API lifecycle test script
├── .venv -> .venv.py311   # Active virtual environment
├── run_test.sh            # Test runner (sources CCP4 env)
└── CLAUDE.md              # Development guide
```

## Verified Working: Full API Lifecycle

The following workflow has been tested end-to-end:

1. Create job: `POST /api/jobs/create_task/`
2. Upload PDB: `POST /api/jobs/<uuid>/upload_file_param/`
3. Upload MTZ: `POST /api/jobs/<uuid>/upload_file_param/`
4. Set parameters: `POST /api/jobs/<uuid>/set_parameter/`
5. Run job: `POST /api/jobs/<uuid>/run/`
6. Poll status: `GET /api/jobs/<uuid>/`
7. Job completes with status 6 (FINISHED)

Reference script: [scripts/run_refmac_via_api.sh](../../scripts/run_refmac_via_api.sh)

## Frontend Integration Tasks

### Phase 1: Execution Environment
- [ ] Update Electron configuration to use `.venv/bin/python` instead of `ccp4-python`
- [ ] Ensure CCP4 environment variables are sourced before job execution
- [ ] Test local job execution from Electron app

### Phase 2: API Response Alignment
- [ ] Audit frontend API calls to identify expected response formats
- [ ] Map backend response fields to frontend expectations
- [ ] Update frontend response handlers as needed
- [ ] Test all API interactions

### Phase 3: End-to-End Validation
- [ ] Test project creation/listing
- [ ] Test job creation workflow
- [ ] Test file upload
- [ ] Test parameter setting
- [ ] Test job execution and status polling
- [ ] Test job output retrieval

## Environment Requirements

For the backend to work:
1. CCP4 environment sourced: `source /path/to/ccp4/bin/ccp4.setup-sh`
2. Virtual environment activated: `source .venv/bin/activate`
3. Django server running: `python server/manage.py runserver`

## Contact Points in Code

| Concern | File |
|---------|------|
| Job execution | `server/ccp4i2/lib/utils/jobs/context_run.py` |
| API endpoints | `server/ccp4i2/api/views/` |
| Job model | `server/ccp4i2/db/models/job.py` |
| Project model | `server/ccp4i2/db/models/project.py` |
| Parameter setting | `server/ccp4i2/lib/utils/jobs/params.py` |

---
*Generated: 2025-11-27*
*Backend Migration: Complete*
*Next: Frontend Integration*
