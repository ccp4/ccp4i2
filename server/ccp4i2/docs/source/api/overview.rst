API Overview
============

Base URL: ``http://localhost:8000/api/``

Core Resources
--------------

Projects
~~~~~~~~

- ``GET /api/projects/`` - List projects
- ``POST /api/projects/`` - Create project
- ``GET /api/projects/{uuid}/`` - Get project

Jobs
~~~~

- ``GET /api/jobs/`` - List jobs
- ``POST /api/jobs/`` - Create job
- ``POST /api/jobs/{uuid}/run/`` - Execute job
- ``GET /api/jobs/{uuid}/validation/`` - Validate job

See ``mddocs/api/API_OVERVIEW.md`` for complete documentation.
