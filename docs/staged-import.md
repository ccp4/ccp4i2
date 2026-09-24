# Staged import for web deployments

Large files (unmerged MTZ 50–200 MB, cryo-EM maps ~250 MB, project zips more)
can't be uploaded directly on a served/cloud deployment: they cross the Azure
Container Apps ingress, the Next middleware body cap (`middlewareClientMaxBodySize`,
100 MB), and Django's `DATA_UPLOAD_MAX_MEMORY_SIZE` (100 MB). Raising the caps is
the wrong fix on a shared deployment (that's what #384 decided). Instead the
browser delivers the file into a server-side **staging directory in chunks** that
each stay under every cap, then imports it by an **owner-bound handle**.

This finishes the a63 import-by-path work (#512/#516/#518): the backend already
had `CCP4I2_IMPORT_STAGING_DIR`, but nothing populated the directory and no client
could name a file in it. It also *replaces* the a63 cloud `local_path` path, which
was containment-only — see Security.

## The flow

```
begin ──▶ server makes  <STAGING_DIR>/<uuid>/  and a StagedUpload row
 ▲                                                    (owner-bound, state=staging)
 │  loop each chunk_bytes slice
 └── PUT chunks/<i> ──▶ <uuid>/part.<i>
finish ──▶ concat parts ─▶ <uuid>/<sanitised name>, verify size (+ sha256),
                            delete parts, state=ready
import  ──▶ upload_file_param / import_project with staged_upload=<uuid>
          ─▶ owner + containment check ─▶ copied into the project
          ─▶ state=consumed, directory deleted (or swept, for detached import)
```

## API (`/api/ccp4i2/staged-uploads/`, `IsAuthenticated`)

| Method | Path | Purpose | Refusals |
|--------|------|---------|----------|
| POST | `staged-uploads/` | **begin** → `{upload_id, chunk_bytes, expires_at}` | 413 over `MAX_BYTES`, 429 over per-owner in-flight |
| PUT | `staged-uploads/<uuid>/chunks/<i>` | a raw-bytes chunk → 204 | 404 unknown/foreign, 410 expired, 413 over `chunk_bytes` |
| POST | `staged-uploads/<uuid>/finish/` | assemble + verify → `{state: ready}` | 409 missing chunks, 422 size/hash mismatch |
| GET | `staged-uploads/<uuid>/` | `{state, present_chunks, ...}` for resume | 404 unknown/foreign |

Then `staged_upload=<uuid>` is accepted as a form field by
`POST jobs/<id>/upload_file_param/` and `POST projects/import_project/`.

The logic is in `server/ccp4i2/lib/utils/files/staged_upload.py`; the viewset
(`api/StagedUploadViewSet.py`) is a thin HTTP layer that maps typed errors to
those status codes.

## Security invariants

- **The server names every path.** Directory = the row's `uuid`; filename is
  sanitised with `get_valid_filename`. No client string is ever joined into a path.
- **Owner-bound.** A row belongs to `str(request.user.pk)` and can only be read,
  finished, or imported by that user. A foreign id is a 404, like an unknown one.
- **`local_path` is dead in a served deployment.** `resolve_importable_path` now
  honours a client-named path *only* for the desktop token; cloud names files by
  handle. This closes the a63 "any file inside the staging dir is importable by
  anyone who can name it" gap — the handle is unguessable and owner-bound.
- **Disk is bounded** by per-upload size (`MAX_BYTES`), per-owner in-flight count
  (`MAX_INFLIGHT`), and expiry (`TTL_HOURS`). `download_file` bounds nothing, so
  this is the only bound.
- **Failure is loud.** Every refusal is a distinct status code; nothing falls
  through to an empty `request.FILES`.
- **Copied, never adopted.** Staged bytes are copied into the project (via the
  existing `_LocalPathUpload`), so the staging directory can be swept freely.

## Consumption and the sweeper

`upload_file_param` consumes on success (mark + delete) — the bytes are copied
synchronously. `import_project` dispatches the importer **detached**, so it marks
the row consumed but keeps the file (`consume(delete=False)`); the sweeper reaps
the directory after the TTL. `manage.py sweep_staged_uploads` deletes rows past
the TTL and consumed rows whose file is gone; schedule it from the deployment's
maintenance job.

## Configuration

| Variable | Default | Meaning |
|----------|---------|---------|
| `CCP4I2_IMPORT_STAGING_DIR` | unset | **Setting it enables everything.** A directory the server can write, outside the projects tree. |
| `CCP4I2_IMPORT_STAGING_CHUNK_BYTES` | 16 MiB | Chunk size. Must stay under the smallest body cap. |
| `CCP4I2_IMPORT_STAGING_MAX_BYTES` | 2 GiB | Largest single staged file. |
| `CCP4I2_IMPORT_STAGING_TTL_HOURS` | 24 | Sweeper horizon. |
| `CCP4I2_IMPORT_STAGING_THRESHOLD_BYTES` | 32 MiB | Below this the client sends bytes as before. |
| `CCP4I2_IMPORT_STAGING_MAX_INFLIGHT` | 8 | Per-owner cap on staging-state uploads. |

When set, `version_info` advertises `import_staging: {chunk_bytes, max_bytes,
threshold_bytes}` and the client stages large files automatically. Absent, the
client uploads bytes exactly as before, and desktop is untouched (it imports by
`local_path`).

## Deploying (e.g. Materia on Azure)

1. Mount a writable directory outside the projects tree (Materia: a folder on the
   Azure Files share the server container already mounts) and set
   `CCP4I2_IMPORT_STAGING_DIR` to it.
2. Ensure `chunk_bytes` stays under the deployment's smallest body cap (16 MiB is
   under the 100 MB Next/Django caps).
3. Schedule `manage.py sweep_staged_uploads` in the maintenance job.
4. That's all — no CORS, no SAS, no direct-to-blob. A future direct-to-blob
   transport can produce the same `staged_upload` handle behind the same `finish`
   contract without touching the import endpoints.

## Tests

- `tests/db/test_staged_upload.py` — the lifecycle and every refusal at the
  library level, owner isolation, consume/no-reuse, sweep.
- `tests/api/unit/test_staged_upload_endpoints.py` — the transport over HTTP.
- `tests/api/unit/test_import_project_by_path.py` — import via a handle, foreign
  handle 404, cloud `local_path` refused.
- `tests/api/unit/test_staging_capability.py` — capability present/absent.
- `tests/db/test_sweep_staged_uploads_command.py` — the sweeper command.
- `client/renderer/__tests__/staged-upload.test.ts` — the client transport.
