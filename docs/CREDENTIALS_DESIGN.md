# Handling Secrets: API Tokens, Passwords and Keys

**Status:** Accepted — implemented for the PDB-REDO case (see *Implementation
status* at the end).
**Date:** 2026-08-21
**Authors:** Martin Noble, Claude

Companion to [preferences-proposal.md](preferences-proposal.md) (which settles
where *non-secret* configuration lives) and
[REMOTE_JOB_EXECUTION_PLAN.md](REMOTE_JOB_EXECUTION_PLAN.md) (which needs the
same machinery for ssh). This document covers the things preferences must *not*
hold: PDB-REDO API tokens, ssh passwords and key passphrases, and — later —
OAuth refresh tokens for facility data access (Diamond/ISPyB).

---

## The problem, concretely

PDB-REDO is the first CCP4i2 task that talks to an authenticated web API, and it
is the case that forces the issue. Today:

```python
# wrappers/pdb_redo_api/script/pdb_redo_api.py
token_id     = str(CCP4Modules.PREFERENCES().PDB_REDO_TOKEN_ID)
token_secret = str(CCP4Modules.PREFERENCES().PDB_REDO_TOKEN_SECRET)
```

Three things are wrong with that, in increasing order of severity:

1. **There is no user interface to set either value.** Not in the task, not in
   Preferences. The only route is hand-editing `~/.ccp4i2-django/preferences.json`.
   The task is effectively unusable out of the box — which is why it has not been
   exercised in alpha.
2. **Unset credentials become the string `"None"`.** `PREFERENCES()` returns
   `None` for an unknown key; `str(None)` is `"None"`, which is then HMAC-signed
   and sent. The user sees an opaque HTTP 401 several minutes into a job rather
   than a clear "you have not configured a token".
3. **At rest the secret is plaintext JSON**, in the same file as `projectsDir`
   and `ccp4Dir` — a file users are routinely asked to paste into support
   threads.

## What the PDB-REDO protocol does and does not require

The token is a pair: a public **token-id** and a private **token-secret**. Each
request is signed AWS-style — `PDBRedoAPIAuth` builds a canonical request,
HMAC-SHA256s it with a key derived from `"PDB-REDO" + token_secret`, and sends
only the signature. **The secret itself never leaves the machine.**

Two consequences follow, and they bound the whole design:

- The signing process must hold the plaintext secret **at the moment of the
  request**. There is no way to keep it off the machine that submits the job. So
  the goal is not "never have the secret in memory" — it is "never have the
  secret *at rest in the clear*, in a log, in the project database, or in an
  exported project".
- Because signing happens in the *job* process, the secret must be reachable
  from a detached `manage.py run_job` subprocess — not only from the GUI.

### Validation

The published API (<https://pdb-redo.eu/api-doc>) documents no token-management
endpoints — only `GET|POST /api/run`, `GET /api/run/{id}`,
`GET /api/run/{id}/output[/{file}|/zipped]`, and `DELETE /api/run/{id}`.

`GET /api/run` is nevertheless a perfectly good validator: it is a signed request
that returns the token's own job list on success and **401** when the token has
been revoked or the secret is wrong. It costs nothing, requires no change at the
PDB-REDO end, and gives us something meaningful to show the user. Per
correspondence with Robbie Joosten (2026-07): expiry is not currently enforced,
but tokens *can* be revoked by PDB-REDO or by the user, so **client-side
validation is worth doing** and a nominal expiry is not something we can read.

---

## Rule zero: a credential is never a task parameter

The tempting shortcut is a `TOKEN_ID` / `TOKEN_SECRET` pair in
`pdb_redo_api.def.xml`, rendered by the generic interface for free.

**This is wrong and must not be done.** Task parameters are persisted into the
job container, written to `input_params.xml` on disk, stored in the project
database, and included verbatim in every exported project zip. A secret placed
there leaks by every one of those routes, silently, and cannot be recalled.

Credentials live in a store that is *orthogonal to projects*: not in the
container, not in the database, not in the project directory. Legacy CCP4i2 got
this right (prefs, not task params) and we keep it.

Corollaries, all of which are testable invariants:

- Credentials never appear in `params_xml`, `diagnostic.xml`, `PROGRAMXML`, or
  any report.
- Credentials are never included in project export/import.
- The API returns credential *descriptors and status*, **never** secret material
  — not even to a loopback client.
- Nothing logs a secret. Validation failures log the status code, not the token.

---

## Where the secret lives at rest

### Decision: an OS keychain, via `keyring`, on the Python side

Checked against `ccp4-20260702`: `keyring` is **not** present, `cryptography`
2.8 and `paramiko` 2.7.1 **are**. `keyring` is pip-installable with pure-Python
platform backends (macOS via `ctypes` to Security.framework, Windows via
`pywin32-ctypes`, Linux via `SecretStorage`/`jeepney`), and the CCP4i2 backend
is already installed into `ccp4-python` by pip — so adding it to
`server/pyproject.toml` is the existing mechanism, not a new one.

**Why the Python side rather than Electron `safeStorage`:**

| | Python `keyring` | Electron `safeStorage` |
|---|---|---|
| GUI works | yes | yes |
| `i2run` / CLI works | yes | **no** — would need env vars only |
| Detached job process can read | yes | only via env plumbing from Electron |
| Number of code paths | one | two (JS vault + Python consumer) |
| Web/cloud deployment | same interface, different backend | irrelevant |

`safeStorage` would make Electron the owner of every secret and force it to
inject them into the Django server's environment at spawn time, which strands
the CLI and makes changing a token require a server restart. One store, read
where it is used, is simpler and strictly more capable.

**Known wart (macOS):** the first keychain read from a given binary prompts the
user for permission. The ACL is per-application, so answering *Always Allow*
once makes subsequent reads by the same `ccp4-python` silent; the prompt may
recur after a CCP4 update changes the binary. As implemented, the job process
reads the keychain directly (twice per job). If that proves irritating in
practice, the mitigation is to read once in the long-lived controlling process
and pass the value to the job in its environment — not done yet, because one
*Always Allow* appears to be enough.

### Fallback where no OS keychain exists

Headless Linux, containers, and CI have no secret service. There the store falls
back to a **mode-`0600` JSON file** at `~/.ccp4i2-django/credentials.json`, kept
deliberately separate from `preferences.json` so the latter stays safe to share.

This fallback is reported to the UI as `secure: false` with a plain-language
explanation. We do **not** encrypt it with a locally-stored key: a Fernet key
sitting next to its own ciphertext is security theatre, and pretending otherwise
is worse than being honest. A `0600` file is roughly what Electron's Linux
`basic_text` backend gives you anyway.

### Resolution precedence

```
environment variable  >  session (in-process)  >  keychain / 0600 file  >  preferences.json (legacy)  >  None
```

- **Env first** preserves the cloud invariance the preferences design depends on:
  containers set `PDB_REDO_TOKEN_SECRET` from Key Vault, and no file or keychain
  is ever consulted.
- **`preferences.json` last** is a read-only compatibility shim for anyone who
  already hand-set `PDB_REDO_TOKEN_ID` there. New writes never go to that file;
  the UI offers to migrate the value into the keychain and remove it.

---

## Three axes, not three systems

The ssh case looks different from the PDB-REDO case, but only along axes that
one store can model. Getting these into the data model on day one is what stops
ssh, and later Diamond, from becoming separate systems.

**1. Kind** — what fields the credential has, and how it is validated.

| Kind | Fields | Example |
|---|---|---|
| `token_pair` | id (public) + secret | PDB-REDO |
| `password` | username + password | ssh password auth |
| `key_file` | path + optional passphrase | ssh key auth |
| `oauth` | access + refresh token, expiry | Diamond/ISPyB (future) |

**2. Scope** — what the credential is attached to: `global` (one PDB-REDO token
for this user), `project` (Robbie's suggestion — a token per CCP4i2 project, so
PDB-REDO's job list stays legible and jobs are attributable), or `target` (this
ssh host, this username). Stored as `(name, scope, scope_id)`; the UI can expose
project-scoped overrides later without a migration.

**3. Persistence** — and this is the axis the ssh case actually needs:

| Mode | Where | Lifetime | For |
|---|---|---|---|
| `keychain` | OS secret store | until cleared | PDB-REDO token, ssh key passphrase |
| `session` | the server process's own environment | until the app quits | an ad-hoc ssh password |
| `none` | nowhere | one request | "don't remember this" |

`session` is implemented by setting the value in the Django process's
`os.environ`, which detached job subprocesses inherit. Nothing touches disk, and
it dies with the process. The remote-execution plan's requirement — "secret sent
once over the POST and never stored" — is exactly `persistence="none"` or
`"session"`, so it needs no separate mechanism.

---

## The API: one choke point

Mirrors the existing `config/program-preferences/` pattern in
[`api/views.py`](../server/ccp4i2/api/views.py), including its `editable` flag,
which is what keeps cloud honest.

| Method | Path | Purpose |
|---|---|---|
| `GET` | `/api/ccp4i2/config/credentials/` | List *descriptors* + status. Never secret material. |
| `POST` | `/api/ccp4i2/config/credentials/{name}/set/` | Store a credential. Loopback-gated. |
| `POST` | `/api/ccp4i2/config/credentials/{name}/validate/` | Server-side signed probe of the live service. |
| `POST` | `/api/ccp4i2/config/credentials/{name}/clear/` | Remove it. |

A descriptor is what the UI needs and nothing more:

```jsonc
{
  "name": "pdb_redo",
  "label": "PDB-REDO web service",
  "signupUrl": "https://pdb-redo.eu/token",
  "fields": [ { "name": "token_id", "label": "Token ID", "secret": false }, … ],
  "isSet": true,
  "source": "keychain",          // or "environment" | "file" | "preferences"
  "secure": true,                // false ⇒ 0600-file fallback; UI says so
  "editable": true,              // false in cloud (env-driven)
  "hint": "…a1b2",               // last 4 of the *public* id only
  "lastValidated": "2026-08-21T09:14:00Z",
  "lastValidationOk": true
}
```

**Write-only secrets.** `set/` takes values in; nothing ever returns them. The
UI can prove a credential is right by pressing *Test*, never by reading it back.

**Desktop gate.** Writes are refused unless `CCP4I2_LOCAL_SESSION_TOKEN` is set
— the signal only the desktop launcher sets, and the same switch `settings.py`
uses to select local-session auth and `set_program_preferences` uses to gate
preference writes. A trusted single-machine server deployment can opt in with
`CCP4I2_ALLOW_CREDENTIAL_WRITE=1`, mirroring the remote-execution plan's
`CCP4I2_ALLOW_INTERACTIVE_SSH`.

Peer address is deliberately **not** consulted. It is tempting to accept
`REMOTE_ADDR in {127.0.0.1, ::1}` as "this machine", but behind a reverse proxy
or a sidecar that is the *proxy's* address, and is `127.0.0.1` for requests that
arrived from the far side of the internet — which would open credential writes
on precisely the multi-tenant deployments the gate exists to protect.

On a shared deployment, credentials arrive as env vars (and, when per-user
credentials are needed there, as encrypted per-user rows — see *Hosted
deployments*).

---

## Validation slots into the existing two-tier architecture

CCP4i2 already distinguishes cheap polled `validity()` from expensive
submit-time `runTimeValidity()` (see CLAUDE.md). Credentials fit exactly:

| Tier | Check | Cost |
|---|---|---|
| `validity()` | Is a token configured at all? | none — a store lookup |
| `runTimeValidity()` | Signed `GET /api/run`: does the service accept it? | one HTTP round trip, once, at submit |

The `validity()` error uses object path `pdb_redo_api.container.inputData.CREDENTIAL_PDB_REDO`.
The `CREDENTIAL_` prefix is a deliberate marker: the run dialog already has a
`QUICK_ACTIONS` table keyed on objectPath suffix, so the same message that says
"no token" can carry the button that fixes it.

The effect is that a missing or revoked token becomes a red indicator *before*
the user presses Run, instead of a 401 twenty minutes later.

---

## The front end

### The error is the entry point

A Preferences → Credentials page alone would be the wrong primary surface, and
building it first would be a mistake. Nobody browses Preferences speculatively;
they discover the need at the moment a task refuses to run. So:

1. **Inline banner in the task.** The PDB-REDO interface shows a compact
   callout at the top when no token is configured: what it is, why it is needed,
   a *Set token…* button, and a link to `pdb-redo.eu/token`. Once the token is
   valid the banner shrinks to a single green line.
2. **Quick action in the Run dialog.** The blocking validation message carries a
   *Set PDB-REDO token…* button, via the existing `QUICK_ACTIONS` mechanism.
   The user fixes it in place and Confirm enables.
3. **Preferences → Credentials is the management surface** — list, re-test,
   clear, see where each value came from. Secondary, not primary.

### The dialog

Generic, driven entirely by the descriptor, so ssh and Diamond reuse it:

- One input per field; `secret: true` renders as a password field with a reveal
  toggle. Existing values are never pre-filled — the field shows a placeholder
  saying a value is stored.
- **A *Test* button is the load-bearing element.** Without it a user cannot
  distinguish wrong-secret from service-down from bad-input-file, and has no way
  to build confidence in something they cannot see. Test runs before Save;
  Save-anyway remains possible with a warning.
- **Paste ergonomics matter more than they sound.** Two long opaque strings,
  copied separately from a web page, is precisely where users fail. The dialog
  accepts a single paste containing both and splits it on whitespace, newline,
  colon, or comma into id and secret.
- **A prominent link to where the token comes from**, opened in the external
  browser (`signupUrl`), because "get an API token first" is a step users have
  usually not taken.
- **Honest storage disclosure**: a line saying where the value will be kept —
  "stored in your macOS Keychain", or "no system secret store available; will be
  written to a file readable only by you".

### What the UI must never do

- Never display, log, or echo a stored secret.
- Never write a secret into a task parameter, a URL, or a query string.
- Never claim encryption the platform is not providing.

---

## Hosted deployments

Keychain storage does nothing for a hosted CCP4i2 (the lab web deployment).
There, one shared token in an env var means every user's PDB-REDO jobs land in
one account — which is exactly the clutter Robbie flagged — *and* means a shared
secret with no attribution.

The correct answer for hosted is **per-user encrypted rows** with the key from
Key Vault, exposed through the same API with `editable: true` and
`source: "database"`. That is a fourth store backend behind an unchanged
interface, deliberately deferred: the desktop alpha does not need it, and the
`(name, scope, scope_id)` key plus the `editable`/`source` fields are what make
it a drop-in later. Until then, hosted deployments report `editable: false` and
take a single token from the environment.

Project-scoped tokens (axis 2) are the cheap partial answer to attribution and
are worth raising with PDB-REDO as a supported pattern.

---

## Implementation status

**Done (this change):**

- `server/ccp4i2/config/credentials.py` — the store: descriptor registry, the
  three axes, keychain/session/file/env backends, resolution precedence, legacy
  `preferences.json` read-through.
- `server/ccp4i2/api/credential_views.py` + routes — list / set / validate /
  clear, loopback-gated, write-only secrets.
- `pdb_redo_api`: reads the store instead of `PREFERENCES()`, gains `validity()`
  and `runTimeValidity()`, and a live validator using signed `GET /api/run`.
  The `str(None)` → `"None"` bug is gone.
- Frontend: generic credential dialog, Preferences → Credentials panel, inline
  banner in the PDB-REDO task, run-dialog quick action.
- `keyring` added to `server/pyproject.toml`.

**Deferred, by design:**

- ssh / remote-execution kinds (`password`, `key_file`) — the store supports
  them; the consumers do not exist yet. See
  [REMOTE_JOB_EXECUTION_PLAN.md](REMOTE_JOB_EXECUTION_PLAN.md).
- Project-scoped credentials — keyed for, not surfaced.
- Per-user encrypted DB backend for hosted deployments.
- OAuth kind for Diamond/ISPyB data access.
