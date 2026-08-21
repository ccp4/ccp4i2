# Writing a Task That Calls an Authenticated Web Service

How to add a task that submits work to a remote API needing a token, password or
key — the `pdb_redo_api` shape. Companion to
[PHIL_TASK_GUIDE.md](PHIL_TASK_GUIDE.md) and
[EXPORT_TASK_GUIDE.md](EXPORT_TASK_GUIDE.md); the reasoning behind the
credential store itself is in [docs/CREDENTIALS_DESIGN.md](../../../docs/CREDENTIALS_DESIGN.md).

Worked example throughout: `wrappers/pdb_redo_api/`.

---

## Rule zero: a credential is never a task parameter

Do **not** put a token, password or passphrase in your `def.xml`, however
convenient it looks. Task parameters are persisted into the job container,
written to `input_params.xml` on disk, stored in the project database, and
included verbatim in every **exported project zip**. A secret placed there leaks
by all four routes, silently, and cannot be recalled.

Credentials live in a store that is orthogonal to projects:
`ccp4i2/config/credentials.py`.

---

## Step 1 — Register the credential

Add a `CredentialSpec` to `CREDENTIALS` in
[`config/credentials.py`](../config/credentials.py). That is the whole
registration: the REST API and the entire frontend are descriptor-driven, so no
UI code is needed for a new credential.

```python
MY_SERVICE = CredentialSpec(
    name="my_service",                 # used in URLs and storage keys
    label="My Service",
    description="One or two sentences shown in the dialog.",
    fields=(
        CredentialField(name="token_id", label="Token ID", secret=False,
                        env="MY_SERVICE_TOKEN_ID"),
        CredentialField(name="token_secret", label="Token secret", secret=True,
                        env="MY_SERVICE_TOKEN_SECRET"),
    ),
    signup_url="https://example.org/token",
    validator=_validate_my_service,    # see step 2
    kind="token_pair",                 # token_pair | password | key_file | oauth
    paste_split=True,                  # accept both values in one paste
    hint_field="token_id",             # non-secret field the UI may hint from
)

CREDENTIALS = {PDB_REDO.name: PDB_REDO, MY_SERVICE.name: MY_SERVICE}
```

Notes on the fields that are easy to get wrong:

| Field | Why it matters |
|---|---|
| `secret=True` | The value is never returned by the API, never rendered, never logged. Mark anything that is not *public* as secret. |
| `env` | Also the environment-variable name. This is how cloud deployments supply the value (env beats every other layer), and how CI runs your task. |
| `hint_field` | Must be a **non-secret** field. The UI shows its last 4 characters so a user can tell two credentials apart. |
| `paste_split` | Set it whenever there is more than one field. Copying two long opaque strings separately from a web page is where users actually fail. |

## Step 2 — Write a validator

A validator turns "is this credential any good?" into one cheap call. It takes
the resolved `{field: value}` dict and returns `(ok, message)`.

```python
def _validate_my_service(values):
    import requests
    r = requests.get("https://example.org/api/whoami",
                     auth=_auth(values), timeout=20)
    if r.status_code == 200:
        return True, "Credential accepted."
    if r.status_code in (401, 403):
        return False, "The service rejected this credential."
    return False, f"The service returned HTTP {r.status_code}."
```

Three requirements:

- **Distinguish rejection from unreachability.** A user told "token rejected"
  when the network is down will hunt for a typo that is not there. Catch
  `requests.RequestException` separately and say so.
- **Never put secret material in the message.** It is displayed and may be
  logged. The status code is the whole story.
- **Any endpoint will do.** Services rarely publish a "validate token" call.
  PDB-REDO has none, so we use `GET /api/run` — the caller's own job list, which
  200s for a good token and 401s for a revoked one, and costs nothing.

Import the service's auth helper *lazily* inside the function so
`config/credentials.py` stays importable without the wrapper tree.

## Step 3 — Read it in the wrapper

```python
from ccp4i2.config import credentials

def _token(self):
    return credentials.get_credential("my_service")
```

`get_credential` returns `None` unless **every** field resolved — never a
half-filled dict. That matters: signing a request with half a token pair
produces an opaque HTTP 401 several minutes into a job, which is exactly the
failure this API exists to prevent.

Do **not** reach for `CCP4Modules.PREFERENCES()`. It returns `None` for an unset
key, so the once-idiomatic `str(PREFERENCES().MY_TOKEN)` yields the literal
string `"None"`, which then gets signed and sent.

Resolution order, for reference:

```
env var  >  session  >  OS keychain / 0600 file  >  preferences.json (legacy)  >  None
```

## Step 4 — Validate in both tiers

CCP4i2 polls `validity()` while the user edits and calls `runTimeValidity()`
once at submission. Credentials use both:

```python
CREDENTIAL_ERROR_PATH = 'CREDENTIAL_MY_SERVICE'   # note the prefix

def validity(self):
    """Cheap, polled. No network."""
    error = super().validity()
    if self._token() is None:
        error.append(
            klass=self.TASKNAME, code=200,
            details='My Service needs a token before it can submit a job. ...',
            name=f'{self.TASKNAME}.container.inputData.{CREDENTIAL_ERROR_PATH}',
            severity=CCP4ErrorHandling.SEVERITY_ERROR)
    return error

def runTimeValidity(self):
    """Expensive, once, at submission."""
    error = super().runTimeValidity()
    if error.maxSeverity() >= CCP4ErrorHandling.SEVERITY_ERROR:
        return error                      # already failing; skip the round trip
    ok, message = credentials.validate_credential('my_service')
    if not ok:
        error.append(klass=self.TASKNAME, code=201, details=message,
                     name=f'{self.TASKNAME}.container.inputData.{CREDENTIAL_ERROR_PATH}',
                     severity=CCP4ErrorHandling.SEVERITY_ERROR)
    return error
```

> **The `CREDENTIAL_` prefix is a contract with the frontend.** The run dialog
> (`client/renderer/providers/run-check-provider.tsx`) matches any objectPath
> containing `.CREDENTIAL_<NAME>` and renders a *"Set token…"* button beside the
> message, opening the dialog for `<name>` (lowercased). Use the prefix and the
> error that blocks the run carries its own fix. The path does not need to
> correspond to a real container item.

Together these turn "opaque 401 twenty minutes in" into a red indicator before
the user presses Run.

## Step 5 — Surface it in the task interface

Add the banner to your React interface. It is descriptor-driven, so this is the
whole of the frontend work:

```tsx
import { CredentialAlert } from "../task-elements/credential-alert";

<CredentialAlert name="my_service" onChanged={() => mutateValidation?.()} />
```

It renders a *"Set token…"* prompt when nothing is configured, an error when the
service last rejected the credential, and a single quiet line with a *Test*
button when all is well.

This is the **primary** surface. Preferences → Credentials exists, but nobody
browses Preferences speculatively — users meet the need when a task will not
run, so the fix belongs where the problem appears.

---

## Polling a long-running remote job

If your task submits and waits, the waiting loop is where things go wrong.
`wrappers/pdb_redo_api/script/test_api.py` is the reference implementation.

**Be a good citizen.** A flat 5-second poll is ~720 requests an hour *per job*
against someone else's server. Start responsive and back off:

```python
POLL_INITIAL_SECONDS = 5      # short runs still feel immediate
POLL_MAX_SECONDS = 60         # ~65 requests in the first hour, not 720
POLL_BACKOFF = 1.5
```

Reset the interval when the reported status changes — transitions cluster, and
that is when responsiveness is worth paying for.

**A blip must not fail the job.** Over a multi-hour run on a laptop, a 502, a
dropped connection or a sleeping lid is likely. Retry transient failures and
require *consecutive* failures before giving up; a single success resets the
count. Always pass `timeout=` to `requests`, or a hung connection wedges the job
forever.

**Print on change, not on every poll.** Plus a heartbeat every ~10 minutes, so a
long wait does not look like a hang and the job log stays readable.

## Failure taxonomy

Distinguish the ways waiting can end, because a bare `except:` cannot tell the
user anything useful — and because the three cases need three different actions:

| Exception | Means | What the user should do |
|---|---|---|
| `…AuthError` | Credentials rejected mid-run (401/403) | Set a working credential; **do not retry** — it will not fix itself |
| `…JobStopped` | The remote run failed at their end | Read the logs — salvage the results archive, which usually contains them |
| `…Unreachable` | We lost contact | **Nothing is lost remotely.** The run is probably still going |

That last row is the one that matters most. When the local job fails because we
gave up waiting, the remote run is very likely fine — so **carry the remote run
id on the exception** and put it, plus the URL of the service's job list, in the
failure message. Being told "look for run 4213 at example.org/jobs" is the
difference between a lost afternoon and a click.

Report failures through both channels:

```python
def _fail(self, code, message):
    print("ERROR:", message)          # lands in the job log, in context
    sys.stdout.flush()
    self.appendErrorReport(code, message,
                           severity=CCP4ErrorHandling.SEVERITY_ERROR)
    return CPluginScript.FAILED
```

## Never

- Put a credential in a `def.xml`, a URL, a query string, or a log line.
- Write a credential into `PROGRAMXML`, `diagnostic.xml`, or a report.
- Return secret material from an API endpoint — not even to a local client.
- Report "rejected" for what was actually a network failure.

## Testing

Credential-dependent tasks must be testable without a credential and without the
network:

- **Store behaviour** — `tests/unit/config/test_credentials.py`. Stub
  `credentials._keyring` to `lambda: None` so tests never touch the developer's
  real keychain.
- **Validators** — stub `requests.get`; assert that rejection, unreachability
  and success produce distinguishable messages, and that no message contains the
  secret.
- **Polling** — `tests/unit/plugins/test_pdb_redo_monitor.py` drives the loop
  from a scripted list of responses and exceptions with `time.sleep` patched
  out. Assert that a transient blip is survived and that backoff is capped.
- **The task itself** — an i2run test can supply the credential from the `env`
  names, which is why every field must declare one.

## Checklist

- [ ] `CredentialSpec` added to `CREDENTIALS`, every field with an `env` name
- [ ] Validator that separates rejection from unreachability and leaks nothing
- [ ] Wrapper reads via `credentials.get_credential()`, never `PREFERENCES()`
- [ ] `validity()` (is it set) and `runTimeValidity()` (does it work)
- [ ] Error `name` uses the `.CREDENTIAL_<NAME>` marker
- [ ] `<CredentialAlert name="…">` in the task interface
- [ ] Polling backs off, tolerates transient failures, passes `timeout=`
- [ ] Typed exceptions carrying the remote run id; no bare `except:`
- [ ] Unit tests that need neither a credential nor the network
