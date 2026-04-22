#!/usr/bin/env python3
"""
Interactive backfill of Target.genes against a deployed ccp4i2 instance.

Runs locally on the operator's machine; talks to the REST API over HTTPS.
Authentication uses the Azure CLI's cached token for the instance's AAD
app registration. No container exec, no web UI, no secrets on disk.

Auth follows the same conventions as `Docker/cli/i2remote.py`:

  - If `CCP4I2_API_URL` + `CCP4I2_API_TOKEN` are set in the environment, those
    are used directly (easiest: `i2remote login` once, then every subsequent
    CLI call just works).
  - Otherwise the script falls back to deriving the URL from CUSTOM_DOMAIN
    in the instance env file, and fetches a token by shelling out to
    `az account get-access-token --resource api://<AZURE_AD_CLIENT_ID>`.

Usage
-----

    # Option 1: use an existing i2remote login
    i2remote login
    export CCP4I2_API_URL=https://ddudatabase.ncl.ac.uk
    export CCP4I2_API_TOKEN=$(jq -r .api_token ~/.ccp4i2remote.json)
    ./scripts/backfill_target_genes.py --dump plan.yaml

    # Option 2: use az CLI + an instance env file
    az login
    ./scripts/backfill_target_genes.py \\
        --env Docker/azure-uksouth/.env.deployment \\
        --dump plan.yaml

    # Edit plan.yaml — set `action:` on each target to either
    #   - a list of gene symbols (will be applied)
    #   - the string "skip" (target left unchanged)

    # Apply
    ./scripts/backfill_target_genes.py \\
        --env Docker/azure-uksouth/.env.deployment \\
        --apply plan.yaml

Flags
-----

    --env <path>        .env file for the target instance (optional if
                        CCP4I2_API_URL + CCP4I2_API_TOKEN are set)
    --dump <path>       stage 1: write plan file
    --apply <path>      stage 2: read plan file and POST updates
    --base-url <url>    override URL (or CCP4I2_API_URL)
    --client-id <guid>  override AAD client id from env file
    --include-genned    include targets that already have genes in the plan
    --candidates <n>    number of HGNC candidates per target (default: 5)

Requirements
------------

    - Python 3.9+
    - PyYAML  (pip install pyyaml)
    - Either: a valid CCP4I2_API_TOKEN already in env
      Or:     Azure CLI (`az login`) with access to the instance's AAD group
"""

from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys
from pathlib import Path
from typing import Any, Optional
from urllib.error import HTTPError, URLError
from urllib.parse import quote, urljoin
from urllib.request import Request, urlopen

try:
    import yaml  # type: ignore
except ImportError:
    sys.stderr.write(
        "ERROR: PyYAML is required. Install with: pip install pyyaml\n"
    )
    sys.exit(2)


HGNC_REST_BASE = "https://rest.genenames.org"
DEFAULT_TIMEOUT = 15

# Frontend (Next.js) proxy prefix for compounds-app REST calls. The public
# host exposes the web tier; authenticated requests to /api/proxy/compounds/<...>
# are forwarded to the Django backend at /api/compounds/<...>. Hitting the
# direct Django path through the public host returns the Next.js SPA shell
# or redirects to /auth/login.
COMPOUNDS_PROXY_PREFIX = "/api/proxy/compounds"


# ---------------------------------------------------------------- #
# env file parsing
# ---------------------------------------------------------------- #


def parse_env_file(path: Path) -> dict[str, str]:
    """Minimal .env parser: KEY=value lines, ignores comments and blanks."""
    env: dict[str, str] = {}
    for line_no, raw in enumerate(path.read_text().splitlines(), start=1):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        if "=" not in line:
            continue
        key, _, value = line.partition("=")
        key = key.strip()
        value = value.strip().strip('"').strip("'")
        env[key] = value
    return env


def derive_base_url(env: dict[str, str], override: Optional[str]) -> str:
    """Precedence: --base-url > CCP4I2_API_URL env > CUSTOM_DOMAIN in env file."""
    if override:
        return override.rstrip("/")
    from_env = os.environ.get("CCP4I2_API_URL", "").strip()
    if from_env:
        return from_env.rstrip("/")
    custom = env.get("CUSTOM_DOMAIN") or ""
    if custom:
        return f"https://{custom}".rstrip("/")
    raise SystemExit(
        "ERROR: cannot determine base URL — set CCP4I2_API_URL, "
        "pass --base-url, or provide an --env file with CUSTOM_DOMAIN."
    )


def derive_client_id(env: dict[str, str], override: Optional[str]) -> str:
    """Precedence: --client-id > AZURE_CLIENT_ID env > AZURE_AD_CLIENT_ID in env file."""
    if override:
        return override
    from_env = os.environ.get("AZURE_CLIENT_ID", "").strip()
    if from_env:
        return from_env
    cid = env.get("AZURE_AD_CLIENT_ID", "").strip()
    if not cid:
        raise SystemExit(
            "ERROR: cannot determine AAD client id — set AZURE_CLIENT_ID, "
            "pass --client-id, or provide an --env file with AZURE_AD_CLIENT_ID."
        )
    return cid


# ---------------------------------------------------------------- #
# auth
# ---------------------------------------------------------------- #


def get_bearer_token(client_id: Optional[str]) -> str:
    """
    Resolve a bearer token. Precedence:
      1. CCP4I2_API_TOKEN env var (i2remote convention — set after `i2remote login`)
      2. ~/.ccp4i2remote.json's api_token (i2remote config file)
      3. Shell out to `az account get-access-token --resource api://<client_id>`
    """
    # (1) env var
    env_token = os.environ.get("CCP4I2_API_TOKEN", "").strip()
    if env_token:
        return env_token

    # (2) i2remote config
    cfg = Path.home() / ".ccp4i2remote.json"
    if cfg.exists():
        try:
            data = json.loads(cfg.read_text())
            disk_token = (data.get("api_token") or "").strip()
            if disk_token:
                return disk_token
        except (json.JSONDecodeError, OSError):
            pass  # fall through

    # (3) az CLI
    if not client_id:
        raise SystemExit(
            "ERROR: no token available — set CCP4I2_API_TOKEN, run "
            "`i2remote login`, or provide --client-id / --env to fetch via az."
        )
    resource = f"api://{client_id}"
    try:
        raw = subprocess.check_output(
            ["az", "account", "get-access-token", "--resource", resource],
            stderr=subprocess.PIPE,
        )
    except FileNotFoundError:
        raise SystemExit(
            "ERROR: `az` CLI not found. Install Azure CLI and run `az login`, "
            "or run `i2remote login` and set CCP4I2_API_TOKEN."
        )
    except subprocess.CalledProcessError as exc:
        stderr = exc.stderr.decode() if exc.stderr else ""
        hint = ""
        if "AADSTS9002313" in stderr or "malformed" in stderr:
            hint = (
                "\nHint: AADSTS9002313 usually means the token-acquisition "
                "scope hasn't been consented for this AAD app in your az "
                "session. Easiest workaround: run `i2remote login` (which "
                "uses MSAL directly) and re-run this script with "
                "CCP4I2_API_TOKEN set from the login output."
            )
        raise SystemExit(
            f"ERROR: `az account get-access-token --resource {resource}` failed.\n"
            f"Are you logged in (`az login`)?\n{stderr}{hint}"
        )
    return json.loads(raw)["accessToken"]


# ---------------------------------------------------------------- #
# REST client (thin)
# ---------------------------------------------------------------- #


def api_get(base_url: str, path: str, token: str) -> Any:
    # NB: strip trailing slashes. The Next.js proxy appends one before
    # forwarding to Django; if we send one too, Django responds 308 →
    # redirect URL *without* trailing slash, and urllib silently drops
    # the Authorization header on the redirect (security default), giving
    # us an HTML login page on the second hop. See Docker/cli/i2remote.py.
    path = path.rstrip("/")
    url = urljoin(base_url + "/", path.lstrip("/"))
    req = Request(
        url,
        headers={
            "Authorization": f"Bearer {token}",
            "Accept": "application/json",
        },
    )
    try:
        with urlopen(req, timeout=DEFAULT_TIMEOUT) as resp:
            return _decode_json(url, resp.status, resp.headers, resp.read())
    except HTTPError as exc:
        body = exc.read().decode("utf-8", errors="replace")[:800]
        raise SystemExit(
            f"GET {url} -> HTTP {exc.code}\n"
            f"Content-Type: {exc.headers.get('Content-Type')}\n"
            f"Body (first 800 chars):\n{body}"
        )
    except URLError as exc:
        raise SystemExit(f"GET {url} network error: {exc.reason}")


def api_patch(base_url: str, path: str, token: str, body: dict) -> Any:
    path = path.rstrip("/")  # see note in api_get
    url = urljoin(base_url + "/", path.lstrip("/"))
    data = json.dumps(body).encode("utf-8")
    req = Request(
        url,
        data=data,
        headers={
            "Authorization": f"Bearer {token}",
            "Accept": "application/json",
            "Content-Type": "application/json",
        },
        method="PATCH",
    )
    try:
        with urlopen(req, timeout=DEFAULT_TIMEOUT) as resp:
            return _decode_json(url, resp.status, resp.headers, resp.read())
    except HTTPError as exc:
        body_text = exc.read().decode("utf-8", errors="replace")[:800]
        raise SystemExit(
            f"PATCH {url} -> HTTP {exc.code}\n"
            f"Content-Type: {exc.headers.get('Content-Type')}\n"
            f"Body (first 800 chars):\n{body_text}"
        )
    except URLError as exc:
        raise SystemExit(f"PATCH {url} network error: {exc.reason}")


def _decode_json(url, status, headers, raw_bytes):
    """Decode a JSON response body, or raise with diagnostic context."""
    body_text = raw_bytes.decode("utf-8", errors="replace")
    content_type = headers.get("Content-Type", "")

    # A 2xx with HTML means we hit an auth redirect / wrong gateway.
    if "json" not in content_type.lower():
        hint = ""
        if "text/html" in content_type.lower():
            hint = (
                "\nHint: HTML response on a JSON endpoint usually means the "
                "token was not accepted (wrong audience, expired, or the "
                "instance's AAD guard rejected it). Quick diagnostics:\n"
                "  curl -sS -D- -H \"Authorization: Bearer $CCP4I2_API_TOKEN\" "
                f"-H \"Accept: application/json\" \"{url}\" | head -20\n"
                "  # Decode your token's aud / groups claims:\n"
                "  python3 -c 'import base64, json, os, sys; "
                "t=os.environ[\"CCP4I2_API_TOKEN\"].split(\".\"); "
                "p=t[1]+\"=\"*(-len(t[1])%4); "
                "print(json.dumps(json.loads(base64.urlsafe_b64decode(p)), indent=2))'"
            )
        raise SystemExit(
            f"GET/PATCH {url} -> HTTP {status}\n"
            f"Content-Type: {content_type or '(unset)'}\n"
            f"Body (first 800 chars):\n{body_text[:800]}{hint}"
        )

    try:
        return json.loads(body_text)
    except json.JSONDecodeError as exc:
        raise SystemExit(
            f"GET/PATCH {url} -> HTTP {status} returned non-JSON body\n"
            f"Content-Type: {content_type}\n"
            f"JSON error: {exc}\n"
            f"Body (first 800 chars):\n{body_text[:800]}"
        )


def list_all(base_url: str, path: str, token: str) -> list[dict]:
    """Fetch a paginated DRF list, returning the flat result list."""
    out: list[dict] = []
    next_url: Optional[str] = path
    while next_url:
        payload = api_get(base_url, next_url, token)
        if isinstance(payload, list):
            return payload  # unpaginated
        out.extend(payload.get("results", []))
        next_url = payload.get("next")
        if next_url and next_url.startswith(base_url):
            next_url = next_url[len(base_url):]  # make relative
    return out


# ---------------------------------------------------------------- #
# HGNC client (candidate search)
# ---------------------------------------------------------------- #


_SYMBOL_LIKE = re.compile(r"\b[A-Z][A-Z0-9]{1,7}\b")


def hgnc_search(query: str, limit: int = 10) -> list[dict]:
    """Text-search HGNC; return the top `limit` docs (symbol, hgnc_id, name…)."""
    if not query.strip():
        return []
    url = f"{HGNC_REST_BASE}/search/{quote(query.strip())}"
    req = Request(url, headers={"Accept": "application/json"})
    try:
        with urlopen(req, timeout=DEFAULT_TIMEOUT) as resp:
            payload = json.loads(resp.read().decode("utf-8"))
    except (HTTPError, URLError):
        return []
    docs = payload.get("response", {}).get("docs", [])
    return docs[:limit]


def hgnc_fetch_by_symbol(symbol: str) -> Optional[dict]:
    url = f"{HGNC_REST_BASE}/fetch/symbol/{quote(symbol.strip().upper())}"
    req = Request(url, headers={"Accept": "application/json"})
    try:
        with urlopen(req, timeout=DEFAULT_TIMEOUT) as resp:
            payload = json.loads(resp.read().decode("utf-8"))
    except (HTTPError, URLError):
        return None
    docs = payload.get("response", {}).get("docs", [])
    return docs[0] if docs else None


def candidates_for_target(target_name: str, limit: int) -> list[dict]:
    """
    Build a ranked list of candidate gene symbols for a target.

    Strategy (v1, HGNC-only, no LLM):
      (1) Direct HGNC search on the whole name.
      (2) Extract symbol-like ALL-CAPS tokens from the name; exact-fetch each.
      (3) Deduplicate by symbol, preserving first-seen order.
    """
    seen: set[str] = set()
    results: list[dict] = []

    # (1) text search
    for doc in hgnc_search(target_name, limit=limit):
        sym = (doc.get("symbol") or "").upper()
        if sym and sym not in seen:
            seen.add(sym)
            results.append({
                "symbol": sym,
                "hgnc_id": doc.get("hgnc_id", ""),
                "name": doc.get("name", ""),
                "source": "hgnc-search",
            })

    # (2) token lookups
    for token in _SYMBOL_LIKE.findall(target_name.upper()):
        if token in seen:
            continue
        rec = hgnc_fetch_by_symbol(token)
        if rec:
            seen.add(token)
            results.append({
                "symbol": token,
                "hgnc_id": rec.get("hgnc_id", ""),
                "name": rec.get("name", ""),
                "source": "hgnc-fetch-by-token",
            })

    return results[:limit]


# ---------------------------------------------------------------- #
# dump / apply stages
# ---------------------------------------------------------------- #


def cmd_dump(args, base_url: str, token: str):
    targets = list_all(base_url, f"{COMPOUNDS_PROXY_PREFIX}/targets", token)

    if not args.include_genned:
        targets = [t for t in targets if not t.get("genes")]

    if not targets:
        print("No targets need backfill.", file=sys.stderr)
        return

    plan = {
        "instance": base_url,
        "schema_version": 1,
        "_usage": (
            "For each target, set `action:` to EITHER a list of symbols "
            "like ['EGFR','ERBB2'], OR the string 'skip' to leave it alone. "
            "Candidates are hints; edit freely. Run with --apply when done."
        ),
        "targets": [
            {
                "id": t["id"],
                "name": t["name"],
                "current_genes": [g["symbol"] for g in t.get("genes", [])],
                "candidates": candidates_for_target(t["name"], args.candidates),
                "action": "?",  # user fills this in
            }
            for t in targets
        ],
    }

    out_path = Path(args.dump)
    out_path.write_text(yaml.safe_dump(plan, sort_keys=False, allow_unicode=True))
    print(
        f"Wrote plan with {len(plan['targets'])} target(s) to {out_path}",
        file=sys.stderr,
    )


def cmd_apply(args, base_url: str, token: str):
    path = Path(args.apply)
    plan = yaml.safe_load(path.read_text())
    if not isinstance(plan, dict) or "targets" not in plan:
        raise SystemExit(f"ERROR: {path} is not a valid plan file.")

    # Safety: warn if plan was generated against a different instance.
    plan_instance = plan.get("instance")
    if plan_instance and plan_instance != base_url:
        sys.stderr.write(
            f"WARNING: plan was generated against {plan_instance} but you are "
            f"applying to {base_url}. Proceed only if you're sure "
            f"(Ctrl-C to abort, Enter to continue). "
        )
        input()

    applied = skipped = unresolved = errors = 0
    for entry in plan.get("targets", []):
        tid = entry.get("id")
        name = entry.get("name", "?")
        action = entry.get("action")

        if action == "skip" or action is None:
            skipped += 1
            continue
        if action == "?":
            unresolved += 1
            sys.stderr.write(f"  [{name}] unresolved (action: '?') — skipped\n")
            continue
        if not isinstance(action, list) or not all(isinstance(s, str) for s in action):
            errors += 1
            sys.stderr.write(
                f"  [{name}] ERROR: action must be a list of symbol strings "
                f"or 'skip', got {action!r}\n"
            )
            continue

        try:
            api_patch(
                base_url,
                f"{COMPOUNDS_PROXY_PREFIX}/targets/{tid}",
                token,
                {"gene_symbols": action},
            )
            sys.stderr.write(f"  [{name}] -> {action}\n")
            applied += 1
        except SystemExit as exc:
            errors += 1
            sys.stderr.write(f"  [{name}] ERROR: {exc}\n")

    print("", file=sys.stderr)
    print(
        f"Applied: {applied}  Skipped: {skipped}  Unresolved: {unresolved}  Errors: {errors}",
        file=sys.stderr,
    )


# ---------------------------------------------------------------- #
# main
# ---------------------------------------------------------------- #


def main():
    parser = argparse.ArgumentParser(
        description="Interactive backfill of Target.genes via REST API.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("--env", help=".env file for the target instance "
                                      "(optional if CCP4I2_API_URL + CCP4I2_API_TOKEN are set)")
    parser.add_argument("--base-url", help="Override base URL from env file / CCP4I2_API_URL")
    parser.add_argument("--client-id", help="Override AAD client id from env file / AZURE_CLIENT_ID")
    parser.add_argument("--include-genned", action="store_true",
                        help="Include targets that already have genes in the plan")
    parser.add_argument("--candidates", type=int, default=5,
                        help="Number of HGNC candidates per target (default: 5)")

    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--dump", metavar="PATH", help="Stage 1: write plan file")
    mode.add_argument("--apply", metavar="PATH", help="Stage 2: read plan and apply")

    args = parser.parse_args()

    env: dict[str, str] = {}
    if args.env:
        env_path = Path(args.env)
        if not env_path.exists():
            raise SystemExit(f"ERROR: env file not found: {env_path}")
        env = parse_env_file(env_path)

    base_url = derive_base_url(env, args.base_url)
    # client_id only needed if we end up falling back to az for a token.
    client_id: Optional[str] = None
    try:
        client_id = derive_client_id(env, args.client_id)
    except SystemExit:
        # No client id; that's fine if CCP4I2_API_TOKEN is set. We'll
        # re-raise inside get_bearer_token if it's actually needed.
        pass

    sys.stderr.write(f"Instance: {base_url}\n")
    if client_id:
        sys.stderr.write(f"AAD client id (for token fallback): {client_id}\n")

    token = get_bearer_token(client_id)

    if args.dump:
        cmd_dump(args, base_url, token)
    else:
        cmd_apply(args, base_url, token)


if __name__ == "__main__":
    main()
