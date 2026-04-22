# Natural-Language Query Proposal (Compounds App)

**Status**: Draft for review. Design captured from dialogue; implementation not yet started.
**Author**: Martin Noble (with Claude)
**Date**: 2026-04-22
**LLM provider**: Azure OpenAI (not Claude API) — same Azure UK South tenant as the rest of the deployment; keeps data in-region and on existing billing.

## 1. Motivation

Medicinal chemists working with the compounds app today compose queries through the aggregation-table UI: pick a project, pick filters, pick columns, export. For many real questions — SAR triage, hit-list construction, ad-hoc "how are my series doing against X assay" lookups — this is fine but costs clicks. A natural-language front door would let a user type:

> *"Tabulate the phys chem properties of all compounds in the ARd project with an IC50 in HTRF less than 10 uM"*

and get the same table without navigating the UI, and with the LLM helping resolve ambiguities (project codenames, protocol name matches, column choices).

The goal is **not** to replace the UI or to generate open-ended analytical prose. The goal is a thin translation layer: **English → structured query spec → existing Django query execution → existing table renderer.** The LLM's job is small and well-bounded (parse intent, propose entity matches, ask for clarification when needed). Everything that actually touches compound data runs in deterministic server-side code.

## 2. Scope

| In scope (v1) | Out of scope (v1, revisit later) |
|---------------|----------------------------------|
| Tabular queries against registered compounds and their assay results | Generative/analytical prose ("summarise the SAR of my series") |
| Entity resolution: target, protocol, metric, units | Writes of any kind (no register / edit / delete via LLM) |
| Any AnalysisResult where `results['kpi_unit']` is populated, `results[KPI]` is numeric, and `status='valid'` — **regardless of `Protocol.import_type`** (Q9 locked; see §17 for evidence) | Rows with missing `kpi_unit` (surfaced as excluded-count in the footer); MS-Intact (no KPI produced); `Protocol.import_type='pharmaron_adme'` or `'table_of_values'` imports made before the import-workflow tightening that lack a `kpi_unit` — curator cleanup |
| Phys-chem columns from `MolecularProperties` | Cross-assay aggregations (best-IC50-across-protocols, etc.) |
| Clarification loop for ambiguous protocol matches | Free-text search of compound comments / notebook entries |
| Backend-side entity resolution (normalized match + distance-1 typo) | Multi-turn conversational memory beyond a single query |
| Rendering into the existing aggregation-table component | A dedicated chat UI (v1 is a command-bar-style input) |
| Per-instance Azure OpenAI binding (demo / kawamura / future) | Streaming token-by-token output |

**v1 is one prompt → one spec → one table.** The dialogue is bounded to clarification only (e.g. "which HTRF protocol did you mean").

## 3. Architectural principle

The LLM is a **parser**, not a query engine. Concretely:

- The LLM **never** sees compound data, assay results, phys-chem values, or anything user-specific beyond the target catalog name strings.
- The LLM emits a **JSON spec** conforming to a schema the backend owns.
- The backend runs **entity resolution** deterministically against the database.
- The backend runs the **query** through the Django ORM; results are rendered through the existing aggregation-table component.
- When resolution is ambiguous, the backend returns a structured "clarify" response and the UI presents a picker. The LLM is not in the loop for clarification turns.

This keeps the LLM's surface area small, makes the system unit-testable end-to-end without an LLM call, and eliminates the class of failures where the LLM fabricates a compound ID or a numeric result.

## 4. End-to-end flow

```
 user prompt (English)
     │
     ▼
 LLM (Azure OpenAI, structured output)
     │
     ▼
 QuerySpec (JSON)
     │
     ▼
 Backend resolver ───► resolve target    ─┐
                  ───► resolve protocol  ─┤── any ambiguity? ──► clarify response ──► UI picker ──► resolved spec
                  ───► resolve metric    ─┘                                                              │
     │                                                                                                   │
     ▼ (all resolved)                                                                                    │
 Backend executor (Django ORM)                                                                           │
     │                                                                                                   │
     ▼                                                                                                   │
 Table payload ──────────────────────► aggregation-table component ◄──────────────────────────────────────
```

The resolver → clarify → resolver loop can iterate more than once (e.g. protocol first, then metric) but each turn is cheap and LLM-free.

## 5. The QuerySpec schema

```jsonc
{
  // Target fields are separate to express cross-project queries
  // ("ARd compounds tested against AKT"). Either can be null.
  "registration_target_as_typed": "ARd",       // Compound.target
  "assay_target_as_typed": "ARd",              // Assay.target
  // target_scope is derived, not emitted — see §6.5

  "protocol_hint": "HTRF",                     // free-text; server token-matches
  "metric": "IC50",                            // must equal results['KPI']
  "threshold": {
    "op": "<",                                 // <, <=, >, >=, =, !=
    "value": 10,
    "unit": "uM"                               // server converts per-row
  },
  "columns": "phys_chem",                      // keyword; expands server-side
  "columns_explicit": null                     // or a list like ["molecular_weight","clogp"]
}
```

The LLM produces this via Azure OpenAI **structured output / JSON schema mode** — not free-form text that we then parse. Fields absent from the prompt are left as nulls/defaults rather than fabricated.

**Two target fields, not one.** The schema separates `registration_target_as_typed` (who designed the compound — `Compound.target`) from `assay_target_as_typed` (what was being tested — `Assay.target`). This lets cross-project queries (*"ARd compounds tested against AKT"*) be expressed without a special syntax, and lets single-target queries set both to the same value. Scope semantics fall out of which fields are populated (§6.5).

## 6. Entity resolution rules

### 6.1 Target (applies to both target fields)

The same resolution rule is applied independently to `registration_target_as_typed` and `assay_target_as_typed`.

- **Data source**: `Target` model in `apps/compounds/registry/models.py`, plus `Gene` metadata when present (see `TARGET_MODEL_PROPOSAL.md`).
- **Matching pool** (when the Target model proposal has landed):
  ```
  normalize(Target.name)
    ∪ {normalize(g.symbol) for g in Target.genes.all()}
    ∪ {normalize(a)         for g in Target.genes.all() for a in g.aliases}
    ∪ {normalize(g.name)    for g in Target.genes.all()}
  ```
  For Targets that haven't yet had gene symbols annotated, the pool degrades gracefully to `{normalize(Target.name)}` only — NLP matching still works, just with fewer hits.
- **Match rule**: normalize both sides (lowercase, strip non-alphanumerics, collapse whitespace) and compare against the pool. On miss, allow Levenshtein distance 1 **only if the normalized query is ≥4 characters** (prevents "AR" matching "AKT"/"AXL"/"ABL"/…).
- **Where**: deterministic, backend-side. The LLM does not see the target list; it emits the string as typed.
- **Ambiguity**: if >1 target matches, return `clarify` with the candidate names. When both target fields are populated and both are ambiguous, clarify one at a time (avoids a combinatorial UI).
- **Miss**: return a structured error with the top 5 nearest matches as suggestions.

**Dependency**: richer matching depends on `TARGET_MODEL_PROPOSAL.md` having landed, but the NLP feature can ship in a degraded form (name-only matching) before that — they are decoupled in rollout order.

### 6.2 Protocol

- **Data source**: `Protocol.name` (CharField, free text, 256 chars). There is **no structured `readout_technology` field** — "HTRF" / "TR-FRET" / "AlphaLISA" / etc. all live as tokens inside the name string.
- **Scope**: protocols that have been run on compounds registered to the resolved target. This narrows a large catalog to a tractable set before matching.
- **Match rule**: tokenize name and query (lowercase, strip punctuation); require all query tokens to appear in the protocol name. Rank ties by (a) fewer extra tokens, (b) recency of last run.
- **Ambiguity is the common case, not an edge.** Multiple `*HTRF*` protocols per target is expected. The clarification picker must be first-class UI (chips with protocol name, number of runs, last-run date, number of compounds tested).
- **Future improvement (not a blocker)**: add a structured `readout_technology` enum to `Protocol` and backfill from name strings. This is a separate migration; v1 does not depend on it.

### 6.3 Metric

- **Data source**: `AnalysisResult.results` (JSONField). Contains a `KPI` string naming the primary metric and one or more numeric values keyed by metric name.
- **Match rule**: the LLM-emitted `metric` must equal `results['KPI']`, and `results[metric]` must be a number, and `status='valid'`. This filters out fits that happen to carry a side-effect value for a different metric than was being targeted.
- **Enumerated metrics today**: IC50, EC50, Ki, Kd, pIC50, % inhibition, DC50, CLint, … (metric is not pre-enumerated in a model; it comes from whatever fitting script produced the result).

### 6.4 Units

**The canonical source is `AnalysisResult.results['kpi_unit']`** — a convention already established in the codebase. The field holds the unit string for whatever `results[KPI]` is, normalized to the form understood by `apps/compounds/assays/kpi_utils.py` (`nM`, `uM`, `mM`, `pM`, `M`, `min`, `s`, `h`, `%`, `uL/min/mg`, `mL/min/kg`, `1e-6 cm/s`, `cm/s`).

- **Data source**: `results['kpi_unit']`. Read directly at query time.
- **Where it comes from**: populated at import time by the fitting/import pipeline; legacy rows can be backfilled with the existing `populate_kpi_units` management command, which walks `DataSeries.dilution_series.unit` → `Protocol.preferred_dilutions.unit` → parsing the KPI field-name suffix (e.g. *"IC50 (nM)"*) in that priority order.
- **Conversion**: per-row — read each row's `kpi_unit`, convert the threshold into that unit once per distinct unit found in the result set, apply the comparison. Use `kpi_utils.normalize_unit()` for canonicalisation.
- **Edge case**: rows where `results['kpi_unit']` is missing or empty are excluded from threshold comparisons with a footer count (*"N rows excluded: no unit recorded"*). This is strictly better than silently comparing numbers of unknown scale.
- **Log-scale metrics (pIC50, etc.)**: handled at metric level, not unit level — if the metric is a `p*` value, the query uses the inverse relation (`pIC50 > 5` ≡ `IC50 < 10 µM`) with the conversion applied in code, not the LLM. The log flag is derivable from the metric name prefix.
- **LLM behaviour unchanged**: LLM still emits `{"value": N, "unit": "..."}` verbatim as the user typed; all unit normalisation and per-row conversion happens server-side.

**Relationship to the DilutionSeries walk from earlier drafts**: that walk is the *backfill* strategy (via `populate_kpi_units`), not the query-time strategy. Query-time always reads `kpi_unit` directly. This was a material simplification over earlier drafts of this proposal.

### 6.5 Target scope (derived from the two target fields)

`Compound.target` (registration, PROTECT, required) and `Assay.target` (experiment target, SET_NULL, nullable) are distinct. Rather than an explicit `target_scope` enum, the scope is **derived** from which of the two target fields are populated and whether they're equal:

| `registration_target_as_typed` | `assay_target_as_typed` | Derived scope | Filter |
|---|---|---|---|
| `"ARd"` | `"ARd"` | Both (default for single-target prompts) | `Compound.target=ARd ∧ Assay.target=ARd` |
| `"ARd"` | `null` | Registered to ARd, any assay target | `Compound.target=ARd` |
| `null` | `"ARd"` | Tested against ARd, any registration | `Assay.target=ARd` |
| `"ARd"` | `"AKT"` | Cross-project: my ARd compounds' AKT potency | `Compound.target=ARd ∧ Assay.target=AKT` |
| `null` | `null` | Error — must have at least one | N/A |

**Default for a single-target prompt** like *"compounds in the ARd project with …"*: LLM emits both target fields equal. This matches the chemist's usual intent (*"my ARd compounds' ARd potency"*). The echo-back names the scope explicitly so the user can override (*"include compounds tested against ARd but registered elsewhere"* → LLM drops `registration_target_as_typed` to `null`).

**"Either" / union queries** (rare in practice) are out of scope for v1. Users who genuinely want the union can run two queries.

**Edge case**: `Assay.target=NULL` is common in older data. Under any scope that filters on `Assay.target`, these rows are excluded with a footer note.

### 6.6 Columns (phys-chem)

- **Data source**: `MolecularProperties` (OneToOne with `Compound`). Fields: `molecular_weight`, `heavy_atom_count`, `hbd`, `hba`, `clogp`, `tpsa`, `rotatable_bonds`, `fraction_sp3`.
- **Default expansion**: `columns: "phys_chem"` → all eight fields.
- **Narrowing vocabulary**: `"lipinski"` → `{MW, HBD, HBA, clogP}`. Specific field names can be listed in `columns_explicit`.
- **Miss**: if `MolecularProperties` does not exist for a compound (SMILES was empty, RDKit failed, etc.), that row shows nulls — not an error.

## 7. Backend components (Django)

Three new pieces under `apps/compounds/` (probably a new `nlp/` app or a module under `assays/`):

1. **`spec.py`** — Pydantic (or `dataclass` + jsonschema) definition of `QuerySpec` and the clarify/error response shapes. Validated on the way in from the LLM, on the way out to the UI.

2. **`resolver.py`** — pure-Python resolution of each spec field against the database. Returns either a fully-resolved `ResolvedQuerySpec` (with `Target` and `Protocol` instances) or a `ClarifyResponse` listing candidates per ambiguous field. No LLM calls. Unit-testable with fixtures only.

3. **`executor.py`** — takes a `ResolvedQuerySpec`, builds the Django ORM query, runs it, shapes the output into the row/column payload the aggregation-table expects. Also computes the footer counts (excluded-for-no-DS, other-scope counts).

Plus one view (or DRF viewset):
- `POST /api/compounds/nlp/query/` — body `{prompt: str}` for fresh queries, or `{spec: QuerySpec, clarifications: {...}}` for continuations. Returns either a `table` payload, a `clarify` payload, or an `error` payload.

## 8. LLM contract

- **Provider**: Azure OpenAI. Model: **`gpt-4o`** during dev/eval — "better rather than cheaper" while the golden set is being built and the failure modes are still being discovered. A downgrade to `gpt-4o-mini` can be revisited post-v1 with eval data showing whether the cheaper model would pass the same bar. Dev/eval is effectively single-user (Martin), so cost during this phase is trivial at either tier.
- **Client library**: the `openai` Python package's `AzureOpenAI` client, configured with an Azure AD token provider.
- **Auth**: **user-assigned managed identity**, reusing the existing `containerAppsIdentity` that already backs Key Vault access, ACR pulls, storage, and queue access (see `Docker/azure-uksouth/infrastructure/infrastructure.bicep:681`). A new role assignment — Cognitive Services OpenAI User (role GUID `5e0bd9bd-7b93-4f28-af87-19fc36ad61bd`) — is added on the Azure OpenAI resource for this identity. **No API keys in env, no secrets to rotate, consistent with every other Azure resource this deployment talks to.**
- **Why not the aacr-abstracts API-key pattern?** aacr-abstracts is a lightweight single-purpose deployment; ccp4i2 is a multi-service production deployment where managed identity is already the project-wide convention for every other secured Azure resource. Matching the aacr-abstracts pattern would introduce a one-off exception to the ccp4i2 security posture. The additional code is small (see §11.1).
- **Mode**: structured output / JSON schema (not free-form; not function-calling). The schema is the `QuerySpec` above minus resolver-only fields.
- **System prompt skeleton**:

  ```
  You translate a user's English request about compounds and assay results
  into a QuerySpec JSON object. You do NOT answer the question; you do NOT
  fetch data; you only parse intent into structured fields.

  Rules:
  - Emit strings verbatim as typed by the user for target_as_typed and
    protocol_hint; do not canonicalise or guess IDs.
  - If the user specifies a unit, echo it verbatim. If not, leave unit as null.
  - Default target_scope to "registered_and_tested" unless the user's phrasing
    clearly indicates otherwise (examples of indicators below).
  - If the request is not a tabular query, emit {"not_a_query": true, "reason": "..."}.
  - Never emit values, IDs, numbers, or compound counts — only the spec.
  ```

- **No catalogs in the prompt** in v1. Target / protocol / metric resolution is backend-side, which keeps the prompt short and avoids token costs that grow with the database.
- **Temperature 0** for determinism.
- **Prompt caching**: system prompt is static ~1k tokens; cache it. Per-call cost is dominated by the user's short prompt.

## 9. Clarification loop

When the resolver returns a `ClarifyResponse`:

```jsonc
{
  "status": "clarify",
  "field": "protocol_hint",
  "question": "Multiple HTRF protocols have been run on AR degraders compounds. Which did you mean?",
  "candidates": [
    {"id": "...", "name": "AR binding HTRF", "n_runs": 12, "last_run": "2026-02-14", "n_compounds": 87},
    {"id": "...", "name": "AR degradation HTRF", "n_runs": 8, "last_run": "2026-03-02", "n_compounds": 54},
    {"id": "...", "name": "AR counter-screen HTRF", "n_runs": 3, "last_run": "2025-11-18", "n_compounds": 31}
  ],
  "partial_spec": { /* the QuerySpec with unresolved field still as hint */ }
}
```

The UI renders a chip picker. User click → client re-POSTs with `{spec: partial_spec, clarifications: {protocol_id: "..."}}`. Backend re-runs resolution with the pinned choice. No LLM call on clarification turns.

## 10. UI surface

- **v1**: a command-bar-style input exposed in **two places**:
  - On the compounds app **landing page** — target must be named in the prompt; enables cross-project queries.
  - On a **per-project page** — target is implicit from the URL route; prompts without a named target inherit the page's target.
- **Rationale for both**: LLM call cost is negligible relative to the rest of the Azure compute budget (see §11), so the per-call frugality argument for "landing page only" or "per-project only" doesn't apply. UX discoverability + cross-project capability justify carrying both surfaces.
- **Placeholder text**: *"Ask a question (e.g. 'phys chem for ARd compounds with HTRF IC50 < 10 µM')"*.
- **Response**: the aggregation-table re-renders with the returned payload. Above it: the echo-back scope sentence (*"Showing compounds registered to AR degraders AND tested in an HTRF protocol against AR degraders, with IC50 < 10 µM"*). Below it: footer counts (excluded-no-DS, counts under other scopes).
- **Clarification**: inline chip picker between input and table, not a modal.
- **Dedicated chat pane is out of scope for v1** — deliberately. Getting the command-bar flow right first lets us learn before adding conversational complexity.

### 10.1 Target context from per-project pages

When the command-bar is used on a project page, the user's prompt may omit the target ("*IC50 < 10 µM for HTRF*"). **Decision: backend post-hoc fill.**

- LLM stays stateless; system prompt is fully static (prompt-caching intact).
- LLM emits `registration_target_as_typed: null` and/or `assay_target_as_typed: null` when no target is named in the prompt.
- Backend checks the request's originating URL route; if it's a per-project page and either target field is null, fills that field with the page's target before entity resolution.
- If the user explicitly names a target in the prompt, that wins — backend does not overwrite a populated field.
- Cross-project phrases like "*also compare against AKT*" from a per-project page won't resolve correctly in v1: the LLM doesn't know there's an implicit primary target, so it can't interpret "also". This is an accepted loss for v1; upgrade to LLM-aware context can be revisited post-v1 if the eval shows this phrasing class matters enough to warrant breaking the static-prompt cache.

**Alternative (a) — LLM-aware context injection — considered and rejected for v1** because: (i) it weakens prompt caching on every call, not just cross-project ones; (ii) the natural-language nuance it buys is narrow and addressable at the landing page where the target would just be named explicitly.

## 11. Azure OpenAI deployment

Reuses the managed-identity pattern already in production across all ccp4i2 container apps (see `Docker/azure-uksouth/infrastructure/infrastructure.bicep:681` for the `containerAppsIdentity` definition, and `applications.bicep` for how server/web/worker consume it).

- **Client library**: `openai` Python package, `AzureOpenAI` client, configured with an Azure AD token provider.
- **Config via non-secret environment variables**:
  - `AZURE_OPENAI_ENDPOINT` — e.g. `https://<resource-name>.openai.azure.com/`
  - `AZURE_OPENAI_MODEL` — deployment name, defaulting to `gpt-4o`
  - `AZURE_OPENAI_API_VERSION` — defaulting to `2024-10-21`
  - (No `AZURE_OPENAI_API_KEY` — auth is via managed identity, below.)
- **Managed identity**: reuses the existing user-assigned identity `containerAppsIdentity`. Needs one new role assignment:
  - **Role**: Cognitive Services OpenAI User
  - **GUID**: `5e0bd9bd-7b93-4f28-af87-19fc36ad61bd`
  - **Scope**: the Azure OpenAI resource
  - Added in the same Bicep module that currently assigns Storage/Queue/KeyVault roles to this identity.
- **Per-instance binding.** Each ccp4i2 Container Apps instance (demo, kawamura, future deployments) can either share one Azure OpenAI resource or have its own, set via `AZURE_OPENAI_ENDPOINT` at deploy time. The identity assignment must exist on whichever resource(s) each instance targets. Sharing one resource across instances is the simpler default; separating per-instance is available if rate-limit isolation or audit separation is wanted later.
- **Region**: UK South (data residency; matches the rest of the deployment).
- **Model choice**: `gpt-4o` during dev/eval (see §8). Downgrade to `-mini` is a post-v1 decision, made with eval data in hand.
- **Rate limiting / quota**: small — v1 is one call per user prompt, temperature 0, short system prompt. Real cost driver is the clarification round-trip count, which is LLM-free.
- **Per-user runaway-loop cap**: 500 LLM calls per user per day. Purpose is to catch retry storms / buggy integrations; far above any plausible human workload. Tripping the cap logs an alert and returns a friendly error for the rest of the day.
- **Disabled in test settings**: tests mock the LLM call; resolver + executor are tested directly.
- **Local development**: `DefaultAzureCredential` falls back to `az login` credentials, so running the Django server locally with `az login` gets the same identity flow without local API keys. If the developer's AAD account doesn't have the OpenAI role, they can either (a) be granted it for dev purposes, or (b) set an optional `AZURE_OPENAI_API_KEY` env var that the client picks up as a fallback, for local-only use.

### 11.1 Client factory

Small shared module, reusable from the NLP view, the eval harness, and any future Azure-OpenAI consumers in ccp4i2:

```python
# apps/compounds/nlp/azure_client.py
import os
from openai import AzureOpenAI
from azure.identity import DefaultAzureCredential, get_bearer_token_provider

_OPENAI_SCOPE = "https://cognitiveservices.azure.com/.default"


def get_azure_openai_client() -> AzureOpenAI:
    endpoint = os.environ["AZURE_OPENAI_ENDPOINT"]
    api_version = os.environ.get("AZURE_OPENAI_API_VERSION", "2024-10-21")

    # Dev-only fallback: if an API key is set, use it. Production paths
    # should not set this var and will fall through to managed identity.
    api_key = os.environ.get("AZURE_OPENAI_API_KEY")
    if api_key:
        return AzureOpenAI(
            api_version=api_version,
            azure_endpoint=endpoint,
            api_key=api_key,
        )

    credential = DefaultAzureCredential()
    token_provider = get_bearer_token_provider(credential, _OPENAI_SCOPE)
    return AzureOpenAI(
        api_version=api_version,
        azure_endpoint=endpoint,
        azure_ad_token_provider=token_provider,
    )


def get_model_name() -> str:
    return os.environ.get("AZURE_OPENAI_MODEL", "gpt-4o")
```

`DefaultAzureCredential` transparently picks up:
- The user-assigned managed identity when running in the Container App (via `AZURE_CLIENT_ID` + the identity binding already set on the app).
- `az login` credentials on a developer machine.
- Workload identity, environment-variable service principals, etc. — standard chain.

Requires adding `azure-identity` to `requirements.txt` (alongside the already-planned `openai` addition).

## 12. Evaluation

- **Golden set**: a YAML file under `apps/compounds/assays/tests/nlp/` with ~30 hand-crafted (prompt, expected_spec) pairs plus ~10 ambiguous prompts and their expected `ClarifyResponse` shape.
- **Test**: for each golden pair, call the LLM (or a mocked deterministic stub in CI) and assert the emitted spec equals expected. Resolver + executor tested independently without the LLM.
- **Success bar for v1**: ≥90% exact match on the golden set at temperature 0; 100% of ambiguous prompts produce a `clarify` response (not a guessed spec).

## 13. Locked decisions (running log)

| # | Decision | Rationale |
|---|----------|-----------|
| 1 | LLM is a parser, not a query engine | Small surface area, unit-testable, no fabrication risk |
| 2 | Backend-side entity resolution, not catalog-in-prompt | Cheaper tokens, deterministic, testable, scales with DB |
| 3 | Target match: normalize + Levenshtein-1 on ≥4 chars | Catches case/whitespace/typos without false hits on short codenames |
| 4 | Protocol match: token-match on `Protocol.name`, scoped to target | No structured readout field exists; ambiguity is the common case |
| 5 | Clarification is UI-mediated (picker), not LLM-mediated | Faster, more reliable, LLM can't garble candidate names |
| 6 | Metric semantics: `KPI == metric ∧ results[KPI] is numeric ∧ status=valid` | Avoids false matches on side-effect values from fits with different KPIs |
| 7 | Units: read `results['kpi_unit']` directly (convention already established in the codebase); convert threshold per-row via `kpi_utils.normalize_unit` | Canonical unit field already exists (see §6.4); earlier DilutionSeries-walk approach is now the backfill strategy (via `populate_kpi_units`), not the query-time strategy |
| 8 | ~~v1 scope: `Protocol.import_type='raw_data'` only~~ — **superseded** — see open question Q9 | Earlier decision assumed units were only recoverable via DilutionSeries; that's wrong — ToV / MS-Intact / ADME imports also populate `kpi_unit` where applicable. v1 scope is being reopened. |
| 9 | Phys-chem default: all 8 `MolecularProperties` fields | Schema-bounded set; narrow via prompt if needed |
| 10 | Target scope default: `registered_and_tested`; echoed back | Matches chemist's usual intent; user can override in the prompt |
| 11 | No writes from LLM (ever) | Minimal trust surface for v1; deliberate |
| 12 | Azure OpenAI auth via the existing **user-assigned managed identity** (`containerAppsIdentity`), with a new Cognitive Services OpenAI User role assignment on the OpenAI resource | Matches ccp4i2's existing pattern (same identity already used for Key Vault / Storage / Queues / ACR); no API-key rotation; single source of truth for container-app permissions. Supersedes an earlier consideration of the `aacr-abstracts` API-key pattern — appropriate there for a thin single-purpose app, not here for a multi-service production deployment. |
| 13 | Structured output / JSON schema mode (not function-calling, not prose) | Task is small and well-bounded; function-calling is overkill for v1 |
| 14 | Command-bar lives on **both** the landing page and per-project pages | LLM call cost is negligible (§11); UX discoverability + cross-project capability justify both surfaces |
| 15 | Two explicit target fields in the spec: `registration_target_as_typed`, `assay_target_as_typed` | Cleanly expresses cross-project queries; `target_scope` becomes derived, not LLM-emitted |
| 16 | Model: `gpt-4o` during dev/eval; downgrade to `-mini` is a post-v1 decision with eval data | "Better rather than cheaper" while failure modes are still being discovered; dev/eval is single-user (Martin) so cost is trivial |
| 17 | Per-user runaway-loop cap: 500 LLM calls/day/user | Purpose is to catch retry storms / buggy integrations; far above any plausible human workload. Not a cost guardrail. |
| 18 | Per-project-page target context: backend post-hoc fill (not LLM-aware injection) | Keeps system prompt fully static and cache-friendly; LLM stays stateless; accepts loss of cross-project "also compare against X" phrasing from per-project pages (user can use landing page instead) |
| 19 | One prompt → one table in v1; multi-table / multi-spec prompts deferred to v2 | Clarification flow across multiple specs is fiddly and worth designing once with real v1 usage data. User can re-ask for the second table. |
| 20 | `not_a_query` UI: show both the LLM-emitted reason and a short list of example prompts | Reason is context-specific (what the user asked for); examples are stable guidance. Combining both is friendlier than either alone. |
| 21 | v1 uses sequential clarification pickers when both target fields are ambiguous; combined-picker option is v1.1 polish | Sequential is simpler and avoids combinatorial UI; the edge case is rare enough that v1 user data will tell us whether a combined picker is worth building. |
| 22 | v1 scope **drops** the `Protocol.import_type='raw_data'` restriction. Any AnalysisResult with numeric `results[KPI]`, `status='valid'`, and valid `results['kpi_unit']` is queryable. Rows failing the `kpi_unit` check are reported as excluded in a footer count. | Audit (see §17) shows raw_data is 100% addressable after `populate_kpi_units`, ToV is 75.5% addressable; excluding ToV entirely would throw away 10,060 genuinely-queryable rows to hide 3,262 "missing unit" ones. Footer is honest and informative. |
| 23 | Import workflows for ToV and Pharmaron ADME will be tightened so new imports cannot create rows with missing `kpi_unit` | Prevents new stuck rows; complements the one-time `populate_kpi_units` run by making the import path correct-by-construction going forward. MS-Intact is skipped because its import does not produce a KPI. |

## 14. Open questions

1. ~~**Where does the command-bar live?**~~ **Resolved (decision 14): both landing page and per-project pages.** See §10.
2. ~~**Exact model choice.**~~ **Resolved (decision 16): `gpt-4o` during dev/eval; `-mini` downgrade is a post-v1 decision with eval data.** See §8.
3. ~~**Cost guardrail**~~ **Resolved (decision 17): 500 LLM calls/day/user as a runaway-loop cap (not cost protection).** See §11.
4. ~~**Does a single prompt ever produce multiple tables?**~~ **Resolved (decision 19): one table per prompt in v1; multi-spec prompts deferred to v2.** See §15.
5. ~~**Are there targets that share a canonical name with different aliases in different projects?**~~ **Resolved: addressed by `TARGET_MODEL_PROPOSAL.md`.** Aliases are hydrated from HGNC via the new `Gene` model and `Target.genes` M2M. NLP matching uses the expanded pool when populated, falls back to `Target.name` only for legacy rows (see §6.1).
6. ~~**What happens when the LLM emits `{"not_a_query": true}`?**~~ **Resolved (decision 20): show both the LLM's reason string and a short list of example prompts.** Exact copy / visual treatment still needs a pass at build time, but the shape is settled.
7. ~~**Per-project page target context — (a) LLM-aware or (b) backend post-hoc fill?**~~ **Resolved (decision 18): (b) backend post-hoc fill.** See §10.1.
8. ~~**When both target fields need clarification, do we clarify sequentially or in one combined picker?**~~ **Resolved (decision 21): sequential for v1; combined-picker revisit is v1.1 polish.** See §15.
9. ~~**(Reopened) Does v1 scope still need the `import_type='raw_data'` restriction?**~~ **Resolved (decision 22): (a) drop the restriction; rows with missing `kpi_unit` surface as an excluded-count in the footer.** Audit evidence in §17 supports this.

## 15. Future work (out of v1)

- ~~Extend scope beyond `raw_data`~~ — potentially already in v1 scope pending Q9 lock (see §14).
- Add `readout_technology` enum to `Protocol`; backfill; replace token-match with attribute filter.
- Cross-assay aggregations: *"best IC50 across any HTRF protocol for each ARd compound"*.
- Multi-table prompts (two metrics, grouped summaries).
- Conversational memory across turns ("now filter to MW < 450") — requires a proper chat UI and per-session state.
- SAR-shaped queries that return structured summaries rather than tables ("what's the MW trend as IC50 improves in my ARd series?").
- Write-capable commands (register, annotate, flag) — only after a clear trust/audit model is in place.

## 16. Migration / rollout

- Ship behind a feature flag / env var (`COMPOUNDS_NLP_ENABLED`) so it can be enabled per-instance independently of the code deploy.
- Demo instance first, with eval golden set passing.
- Kawamura instance second, with instance-specific targets/protocols exercised in the golden set.
- No database migration needed for v1 (purely additive, read-only).

## 17. Appendix: audit evidence for Q9 decision

The `audit_kpi_unit_coverage` management command was run against the DDU Database production instance (ddudatabase.ncl.ac.uk) on 2026-04-22, before and after a one-off `populate_kpi_units` run.

### 17.1 Before `populate_kpi_units`

| import_type | total | valid | +num | +unit | +valid | ready-of-`+num` |
|---|---:|---:|---:|---:|---:|---:|
| `raw_data` | 21,848 | 11,502 | 11,499 | 10,347 | 10,347 | **90.0%** |
| `table_of_values` | 17,429 | 16,884 | 13,322 | 10,060 | 10,060 | **75.5%** |
| `ms_intact` | 1 | 1 | 1 | 1 | 1 | 100% |
| `__no_protocol__` (orphans) | 16 | 16 | 0 | 0 | 0 | — |

- raw_data: all 1,152 missing-unit rows recoverable via DilutionSeries walk.
- ToV: 3,262 missing-unit rows, **0 recoverable** — their KPI strings (e.g. `'EC50'`, `'SPR solubility'`, `'Modification'`) have no unit suffix to parse and no DilutionSeries linkage.
- ms_intact: single row, not statistically meaningful; MS-Intact imports do not produce a KPI string in general (see §17.3).

### 17.2 After `populate_kpi_units`

| import_type | total | valid | +num | +unit | +valid | ready-of-`+num` |
|---|---:|---:|---:|---:|---:|---:|
| `raw_data` | 21,848 | 11,502 | 11,499 | 11,499 | 11,499 | **100.0%** |
| `table_of_values` | 17,429 | 16,884 | 13,322 | 10,060 | 10,060 | **75.5%** |

- raw_data: every valid, numeric-KPI row is now fully queryable.
- ToV: unchanged, as predicted — recovery required information that was never captured at import time.

### 17.3 Per-bucket KPI distributions (top entries)

**`raw_data`** (all clean, dose-response): EC50 (10,347), IC50 (969), KI (184).

**`table_of_values`**: mixture of unit-suffixed keys (which are clean) and bare keys (which are the stuck 3,262):

- `'KD (M)'` (2,868) — unit in name → clean
- `'EC50'` (2,456) — **no unit → stuck**
- `'EC50 (uM)'` (1,264), `'A431 IC50(nM) '` (756), `'H1975 IC50(nM)'` (495) — unit in name → clean
- `'SPR solubility'` (577) — **no unit → stuck**
- `'Modification'` (493) — **not numeric, outside NLP scope**
- `'EC50 (M)'` (475), `'OVCAR3 GI50 (µM)'` (456), `'WT GI50(nM)'` (449) — clean

The lesson: the import-time parser catches unit suffixes when they exist; the 3,262 stuck rows correspond to spreadsheet columns whose headers lacked a suffix and where no override was supplied. This is exactly the failure mode the §23 import-workflow tightening addresses going forward.

### 17.4 What this means for Q9

- v1 scope drops the `Protocol.import_type='raw_data'` restriction: raw_data is fully solved, ToV is 75.5%-solvable with the remaining 25% honestly surfaced as excluded-count in the footer.
- The alternative (`raw_data`-only) would throw away 10,060 genuinely queryable ToV rows to avoid showing 3,262 exclusions. Not a trade that pays off.
- The import-workflow tightening (§13 decision 23) prevents new stuck rows; `populate_kpi_units` handles what's recoverable in existing data; curator intervention is the only path for the residual 3,262 stuck ToV rows, and it's optional (they'll simply be footer-counted).
