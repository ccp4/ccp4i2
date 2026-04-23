# Natural-Language Query Proposal (Compounds App)

**Status**: In implementation. Backend slices 1–5 shipped (target resolution, protocol resolution + two-target orchestrator, row-level evaluator with unit tri-state, end-to-end executor with phys-chem column expansion + scope-sentence generation, Azure OpenAI client + parse_prompt with structured-output extraction). DRF view (§7 view, §9), UI (§10), and evaluation harness (§12) remain. **See §18 for slice-by-slice progress and the pick-up notes for the next slice** — that section is the canonical resumption anchor if this work is being resumed in a fresh conversation.
**Author**: Martin Noble (with Claude)
**Date**: 2026-04-22 (drafted); §§18–19 added 2026-04-23; §18 slices 4–5 updates 2026-04-23
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

- **Data source**: `Target` model in `apps/compounds/registry/models.py`, plus `Gene` metadata (see `TARGET_MODEL_PROPOSAL.md` — shipped).
- **Matching pool**:
  ```
  normalize(Target.name)
    ∪ {normalize(g.symbol) for g in Target.genes.all()}
    ∪ {normalize(a)         for g in Target.genes.all() for a in g.aliases}
    ∪ {normalize(g.name)    for g in Target.genes.all()}
  ```
  Targets with no gene symbols annotated (legacy rows, non-gene modality tags, genuinely ambiguous programmes left un-genned during backfill) fall back to `{normalize(Target.name)}` only — still matches on the programme name verbatim.
- **Match rule**: normalize both sides (lowercase, strip non-alphanumerics, collapse whitespace) and compare against the pool. On miss, allow Levenshtein distance 1 **only if the normalized query is ≥4 characters** (prevents "AR" matching "AKT"/"AXL"/"ABL"/…).
- **Where**: deterministic, backend-side. The LLM does not see the target list; it emits the string as typed.
- **Ambiguity**: if >1 target matches, return `clarify` with the candidate names. When both target fields are populated and both are ambiguous, clarify one at a time (avoids a combinatorial UI).
- **Miss**: return a structured error with the top 5 nearest matches as suggestions.

**Implications of hydrated aliases (post-ship affordances):**

- User can type any of a target's gene symbols, its aliases, its prev-symbols, or its canonical HGNC name and hit the same target. *"SKP2"*, *"CKS1B"*, and *"Skp2-Cks1"* all resolve to the Skp2-Cks1 PPI programme.
- Pan-family queries ("*Aurora kinase*", "*cyclin D*") resolve naturally via the multi-gene Target.genes links (Myc-Aur carries AURKA/AURKB; CDK4 carries CCND1/CCND2/CCND3).
- **Ambiguity is expected to rise** — typing "MYC" now hits both Myc-Aur and Myc RNA. That's the clarification picker's happy path (§9, sequential); no new code needed. Worth calling out in eval so we don't treat it as a regression.

**No longer a dependency, now a prerequisite that has landed.** The rollout order in §16 reflects this.

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

**The canonical source is `AnalysisResult.results['kpi_unit']`** — a convention already established in the codebase. The field holds the unit string for whatever `results[KPI]` is, normalized to the form understood by `apps/compounds/assays/kpi_utils.py`:

- **Concentrations**: `nM`, `uM`, `mM`, `pM`, `M`
- **Time**: `min`, `s`, `h`
- **Fraction**: `%`
- **Clearance / permeability**: `uL/min/mg`, `mL/min/kg`, `1e-6 cm/s`, `cm/s`
- **`unitless`** — the explicit "no unit applies" sentinel for dimensionless KPIs (efflux ratio, Hill coefficient, Fsp3, fold-shifts, any ratio or fraction). *Distinct from unknown/missing.*

Three states matter when reading a row's `kpi_unit`:

| State | `results['kpi_unit']` | Behaviour under a threshold query |
|---|---|---|
| Known real unit | `"nM"`, `"uM"`, `"%"`, etc. | Compare directly (unit-match) or convert (unit-compatible). Include in results. |
| Deliberately unit-less | `"unitless"` (or, for backward compat, `None`) | Compare directly *only if the query supplied no unit*. Query with a real unit → mismatch → row excluded with "unit-type mismatch" footer note. |
| Unknown / missing | Absent key, empty string | Exclude from threshold comparisons with a footer count (*"N rows excluded: no unit recorded"*). Strictly better than silently comparing numbers of unknown scale. |

- **Data source**: `results['kpi_unit']`. Read directly at query time.
- **Where it comes from**: populated at import time by the fitting/import pipeline; legacy rows can be backfilled with the existing `populate_kpi_units` management command, which walks `DataSeries.dilution_series.unit` → `Protocol.preferred_dilutions.unit` → parsing the KPI field-name suffix (e.g. *"IC50 (nM)"*) in that priority order. Unit-less KPIs (Caco-2 efflux ratio etc.) are set explicitly by their importers to `"unitless"`; they never need backfill guessing.
- **Conversion**: per-row — read each row's `kpi_unit`, convert the threshold into that unit once per distinct unit found in the result set, apply the comparison. Use `kpi_utils.normalize_unit()` for canonicalisation.
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
| 24 | Resolver uses the Gene-hydrated matching pool from day one (symbol + aliases + prev_symbols + HGNC name, per-target). Name-only is the fallback for rows without genes. | `TARGET_MODEL_PROPOSAL.md` has shipped: Gene model, Target.genes M2M, HGNC hydration, TargetSerializer read/write of genes, and a backfill script. The DDU instance's targets are hydrated. NLP v1 uses this matching surface from first deploy rather than waiting for a v1.1. |
| 25 | `kpi_unit` is tri-state: known real unit, `"unitless"` sentinel, or unknown (absent/empty). Threshold comparison treats each state differently: known → direct/convert; unitless → only matches unit-less queries; unknown → excluded with footer count. | Caco-2 efflux ratio, Hill coefficient, Fsp3 and other dimensionless KPIs are legitimate query targets. Conflating them with "import was sloppy" (the earlier decision) would silently drop all Caco-2 rows. The explicit sentinel distinguishes the two. |

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

The bullets below are the short list of v1 carve-outs. **See §19 for a longer-form treatment of v2+ directions** — including chemistry-aware substructure queries and the architectural framing that unifies them.

- ~~Extend scope beyond `raw_data`~~ — potentially already in v1 scope pending Q9 lock (see §14).
- Add `readout_technology` enum to `Protocol`; backfill; replace token-match with attribute filter.
- Cross-assay aggregations: *"best IC50 across any HTRF protocol for each ARd compound"*.
- Multi-table prompts (two metrics, grouped summaries).
- Conversational memory across turns ("now filter to MW < 450") — requires a proper chat UI and per-session state.
- SAR-shaped queries that return structured summaries rather than tables ("what's the MW trend as IC50 improves in my ARd series?").
- Write-capable commands (register, annotate, flag) — only after a clear trust/audit model is in place.

## 16. Migration / rollout

**Step 0 (already done):** the Target-model dependency and hydration pipeline are shipped. `Gene` model + `Target.genes` M2M + `hydrate_genes` command + serializer read/write of genes + `backfill_target_genes.py` are deployed to all three instances. DDU Database is fully backfilled (32 of 44 targets linked to gene_symbols; 1 deferred (GLK — needs disambiguation); 11 marked `skip` as non-gene programmes). Demo and Kawamura have the code but haven't yet been backfilled.

**Remaining NLP-v1 rollout:**

- Ship NLP v1 behind a feature flag / env var (`COMPOUNDS_NLP_ENABLED`) so it can be enabled per-instance independently of the code deploy.
- Demo instance first, with eval golden set passing.
- Kawamura instance second, with instance-specific targets/protocols exercised in the golden set.
- No database migration needed for NLP v1 (purely additive, read-only).
- Before turning it on for Demo/Kawamura: run `backfill_target_genes.py --dump` on those instances and curate the plans, so the matching pool is as rich there as on DDU.

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

## 18. Implementation progress

This section tracks what has been built, what remains, and the decisions made *during* implementation that are not captured in §1–17. **If you are picking up this work in a new conversation, read §1–17 for the design, then read this section for state-of-play.** The memory file `reference_nlp_query_proposals.md` points here; treat this as the canonical progress index, not conversation memory.

### 18.1 Slice plan

Work is delivered as narrow vertical slices that build the backend first and defer the LLM to last. This keeps the resolver + executor fully unit-testable against deterministic fixtures before any Azure OpenAI call enters the picture — the principle of §3.

| Slice | Scope | Status | Proposal refs |
|---|---|---|---|
| 1 | Target field resolution | Shipped | §5 (spec), §6.1, decision 24 |
| 2 | Protocol resolution + two-target orchestrator | Shipped | §6.2, §6.5, Q21 |
| 3 | Row-level evaluator (metric + unit tri-state) | Shipped | §6.3, §6.4, Q25 |
| 4 | Executor (`executor.py`): ORM query + payload + footer counts + scope sentence | Shipped | §7, §6.6, §10 |
| 5 | Azure OpenAI client + prompt + structured-output extraction | Shipped | §8, §11.1 |
| 6 | DRF view with clarify loop | Pending | §7 (view), §9 |
| 7 | Frontend command-bar (landing + per-project) | Pending | §10 |
| 8 | Evaluation harness + golden set | Pending | §12 |

Rollout dependencies: §16 prerequisites (Gene model + DDU hydration) are done. Demo + Kawamura backfill happens after the executor lands but before the feature flag is turned on for those instances.

### 18.2 Landed slices — public API surface

All code lives under `apps/compounds/nlp/`. Tests live under `apps/compounds/nlp/tests/`. Test counts are cumulative.

**Slice 1 — Target resolution** (28 tests)

- Modules: `apps/compounds/nlp/spec.py`, `apps/compounds/nlp/resolver.py`
- Public functions:
  - `normalize(s) -> str` — lowercase + remove all non-alphanumerics
  - `resolve_target(as_typed) -> TargetResolution`
- Return type is a discriminated union:
  - `ResolvedTarget(target, matched_via, query)` — `target` is the ORM instance; `matched_via` ∈ {`name`, `gene_symbol`, `gene_alias`, `gene_name`} plus `*_fuzzy` suffixes for Levenshtein-1 hits
  - `TargetClarify(query, candidates: list[TargetCandidate], field)` — more than one target matched
  - `TargetMiss(query, suggestions: list[TargetCandidate], field)` — no match; top-5 `difflib.SequenceMatcher`-ranked
- `TargetCandidate(id, name, gene_symbols)` — thin JSON-serialisable view used for clarify/miss payloads
- Tests: `apps/compounds/nlp/tests/test_resolver_target.py`. Fixture mirrors DDU shape (single-gene target, PPI with two genes, shared-gene MYC ambiguity across two targets, name-only legacy target).

**Slice 2 — Protocol resolution + two-target orchestrator** (59 tests cumulative)

- Additions to `spec.py`: `ResolvedTargets`, `ScopeError`, `ProtocolCandidate`, `ResolvedProtocol`, `ProtocolClarify`, `ProtocolMiss`; string constants `SCOPE_BOTH_SAME / SCOPE_REG_ONLY / SCOPE_ASSAY_ONLY / SCOPE_CROSS` and `FIELD_REGISTRATION_TARGET / FIELD_ASSAY_TARGET / FIELD_PROTOCOL_HINT`.
- Public functions (in `resolver.py`):
  - `resolve_targets(spec) -> ResolvedTargets | TargetClarify | TargetMiss | ScopeError` — §6.5 scope derivation; Q21 sequential clarify (registration field first, then assay); returns `ScopeError` when both fields empty.
  - `tokenize(s) -> set[str]` — lowercase + split on runs of non-alphanumerics.
  - `resolve_protocol(hint, resolved_targets) -> ProtocolResolution` — scope-filtered pool via `Assay.target` + `DataSeries.compound.target`; set-containment token match; Clarify ranks by (fewer extra tokens, then recency of last run); Miss ranks by token overlap, top 5.
- Tests: `apps/compounds/nlp/tests/test_resolver_orchestrator.py`, `apps/compounds/nlp/tests/test_resolver_protocol.py`. Protocol fixture includes an explicit cross-project selectivity case (an AR compound tested against AKT) to exercise the reg=X/assay=Y scope combination.

**Slice 3 — Row evaluator** (129 tests cumulative)

- Additions to `spec.py`: `RowMatched / RowFiltered / RowExcluded / RowOutcome` union. Reason constants are split into two groups:
  - `FILTER_*` — silent, not footer-noted: `invalid_status`, `kpi_mismatch`, `value_not_numeric`, `threshold_not_met`
  - `EXCLUDE_*` — footer-noted per §6.4: `unit_unknown`, `unit_type_mismatch`, `unit_incompatible`, `query_missing_unit`
- New module: `apps/compounds/nlp/evaluator.py`
- Public functions:
  - `convert(value, from_unit, to_unit) -> Optional[float]` — within-family conversion using a `_UNIT_FAMILIES` factor table (concentration, time, percent, permeability, two separate clearance families). Cross-family returns `None`.
  - `classify_row_unit(kpi_unit) -> str` — Q25 tri-state: a real-unit string, `ROW_UNIT_UNITLESS`, or `ROW_UNIT_UNKNOWN`. Delegates to `compounds.assays.kpi_utils.is_unitless / normalize_unit / VALID_UNITS`.
  - `classify_query_unit(unit) -> str` — a real-unit string, `QUERY_UNIT_NONE`, or `QUERY_UNIT_UNKNOWN`. A query-typed `"unitless"` is treated as `QUERY_UNIT_NONE` — users don't type "unitless" in natural language.
  - `evaluate_row(row, metric, threshold) -> RowOutcome` — duck-typed over anything with `.status` + `.results` (tests use `SimpleNamespace` stand-ins; no DB needed).
- Tests: `apps/compounds/nlp/tests/test_evaluator.py` — exhaustive matrix across row-unit class × query-unit class × threshold op, plus conversion correctness.

**Slice 4 — Executor** (148 tests cumulative)

- Additions to `spec.py`: `TableRow`, `TablePayload`, `SpecError`, `ExecutionResult` union; column-preset constants (`COL_PRESET_PHYS_CHEM`, `COL_PRESET_LIPINSKI`) and canonical ordered tuples `PHYS_CHEM_COLUMNS` (8 fields) and `LIPINSKI_COLUMNS` (4 fields).
- New module: `apps/compounds/nlp/executor.py`
- Public function: `execute(spec: QuerySpec) -> ExecutionResult`. Composes `resolve_targets` → `resolve_protocol` → scoped `AnalysisResult` walk → `evaluate_row` per row. Matched rows become `TableRow`s; `EXCLUDE_*` outcomes aggregate into `footer_excluded`; `FILTER_*` outcomes aggregate into `filtered_silent` (exposed for debug/telemetry, not surfaced in UI per §18.2).
- Scope queryset: `AnalysisResult.objects.filter(status='valid', data_series__assay__protocol=<resolved>, ...)` with optional `data_series__compound__target` and `data_series__assay__target` filters; `select_related('data_series__compound__molecular_properties')` for phys-chem plucking.
- Scope sentence generation for §10 echo-back is per-`scope_kind` via f-strings in `_scope_sentence`; threshold clause (`, with IC50 < 10.0 nM`) or metric clause (`, IC50 values`) appended at the end.
- Tests: `apps/compounds/nlp/tests/test_executor.py` (19 tests). Covers: spec validation (missing metric / protocol_hint → `SpecError`); pass-through of `ScopeError` / `TargetMiss` / `TargetClarify` / `ProtocolMiss` / `ProtocolClarify`; happy-path row emission with phys-chem values and unit conversion; footer + silent-filtered bucket counts; column expansion (default / lipinski / explicit / explicit-with-unknowns); scope sentence across all four `scope_kind`s.

**Slice 5 — LLM parse** (165 tests cumulative)

- New `apps/compounds/nlp/azure_client.py`:
  - `get_azure_openai_client()` factory implementing §11.1 exactly — `DefaultAzureCredential` → managed identity in prod, `az login` locally. Reads `AZURE_OPENAI_API_KEY` first as the documented dev-only fallback.
  - Deferred imports: `openai` and `azure.identity` are imported *inside* the factory, not at module top. Keeps test runs free of the SDK requirement (tests mock the factory), and means the rest of `compounds.nlp` is safe to import in environments where `openai` isn't pip-installed yet.
  - `get_model_name()` reads `AZURE_OPENAI_MODEL`, defaulting to `"gpt-4o"` per decision 16.
- New `apps/compounds/nlp/llm.py`:
  - `SYSTEM_PROMPT` — the static system prompt. Encodes the load-bearing rules (verbatim strings, no IDs/values fabricated, single-target vs cross-project convention, `not_a_query` emission, filler-noun discipline from §19.6). Static so prompt caching works.
  - `PROMPT_SCHEMA` — the JSON Schema given to Azure OpenAI in strict-mode structured-output. Every property is nullable and listed in `required` (strict-mode requirement). `threshold` is a nested object schema with `op/value/unit`.
  - `parse_prompt(prompt: str) -> PromptParseResult` — calls the SDK with `response_format={"type": "json_schema", …, "strict": True}` and `temperature=0`; parses the JSON into a `QuerySpec` or `NotAQuery`. Handles malformed JSON, non-object JSON, and empty content as `ParseError`.
  - `_to_parse_result(data: dict)` — pure translator, exported for direct testing without any mocking.
- `spec.py` additions: `NotAQuery(reason)`, `ParseError(message, raw)`, and the `PromptParseResult = Union[QuerySpec, NotAQuery, ParseError]` union. `NotAQuery` is user-facing (the LLM's decision); `ParseError` is infrastructure-facing (our glue code failed).
- Tests: `apps/compounds/nlp/tests/test_llm.py` (17 tests). Pure translator branches (single-target/threshold, assay-only, cross-project, unitless threshold, not-a-query, malformed-threshold, columns_explicit preservation) plus `parse_prompt` paths with a mocked `AzureOpenAI` client (happy path asserting full SDK call shape, not_a_query, malformed-JSON, non-object-JSON, empty-content, empty-prompt short-circuit), plus schema-drift guards (every property is in `required`; threshold subschema locks `additionalProperties=False`) and a prompt content spot-check.
- Requirements: `openai>=1.50,<2.0` added to `Docker/server/Dockerfile` and documented in `server/requirements-azure.txt`. `azure-identity` was already pinned for the rest of the Azure stack.

### 18.3 Decisions made during implementation (additions to §13)

These are lived-in details that aren't in §1–17. They don't contradict §13 — they extend it.

| Slice | Decision | Rationale |
|---|---|---|
| 1 | `normalize(s)` **removes** non-alphanumerics entirely (does not replace with whitespace) | Makes "Skp2-Cks1", "skp2 cks1", and "skp2cks1" all collapse to the same string. "AR degraders" matches "ardegraders" and "ar_degraders". Verified against §6.1's implied behaviour via tests. Implementation is the shared `compounds.utils.normalize_ref`, reused by both compound-ref resolution and NLP target resolution. |
| 1 | `ResolvedTarget.target` holds the Django model instance directly; `TargetClarify / TargetMiss` candidates use the thin `TargetCandidate` | The downstream executor (slice 4) needs the ORM instance; outbound JSON payloads need the trimmed view. No reason to pick one over both. |
| 1 | Target matching pool is built per call (no cache/inverted index) | DDU has ~44 targets; a full iteration is trivial. Promote to a cached index only if scale demands. |
| 2 | `scope_kind` is a string constant rather than an `Enum` | Cleaner JSON serialisation; matches §6.5 naming conventions. Same applies to `matched_via`, filter/exclude reasons, field-name tags. |
| 2 | The orchestrator tags `TargetClarify.field / TargetMiss.field` by passing the field name into `dataclasses.replace` at construction | Set once, immutably, as part of the returned result — no post-hoc mutation of the returned dataclass. |
| 2 | `tokenize` uses **set** semantics (duplicate tokens in a protocol name are collapsed) | Matches §6.2's "all query tokens appear in the protocol name" literally. Repeated tokens carry no weight. |
| 2 | Scope join for registration-target filter is `data_series__compound__target` | An assay with no data_series is not considered "run on any compound" and is correctly filtered out. Cross-project selectivity (AR compound tested against AKT) is expressible as `reg=AR, assay=AKT`. |
| 3 | New exclusion reason `EXCLUDE_QUERY_MISSING_UNIT` for the case "query has no unit but row has a real unit" | §6.4's table did not cover this cell explicitly. Footer-excluding keeps the feature working in heterogeneous scopes (some rows unitless, some real-unit); the UI can aggregate the count and surface "your threshold needs a unit". A query-level hard reject is also defensible — flagged here so slice 6 (view) can revisit if eval feedback warrants. |
| 3 | An **absent** `kpi_unit` key resolves to `ROW_UNIT_UNITLESS` (because `is_unitless(None) == True` per `kpi_utils`), **not** `ROW_UNIT_UNKNOWN` | Inherited from the backward-compat rule baked into `compounds.assays.kpi_utils.is_unitless`. Consequence: a row with no `kpi_unit` key under a real-unit threshold query → `EXCLUDE_UNIT_TYPE_MISMATCH`, not `EXCLUDE_UNIT_UNKNOWN`. Subtle but consistent with how legacy importers wrote ratio KPIs. An **empty** `kpi_unit` string (`""`) does resolve to `ROW_UNIT_UNKNOWN` — that is what §6.4's "missing" column means. |
| 3 | `evaluate_row` is duck-typed over anything with `.status` + `.results` | Lets evaluator tests use `SimpleNamespace` stand-ins; no DB required, suite runs in ~2s. The real executor will hand `AnalysisResult` instances in at call time. |
| 3 | Two clearance families (`uL/min/mg` microsomal, `mL/min/kg` hepatic) stand alone — they are not interconvertible | Different denominators (mg microsomal protein vs kg body weight) aren't scalar-related; silent conversion would mask a category error. Same-unit comparison only. |
| 3 | Log-scale (`pIC50` ↔ `IC50`) rewriting is deferred | §6.4 mentions it but doesn't mandate v1 support. A row with `KPI='pIC50'` matches a `pIC50` query directly. Revisit as a dedicated slice if eval shows users interchangeably write both forms. |
| 4 | Executor emits **per-measurement** rows, not per-compound aggregations | The evaluator is row-level; emitting matching rows directly gives slice 7 maximum flexibility. Aggregation to per-compound (geomean etc.) can be done by the frontend using existing aggregation helpers, or by a future v1.1 payload variant if chemists prefer that shape. |
| 4 | Scope queryset pre-filters `status='valid'` at the DB level, short-circuiting `FILTER_INVALID_STATUS` | `FILTER_INVALID_STATUS` is silent (not footer-noted) so pre-filtering is observationally equivalent and avoids hauling invalid rows back to Python. The evaluator still classifies it for safety/tests; it just won't fire in normal flow. |
| 4 | `FILTER_*` reason counts are still exposed on `TablePayload.filtered_silent` | Not surfaced in UI, but useful for debug/telemetry ("how many rows did we toss for KPI mismatch?"). Slice 6 can choose whether to log / return / hide it. |
| 4 | `columns_explicit` silently drops unknown field names; fully-unknown → phys-chem default | Prompt → spec → executor is one-way; hard errors on mid-pipeline column names would produce brittle UX. Unknown names usually mean the LLM invented something; phys-chem fallback keeps the query productive. Spec-level validation (missing `metric`, missing `protocol_hint`) *is* a hard `SpecError` because those are load-bearing. |
| 4 | `SpecError` is a distinct response type from `ScopeError` | `ScopeError` is specific to the "neither target field set" violation of §6.5. `SpecError(field, message)` is a generic executor-level validation error. Keeping them separate makes the view's error-to-HTTP mapping (slice 6) cleaner. |
| 4 | Scope-sentence rendering uses literal f-strings per `scope_kind`, not a template library | Four scope kinds × two threshold states = 8 variants, all short. A template engine would add indirection without saving code. `, with {metric} {op} {value} {unit}` or `, {metric} values` appended at the end. |
| 5 | `NotAQuery` (user-facing) and `ParseError` (infra-facing) are separate response types | `NotAQuery` is a considered decision by the LLM ("this isn't a table") and should surface to the user with the reason + example prompts per decision 20. `ParseError` is *our* glue code failing (bad JSON, empty content, wrong shape). Slice 6 view will map them to different HTTP statuses (200-with-explanation vs 502/500). |
| 5 | `openai` / `azure.identity` imports are deferred to inside `get_azure_openai_client()` | Lets `compounds.nlp.llm` and tests import without requiring the SDK installed locally. Tests mock the factory, so even CI with no `openai` on the interpreter path passes cleanly. Production Dockerfile pulls the SDK in. |
| 5 | Strict-mode JSON schema with every property in `required` and nullable union types (`["string", "null"]`) | Azure OpenAI strict-mode structured output requires `additionalProperties=false` and `required` listing every property. Optional fields are expressed as `["T", "null"]` rather than `required` omission. A `test_schema_declares_every_spec_field_as_required` guard catches drift if fields are added without updating `required`. |
| 5 | `not_a_query` is a boolean discriminator on the flat schema, not an `anyOf` union | Keeps the schema strict-mode-compliant and the JSON parse straightforward — a single object with optional fields, `not_a_query=true` selecting the NotAQuery branch. `anyOf` works in structured output but adds complexity the v1 shape doesn't need. Revisit if the schema grows a second-level branch. |
| 5 | Empty/whitespace prompt short-circuits to `ParseError` before calling the SDK | Cheap guard against accidental empty-input burning LLM budget. Tests verify the factory isn't called. |
| 5 | `_to_parse_result` is exported as a public-for-testing helper | Lets tests exhaustively cover the translator without repeatedly mocking the HTTP path. Ten of the seventeen tests hit the translator directly. |

### 18.4 Running the tests

From `server/`, with CCP4 sourced:

```bash
source ../../ccp4-20251105/bin/ccp4.setup-sh
PYTHONPATH="$PWD:$PWD/../apps" DJANGO_SETTINGS_MODULE=compounds.settings \
  ccp4-python -m pytest ../apps/compounds/nlp/tests/ -v
```

Expected: 165 passing as of slice 5 (2026-04-23). Evaluator (slice 3) and LLM (slice 5) tests do not touch the DB — LLM tests mock the Azure OpenAI factory so they also don't require `openai` to be installed locally. Resolver / orchestrator / protocol / executor tests use `@pytest.mark.django_db` via the `db` fixture argument. Fixtures are defined inline in each test module (no JSON fixture files, no separate conftest).

### 18.5 Pending slices — pick-up notes

**Slice 6: DRF view + clarify loop** — the next slice.

- `POST /api/compounds/nlp/query/` — body `{prompt}` for fresh queries, `{spec, clarifications}` for continuations. Continuations must NOT re-call the LLM (see decision 5) — they resume from a prior partial_spec + the UI's chosen clarification value (e.g. a pinned target_id or protocol_id).
- Response union: `TablePayload` (200) | clarify payload derived from `TargetClarify` / `ProtocolClarify` (200 + status discriminator) | error payload derived from `TargetMiss` / `ProtocolMiss` / `ScopeError` / `SpecError` / `NotAQuery` / `ParseError`. `NotAQuery` shows the LLM's reason + example prompts (decision 20). `ParseError` is a 502-ish "LLM output couldn't be parsed — retry?".
- Per-project-page target context (decision 18) is a backend post-hoc fill on null target fields; the request needs to carry the originating URL route or an explicit `project_id` / `project_target` hint.
- Hook in the §11 per-user 500-calls-per-day runaway cap. Counter should be per-user, reset at UTC midnight, and fail closed with a friendly 429.
- **Infra prerequisite** (flag before deploying): the Bicep `applications.bicep` module needs a new role assignment of "Cognitive Services OpenAI User" (GUID `5e0bd9bd-7b93-4f28-af87-19fc36ad61bd`) on the target Azure OpenAI resource, scoped to the existing `containerAppsIdentity`. Without it, `parse_prompt` will fail at token-acquisition in prod.

**Slice 7: Frontend command-bar.**

- Two surfaces (compounds-app landing page and per-project page) per decision 14.
- Inline chip picker for clarify, not a modal.
- Response rendered through the existing aggregation-table component.

**Slice 8: Evaluation.**

- Golden set: `apps/compounds/assays/tests/nlp/*.yaml` with (prompt, expected_spec) pairs and (prompt, expected_clarify_shape) pairs for ambiguous cases.
- Success bar per §12: ≥90% exact spec match at T=0; 100% of ambiguous prompts produce a clarify rather than a guessed spec.

### 18.6 Deliberately out of scope mid-build

These are tracked in §15 as future work but worth calling out so they're not accidentally pulled into slices 4–7:

- Log-scale metric rewriting — deferred pending eval signal.
- Multi-table prompts (decision 19 defers to v2).
- Conversational memory across turns.
- Cross-assay aggregations (best-IC50-across-protocols).
- Write-capable commands (decision 11 locks "never" for v1).

## 19. Direction sketches for v2+

Forward-looking thinking about where this interface could go after v1 ships. §15 is the terse list of v1 carve-outs; this section captures the **architectural consequences** of various v2 directions — which ones fit the v1 shape and which need a distinct layer. Not a commitment to build; a map of where the shape could grow. Preserved here so the thinking survives any session boundary.

### 19.1 What the v1 architecture buys us

The v1 shape — LLM-as-parser → deterministic resolver → ORM — is a **tabular-query door**. It generalises cleanly to anything expressible as "structured query against a known schema with human-friendly input." It does **not** generalise to generative analysis (SAR summaries, novel hypothesis generation, cross-row narrative) without adding a separate layer with different data-visibility rules for the LLM. Keeping those two directions distinct is the load-bearing discipline — conflating them is how interfaces like this turn into poorly-specified chat-over-chemistry.

### 19.2 Direction buckets

| Bucket | What it is | Example prompts |
|---|---|---|
| Deeper query shapes | Same architecture, richer spec — new `QuerySpec` fields, same resolver/executor pattern | *"ARd compounds with no HTRF measurement yet"*; *"compounds where ARd IC50 is 10× better than AKT IC50"*; *"IC50 and MW and CLint for ARd compounds"*; *"compounds registered in the last 6 months"* |
| Adjacent capabilities | Same rails, different operations | Explain-this-query (show resolved `QuerySpec` + scope sentence + generated SQL + footer reasons); saved-query library; notification triggers; project-status dashboards |
| Cross-app bridges | CCP4i2-unique — the compounds app lives next to crystallography and constructs | *"Structures I've solved for AR degrader complexes"* (Compounds ↔ Projects ↔ Moorhen); *"Constructs that produce the protein my AR degraders were tested against"*; *"Stock of my top-10 hit list"* (once `INVENTORY_PROPOSAL.md` ships) |
| Chemistry-aware queries | Substructure filter as a per-row predicate — covered in §19.3–19.4 | *"HTRF values of ARd compounds containing pyrimidine"*; *"compounds structurally similar to NCL-00026042 with IC50 < 10 nM"* |
| Chart-type extensions | Same resolver/executor, richer output formats (`scatter` / `histogram` / …) + multi-axis + highlight rules — wires into the existing aggregation-page scatter component. Covered in §19.5. | *"Plot cell-free HTRF against CDK4 nanobret data for compounds in the CDK4 project, highlighting pyrimidine-containing compounds in green"* |
| Generative analysis | Separate architecture; LLM sees data here, unlike v1 — not a v2 of this proposal | *"Summarise what's driving IC50 in the pyrimidine series"*; *"one-page triage report of top 20 ARd hits and their ADME liabilities"* |
| Interaction modes | Different surface, same backend | Voice input for bench-side use; NLP bar embedded in a lab notebook entry; weekly email digest from a saved query |
| Write operations | Locked out by decision 11 for v1 — listed here so the shape is acknowledged but the lock-out is not forgotten | *"Send my top-10 ARd compounds to the Pharmaron ADME queue"*; *"tag all IC50 < 10 nM compounds as hit series 3"* |

### 19.3 Chemistry-aware queries (v2's most natural first extension)

The most naturally-phrased chemistry query class — *"HTRF values of ARd compounds containing pyrimidine"* — fits the v1 architecture with **no change to the LLM's role**. The LLM still just emits a name (*"pyrimidine"*); the backend owns the authoritative chemistry catalog.

Implementation shape:

- **New resolver.** `resolve_scaffold(hint) -> ScaffoldResolution` follows the `Resolved / Clarify / Miss` triad established in §18.2. Matches against a curated `Scaffold` table (`name`, `aliases`, `smarts`, optional `svg`), edited by a curator the same way `Gene` rows are edited today. Clarify when a fragment name is ambiguous (*"pyrimidine"* vs *"pyrimidinone"* vs *"fused-pyrimidine core"*); miss with visual suggestions.
- **New executor predicate.** Substructure match is a per-row predicate in the same bucket as unit-threshold. `FILTER_SUBSTRUCTURE_ABSENT` slots next to `FILTER_THRESHOLD_NOT_MET`; `EXCLUDE_STRUCTURE_UNPARSEABLE` handles SMILES parse failures at the row level.
- **New `QuerySpec` field.** `scaffold_hint: Optional[str]` — the user's typed fragment name, verbatim.
- **No SMARTS in the prompt, ever.** The LLM never sees SMARTS and never emits SMARTS. Same discipline as "no catalogs in the prompt" for targets/protocols (§8, decision 2). This eliminates the class of failures where the LLM hallucinates structurally-plausible-but-wrong SMARTS.

**Scale path.** Naive per-row `rdkit.Chem.Mol.HasSubstructMatch` is fine at DDU-scale (low thousands of compounds) and should be the first slice. RDKit pattern fingerprints in a sidecar table add a screening layer for larger libraries (SSS via FP bit-subset match before exact substructure). The Postgres RDKit cartridge is the proper endgame for very large libraries but needs a prod/SQLite-tests split story. **Build naive first; add FP screening only when a real library hits the wall** — premature FP infrastructure is the canonical trap in this direction.

**Boundary with generative SAR prose.** *"Compounds containing pyrimidine with IC50 < 10 nM"* is tabular — fits this bucket. *"What's driving IC50 in the pyrimidine series?"* is generative SAR — belongs in bucket 5 (generative analysis) with a different LLM data-visibility model. Keeping the two in separate proposals matters.

### 19.4 Chemistry queries collapse into another resolver

A useful framing: **chemistry-aware queries are structurally another resolver**, not a new subsystem. The substructure resolver carries the same `Resolved / Clarify / Miss` triad, the same field-tagging for the orchestrator, the same LLM-emits-name-backend-owns-catalog discipline. What differs is the matching engine (RDKit `HasSubstructMatch` instead of string normalization or token-set containment) and the candidate payload shape (SVG drawing + SMARTS in place of gene symbols or `n_runs`). The executor doesn't grow new surface; it just handles one more optional `QuerySpec` field.

**The family this seeds.** Compound-by-ID resolver (*"compound NCL-00026042"*); similarity resolver (*"structurally similar to X"* — Tanimoto k-NN with a tuneable threshold); scaffold-by-project resolver (project-local series curated as `Scaffold` rows linked to a `Target`). All cast into the same mould.

**When to abstract.** Four-plus concrete resolvers following the same shape would justify lifting a generic `Resolution[T, C]` protocol. We're at three today (target / protocol / — and arguably metric + unit as a row-level predicate which doesn't quite fit the same triad). Substructure makes the fourth true-triad resolver; if compound-by-ID or similarity also land, the pressure to abstract becomes clear. Default position: **wait for the push, don't pre-empt the abstraction.** The cost of premature generality here is mostly design-lock-in around candidate payload shapes that are genuinely different per resolver, for the benefit of saving a few lines of type declarations.

### 19.5 NLP-driven chart generation (scatter + narrative selections)

The aggregation page already ships a scatter-plot component. A natural v2 extension is to let the NLP interface drive it — same query door, different output format — so prompts like *"Plot cell-free HTRF against CDK4 nanobret data for all compounds in the CDK4 project, highlighting pyrimidine-containing compounds in green"* render directly into that existing component. Visual reasoning over multi-assay data is how chemists actually triage, and the leverage-to-effort ratio here is high because **the viz surface already exists** — this direction is mostly backend payload shaping plus a thin palette layer, not a new frontend build.

The discipline to preserve: the LLM parses intent into a **structured `ChartSpec`**, the backend resolves and runs the query, the frontend renders using the typed component. **The LLM never generates plotting code** — no matplotlib snippets, no Vega-Lite JSON, no D3 code. That is the line between an audit-friendly, reproducible feature and a magical-but-silently-wrong one. Chart-code generation by LLM is where comparable systems fall over; keeping it out is a load-bearing rule.

**What's reused.** Per-axis resolution uses `resolve_targets` and `resolve_protocol` as-is — each axis is an independent (protocol, metric) pair filtered through the same scope. The substructure resolver from §19.3 provides the *"pyrimidine-containing"* predicate for highlights. Scope filtering is shared across axes (one registration target, one assay target, one pair of filters applied to the scoped queryset). The executor walks the scope once and emits one data-point per compound carrying both axis values.

**What's genuinely new.**

- **Output format field.** `output: "table" | "scatter" | ...` on `QuerySpec`. LLM infers from verbs — *"plot"* / *"chart"* / *"compare"* → scatter; *"tabulate"* / *"list"* / *"show"* → table; *"distribution of"* → histogram. A small, bounded vocabulary; each format has its own axis-cardinality rules (scatter = 2, histogram = 1, heatmap = 2 + categorical, …).
- **Multi-axis support.** `axes: list[Axis]`, each `Axis = (protocol_hint, metric?)`. Metric can be null — the backend infers the dominant KPI from the protocol's valid AnalysisResults if there is exactly one, otherwise clarifies. This **also quietly improves the table-output path**: *"cell-free HTRF data for ARd compounds"* no longer forces the user to name the KPI.
- **Highlight rules.** `highlights: list[HighlightRule]`, each a `(predicate, style)` pair. The **predicate reuses existing filter machinery** — substructure, metric threshold, compound-list, scope difference. Nothing new on the resolver side; a highlight predicate has the same shape as a scope predicate, applied at the point-labelling stage rather than the inclusion stage. The **style is a curated vocabulary** — ~8 named colours mapped to colorblind-friendly hex, a few marker shapes. Not arbitrary CSS. Palette resolver follows the §19.4 pattern: LLM emits a name, backend maps to a palette entry, clarifies when ambiguous (*"red"* vs *"red outline"* vs *"red-filled circle"*), misses with suggestions.

**Composability.** Because highlights are `(predicate, style)` pairs, they compose naturally. *"Pyrimidines in green and compounds above 10 µM IC50 in red outlines"* is two rules. First-match-wins in emitted order, so the user's phrasing controls priority. Compound-list highlights (*"highlight NCL-00026042 in red"*) unlock personalised markup once a compound-by-ID resolver lands — another instance of the §19.4 family.

**QuerySpec extension sketch.**

```python
@dataclass
class Axis:
    protocol_hint: str
    metric: Optional[str] = None  # None → infer from protocol's dominant KPI

@dataclass
class HighlightRule:
    style: str                                 # "green" / "red-circles" / etc. — resolved via palette catalog
    # Reuses existing filter machinery — exactly one of the following is set:
    scaffold_hint: Optional[str] = None
    threshold: Optional[Threshold] = None
    compound_ids: Optional[List[str]] = None

@dataclass
class QuerySpec:
    # ... existing fields unchanged ...
    output: str = "table"                      # "table" | "scatter" | "histogram" | ...
    axes: Optional[List[Axis]] = None          # required when output != "table"
    highlights: List[HighlightRule] = field(default_factory=list)
```

**Boundary — what this is not.** The LLM does not generate plotting code, does not choose visual encodings freely, does not invent palette entries. The backend does not compute summary statistics beyond what a scatter plot naturally reveals. Prompts that drift into *"...and tell me what the structure-activity trend is"* cross into generative SAR territory (§19.2 bucket 5) — different architecture, separate proposal, different data-visibility model for the LLM. The scatter-plot door is strictly tabular-query-plus-rendering.

**Integration cost.** Likely small. The existing aggregation-page scatter already accepts a structured payload; wiring the NLP executor to emit that payload is additive. The palette catalog + style resolver is a new thin module. The multi-axis + highlight fields on `QuerySpec` are purely additive — no changes to v1 fields, no changes to the v1 executor-for-tables path. The v1 invariants (resolver pattern, LLM-as-parser, no fabrication) hold unchanged.

### 19.6 Filler nouns and prompt-engineering discipline

Chemistry prompts liberally use conversational nouns — *"hits"*, *"compounds"*, *"molecules"*, *"ligands"*, *"series"*, *"analogues"* — that **carry no structural meaning and should be stripped by the LLM parser without generating spurious filter fields**. *"ARd hits containing pyrimidine"* and *"ARd compounds containing pyrimidine"* must resolve to the same `QuerySpec`.

The discipline: the system prompt teaches the LLM with examples that these words are filler in the absence of explicit cutoff language. A prompt that **does** name a cutoff (*"hits defined as IC50 < 10 nM"*, *"the top-10% series"*) carries structure and should parse accordingly. The LLM is not asked to infer a hit-cutoff from the word *"hits"* alone — doing so would invent filter criteria the user did not specify.

An earlier draft of §19.3 mistakenly conflated *"hits"* with per-protocol `Protocol.target_value` cutoffs. That would be a schema-level inference from conversational filler — wrong, and it would produce result sets the user didn't ask for. The `target_value`/`poor_value` fields are legitimate once the user explicitly invokes them (e.g. *"compounds passing the protocol's own hit cutoff"*), but not implicitly from filler-noun phrasing.

This applies generally: **the LLM must not synthesise filter fields from words that only exist to make the sentence flow**. Good v1 discipline even without v2 extensions; worth baking into the system prompt's few-shot examples when slice 5 lands.

### 19.7 First picks after v1

Two directions have the clearest leverage-to-effort ratio and are worth naming as top candidates for the post-v1 slice queue:

- **Compound-rooted queries** — unifying *"ARd compounds with no HTRF measurement yet"* (missing-data), *"compounds registered in 2025"* (registration date range), and *"compounds registered by Alice"* (provenance). Highest-value single direction — chemistry decision-making is forward-looking, these phrasings are the first thing a user tries on the command bar, and they all share one architectural shift: root the executor's queryset at `Compound` rather than `AnalysisResult`, and relax the v1 hard requirement for `metric` + `protocol_hint`. Concretely this is (a) new `QuerySpec` fields (`registration_date_range`, `registered_by_as_typed`, `missing_measurement_for_protocol_hint`), (b) a compound-only execute path that applies those predicates, and (c) a handful of new system-prompt examples. No new resolver — target resolution is already optional-per-field. Ship the date-range + missing-measurement cases together because they share the same backend lift. **Flag for when this lands**: today's executor returns `SpecError` when `metric`/`protocol_hint` are null — the relaxation must be gated on "at least one compound-level predicate is present" to keep the "empty spec" failure mode. Without that guard the command bar silently dumps the full Compound table.
- **Explain-this-query.** One-click reveal of the resolved `QuerySpec`, scope sentence, generated SQL, and footer reasons. Every chemist who uses an LLM over their data asks *"what did it actually do?"* Shipping transparency alongside the feature (not as a retrofit) is the single highest trust-multiplier. The scope-sentence generation in slice 4 already does 80% of this; exposing it is mostly a UI lift.

Two near-tied thirds, each with a different dependency profile — sequence based on which dependency is cheaper:

- **NLP-driven scatter** (§19.5). Gated on the existing aggregation-page scatter component being easy to drive from a structured payload. If that integration is cheap (likely), this may well be the single most compelling demo — multi-axis data visualisation with narrative highlights — and is probably the cheapest high-impact win.
- **Chemistry-aware substructure** (§19.3). Gated on curator time to seed the initial `Scaffold` catalog. The curator workflow is a real cost but one-time; once seeded, every chemistry-aware prompt benefits.

### 19.8 Deliberate non-goals for v2

To protect against scope creep in the v1 → v2 transition:

- **Generative SAR narrative** needs a separate proposal — different data-visibility model for the LLM, different prompt architecture, different evaluation harness. Do not graft it onto this one.
- **Chart code generation by the LLM** (§19.5 boundary) — the LLM emits a `ChartSpec`, never plotting code. This is an architectural non-negotiable for the scatter extension.
- **Write-capable commands** remain locked out by decision 11 until a trust/audit model is specified in its own proposal.
- **Cross-session conversational memory** is a distinct problem from extending the `QuerySpec` schema and should not be conflated with "more query shapes". Revisit when v1 usage data justifies a chat-style surface.
