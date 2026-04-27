# Natural-Language Query Proposal (Compounds App)

**Status**: v2 **pivoted** 2026-04-23 — the NLP endpoint now produces a *compound selection* and redirects to `/assays/aggregate` rather than rendering its own table. Pivot driven by first-day user feedback: the aggregation page is the richer display surface, and the cross-protocol selectivity failure mode (mEGFR WT/TM) collapsed for free when filters became a list. See §18.2 slice-9 for the pivot diff.
**Author**: Martin Noble (with Claude)
**Date**: 2026-04-22 (drafted); §§18–19 added 2026-04-23; §18 slices 4–8 updates 2026-04-23; §18 slice-9 pivot 2026-04-23

## 0. Pivot note (2026-04-23)

After v1 shipped, direct user feedback on DDU surfaced that:

1. Chemists immediately try **cross-protocol selectivity** queries (filter by WT IC50, view TM IC50). The v1 single-`protocol_hint`/single-`metric` `QuerySpec` couldn't express this — the LLM collapsed both halves, producing indistinguishable results for semantically different queries.
2. The aggregation page already owns the *visual* display — project cards with spider plots auto-rendered from configured data relationships, plus compact / medium / long / pivot formats, plus phys-chem column toggles. The NLP-rendered table (v1) duplicated a thin slice of that machinery.

v2 pivots the subsystem to produce a **compound selection**, not a table:

- `QuerySpec` (single protocol + metric + threshold + columns) → `CompoundSelector` (target + list of AND'd `MeasurementFilter`s).
- Executor output is a list of compound IDs plus the target/protocol names used in the filter, not a row/column payload.
- View returns a `redirect_url` pointing at `/assays/aggregate?targets=…&compound=…&protocols=…&format=cards`. The user lands on project cards with the filtered compound set.
- Cross-protocol selectivity falls out from the list-of-filters shape: *"mEGFR compounds with WT IC50 < 10 uM AND TM IC50 > 1 uM"* → two filters intersected, correct compound subset, protocols pre-populated in the aggregation page's URL.

Sections §5 (schema), §7 (executor), §10 (UI) describe the *v1* shape. Read them in context of the pivot: the structural intent (LLM-as-parser, deterministic backend, clarify-as-UI, filler-noun discipline, managed-identity auth) is unchanged; the schema and output shape are the specific things that changed. §18.2 slice-9 documents the concrete diff; §18.3 slice-9 the decisions made during the pivot.

§15's "cross-assay aggregations" item is now **shipped** — it was the headline of the pivot. §19.7 updated accordingly.
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

### 16.1 Step 0 — Target hydration prerequisite (done on DDU, pending on Demo/Kawamura)

The Target-model dependency and hydration pipeline are shipped. `Gene` model + `Target.genes` M2M + `hydrate_genes` command + serializer read/write of genes + `backfill_target_genes.py` are deployed to all three instances. DDU Database is fully backfilled (32 of 44 targets linked to gene_symbols; 1 deferred (GLK — needs disambiguation); 11 marked `skip` as non-gene programmes). Demo and Kawamura have the code but haven't yet been backfilled.

Before turning NLP on for Demo/Kawamura: run `backfill_target_genes.py --dump` on those instances and curate the plans, so the matching pool is as rich there as on DDU.

### 16.2 Azure infrastructure prerequisites for the LLM path

Survey of the live resource group (2026-04-23) showed the infra was further along than the original draft of this section assumed:

- A shared Azure OpenAI resource `ddu-openai` already exists in `ccp4i2-bicep-rg-uksouth` — provisioned for the aacr-abstracts container app (which authenticates via API key against the regional endpoint).
- A `gpt-4o` model deployment (GlobalStandard, capacity 10, version 2024-11-20) is already live on that resource.
- The resource does **not** have `customSubDomainName` set, which is a hard prerequisite for managed-identity auth (MI token validation uses `https://<name>.openai.azure.com/`, not the regional endpoint).

So the §16.2 work reduced to three concrete code/ops changes and a small one-time `az cli` step:

| # | Task | Status | Notes |
|---|---|---|---|
| 1 | Set `customSubDomainName = "ddu-openai"` on the shared OpenAI resource | Script shipped: `Docker/azure-uksouth/scripts/configure-openai-for-nlp.sh` | Idempotent az-cli one-shot. One-way Azure change (cannot be removed) but additive — aacr-abstracts' existing regional endpoint + API key keep working. |
| 2 | Role assignment `Cognitive Services OpenAI User` (GUID `5e0bd9bd-7b93-4f28-af87-19fc36ad61bd`) on the OpenAI resource, scoped to `containerAppsIdentity` | Script applies it; also mirrored into `infrastructure/infrastructure.bicep` so a fresh full deploy picks it up without running the script | Same GUID, same identity — matches the Storage/KeyVault pattern in the same Bicep module. |
| 3 | Server Container App env vars: `AZURE_OPENAI_ENDPOINT`, `AZURE_OPENAI_MODEL` (default `gpt-4o`), `COMPOUNDS_NLP_ENABLED` (default off) | `applications.bicep` params + env block + `deploy-applications.sh` wiring; documented in `.env.deployment.template` | Per-instance values live in `.env.deployment` / `.env.demo` / `.env.kawamura`. Endpoint string is `https://ddu-openai.openai.azure.com/` once step 1 has run. |

**Runtime dependency** — `openai>=1.50,<2.0` was added to `Docker/server/Dockerfile` during slice 5. Picked up on the next server image build and push; nothing to install on the Container App itself.

**Enablement sequence** (per instance):

1. (Once, shared across instances) — run `Docker/azure-uksouth/scripts/configure-openai-for-nlp.sh`. Idempotent; safe to re-run. Prints the managed-identity endpoint URL to copy into env files.
2. In the instance's `.env.*` file, set:
   - `AZURE_OPENAI_ENDPOINT=https://ddu-openai.openai.azure.com/`
   - `COMPOUNDS_NLP_ENABLED=true`
3. Rebuild + push server image if slice 5 hasn't shipped to prod yet (so `openai` is in the container).
4. Redeploy: `./scripts/deploy-applications.sh [--env .env.<instance>]`.

**Non-blocking for earlier slices.** None of this gates slices 6 (view) or 8 (eval). Slice 6 ships behind the flag off-by-default; eval runs locally against mocks. §16.2 only gates actually *using* the feature in a deployed instance.

**Coexistence with aacr-abstracts.** aacr-abstracts authenticates to the same `ddu-openai` resource via API key against the regional endpoint. Setting `customSubDomainName` is additive — the regional endpoint continues to resolve, and `disableLocalAuth` is deliberately left unset, so API-key auth keeps working. The two apps share the resource's token-per-minute quota; if they ever contend, the `§11` "per-instance separation" escape hatch is straightforward (deploy a second OpenAI resource; change one env var).

### 16.3 App rollout

- Ship NLP v1 behind a feature flag / env var (`COMPOUNDS_NLP_ENABLED`) so it can be enabled per-instance independently of the code deploy.
- Demo instance first, with eval golden set passing.
- Kawamura instance second, with instance-specific targets/protocols exercised in the golden set.
- No database migration needed for NLP v1 (purely additive, read-only).

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
| 6 | DRF view with clarify loop | Shipped | §7 (view), §9, §16.3 feature flag |
| 7 | Frontend command-bar (landing page for v1; per-project deferred) | Shipped | §10 |
| 8 | Evaluation harness + golden set | Shipped | §12 |
| 9 | **Pivot** — QuerySpec → CompoundSelector; output → redirect to `/assays/aggregate` | Shipped | §0 pivot note, §19.7 update |
| 10 | Date filters — `registered_date_range` + per-filter `assay_date_range` | Shipped | §7, §19.7 item 1 partial |
| 11 | User filters — `registered_by_as_typed` + per-filter `assayed_by_as_typed` with User resolver + clarify | Shipped | §6 (user resolution), §19.7 item 1 closing |
| 12 | **Assay-selection queries** — second query family (AssaySelector) alongside CompoundSelector; LLM dispatches via `assay_selector` nested object; lands on `/assays?ids=…` | Shipped | §5 schema, §7 dual-executor, §19.7 update |
| 13 | **Substructure / scaffold queries** — curated name→SMARTS catalog + RDKit on-the-fly `HasSubstructMatch`; `scaffold_hints` list on both selectors; LLM emits typed names verbatim, never SMARTS | Shipped | §19.3, §19.4, §19.7 third-position update |
| 14 | **Compound-ID pinning** — `compound_refs_as_typed` list on `CompoundSelector`; UNION semantics (additive, not narrowing); deterministic resolver via `compounds.formatting.extract_reg_number`; "compound 26007 and recent compounds" lands a known lead alongside a filtered set for comparison | Shipped | §19.7 update |
| 15 | **Diagnostic empty selections — MetricMiss** — lenient KPI match (case + punctuation insensitive), correct `metric=None` semantics, and a `MetricMiss` with `available_metrics` when the user's typed metric isn't recorded in any in-scope row's KPI. *"ARd HTRF IC50 < 50"* with rows recording `pIC50` no longer silently returns "Found 0" — the response answers itself with the available KPIs. | Shipped | §6.3 (revisited), §18.3 slice-3 footnote |

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

**Slice 6 — DRF view** (182 tests cumulative)

- New `apps/compounds/nlp/view.py` — a single `@api_view(['POST'])` function `nlp_query` mounted at `/api/compounds/nlp/query/` via `apps/compounds/urls.py`.
- Two request shapes:
  - Fresh query: `{"prompt": "<english>"}` — runs `parse_prompt` → `execute`.
  - Continuation (clarify picker): `{"spec": {<QuerySpec dict>}}` — skips the LLM entirely (decision 5), parses the round-tripped spec, runs `execute`. The UI fills `registration_target_id` / `assay_target_id` / `protocol_id` on the spec to pin the user's picker choice.
- Response shape: a single-object union discriminated by `status`: `"table" | "clarify" | "miss" | "not_a_query" | "error"`. HTTP codes: 200 for table/clarify/miss/not_a_query, 400 for scope/spec/bad_request errors, 429 for rate limit, 502 for parse errors, 404 for feature-off.
- Auth: `[IsAuthenticated]` — matches the compounds-app convention. DRF's `AzureADAuthentication` handles prod tokens, `DevAuthMiddleware` handles local dev.
- Feature flag `COMPOUNDS_NLP_ENABLED` env var gates the endpoint (returns 404 when off). See §16.3.
- Daily-call cap of 500/user (§11, decision 17) via Django's default cache. Counter keyed `nlp_cap:<user_pk>:<iso_date>`, TTL ~26h. Per-process `LocMemCache` is acceptable because the cap is an anti-runaway guard, not a cost gate — a user hitting N workers could reach N×500 in the worst case, which is still orders of magnitude above real human workload.
- Pinning-ID support added to `resolve_target` / `resolve_protocol` / `execute` as optional `pinned_id` params threaded through from the new `QuerySpec.{registration,assay}_target_id` / `protocol_id` fields. When set, the resolver skips fuzzy matching and fetches the entity directly — with a scope-validity check for protocol (a stale UUID from another target's clarify response can't slip through).
- Tests: `apps/compounds/nlp/tests/test_view.py` (17 tests). Exercises the feature-flag gate (off, explicit-false), body validation (missing, both-prompt-and-spec, non-string prompt, malformed threshold in spec), happy path through parse_prompt (mocked) → execute → table, continuation skipping the LLM, pinning IDs bypassing resolution, clarify/miss/scope-error/spec-error/not_a_query/parse-error HTTP mapping, and the daily cap (per-user isolation + cap-reached 429).

**Slice 8 — Evaluation (golden set + offline + online eval)** (193 tests cumulative; 2 online tests skipped by default)

- New `apps/compounds/nlp/tests/golden/prompts.yaml` — 35 curated (prompt, expected) pairs covering:
  - single-target with/without threshold (various metrics: IC50, EC50, Ki, Kd, CLint, pIC50, % inhibition; units: nM, uM, pM, %, uL/min/mg; operators: <, <=, >, >=, =);
  - cross-project and assay-only prompts;
  - column preset variations (phys_chem default, lipinski, columns_explicit);
  - filler-noun discipline (§19.6) — "hits" / "series" / "analogues" must NOT synthesise a threshold unless the user names a cutoff;
  - ambiguous-target cases where the LLM must emit the user's string verbatim and let the resolver Clarify downstream;
  - not_a_query rejections (SAR summary requests, register/delete commands, meta-questions about assay mechanism).
- New `apps/compounds/nlp/tests/test_golden.py` — two test groups in one module:
  - **Offline** (11 tests, always runs, ~40ms total): structural validation of the YAML (required keys, unique names, valid types, valid threshold ops, valid column presets, every spec entry names a metric + at least one target, columns_explicit shape) plus a round-trip check that converts each entry's expected shape into an LLM-output dict and feeds it through `_to_parse_result`, asserting it parses back to the same QuerySpec or NotAQuery. Protects against YAML typos and drift as the set grows.
  - **Online** (2 tests, skipped unless `RUN_NLP_ONLINE_EVAL=true`): actually calls `parse_prompt` against the live Azure endpoint per golden entry. Measures:
    - **spec match rate** across `type: spec` / `type: spec_with_ambiguity` entries — asserts ≥90% per §12.
    - **not_a_query match rate** — asserts 100% per §12 (any "parses as spec where it shouldn't" is a false positive).
    - **parse-error rate** in a separate test — asserts 0. ParseError means our glue got junk from upstream, distinct from the LLM emitting a spec that disagrees with the golden file. Catching these separately makes failure modes legible.
- Diff helper `_spec_diff(actual, expected)` on mismatch produces a `"; "`-joined human-readable string per field — makes online-eval failures actionable at a glance.

**How to run the online eval:**

```bash
source ../../ccp4-20251105/bin/ccp4.setup-sh
cd server
AZURE_OPENAI_ENDPOINT=https://ddu-openai.openai.azure.com/ \
RUN_NLP_ONLINE_EVAL=true \
PYTHONPATH="$PWD:$PWD/../apps" DJANGO_SETTINGS_MODULE=compounds.settings \
  ccp4-python -m pytest ../apps/compounds/nlp/tests/test_golden.py -v -s
```

Requires `openai` installed in the ccp4-python env (`ccp4-python -m pip install 'openai>=1.50,<2.0'`) and either managed-identity auth via `az login` (developer's AAD account must have the Cognitive Services OpenAI User role on `ddu-openai`) OR `AZURE_OPENAI_API_KEY=<key>` set as the dev-only fallback.

**Also in this slice** — Dockerfile fix: `apps/compounds/nlp` is now in the `COPY` list for the server image (`Docker/server/Dockerfile`), so the NLP backend actually ships in production. Missed in slices 5/6 because local tests passed off the source tree; caught now before the first deployed enablement.

**Slice 9 — Pivot to CompoundSelector + redirect to aggregation page** (192 tests cumulative; net -2 from shape changes, +selectivity cases)

The headline change: QuerySpec (single protocol/metric/threshold) → CompoundSelector (target + list of AND'd MeasurementFilters). The executor no longer produces a TablePayload; it produces a CompoundSelection (list of compound IDs) which the view turns into a redirect URL for `/assays/aggregate`.

- `spec.py` rewrite:
  - **Removed**: `QuerySpec`, `Threshold` kept, `TableRow`, `TablePayload`, `COL_PRESET_PHYS_CHEM`, `COL_PRESET_LIPINSKI`, `PHYS_CHEM_COLUMNS`, `LIPINSKI_COLUMNS`.
  - **Added**: `MeasurementFilter` (protocol_hint / metric / threshold / protocol_id), `CompoundSelector` (registration/assay target fields + filters list + pinning ids), `CompoundSelection` (compound_formatted_ids / target_names / protocol_names / n_matched / scope_sentence).
  - `PromptParseResult = Union[CompoundSelector, NotAQuery, ParseError]` (was QuerySpec).
  - `ExecutionResult` swaps TablePayload for CompoundSelection.
  - `ProtocolClarify` / `ProtocolMiss` gain a `filter_index` field so the UI knows which filter to pin on a protocol-level disambiguation.
- `resolver.py`: unchanged target-resolution logic. `resolve_targets` takes a `CompoundSelector`. `resolve_protocol` gains `filter_index` to tag Clarify/Miss results.
- `executor.py`: rewritten as a compound-intersector. For each filter, computes the set of compound IDs whose AnalysisResult rows pass (valid + metric match + threshold — evaluator logic unchanged). Intersects across filters. No filters = all compounds in scope. Hydrates compound IDs into formatted identifiers (`NCL-00031070` etc.) ordered by `reg_number` for stable URLs.
- `view.py`: response body carries `redirect_url` constructed from the selection. URL uses the `/assays/aggregate` contract (`targets=`, `compound=`, `protocols=`, `format=cards` default). `partial_selector` on clarify responses.
- `llm.py`: new system prompt emphasising *"describe WHICH compounds, not WHAT columns"*; new schema with `measurement_filters` as a list of `{protocol_hint, metric, threshold}` objects. Strict-mode schema guards updated.
- Golden set rewritten: 35 entries including new cross-protocol-selectivity cases (mEGFR WT/TM, CDK4/6), concordance (HTRF + AlphaLISA on same metric), no-threshold filter forms, display-preference-ignoring entries ("phys chem for X", "show me spider plots of Y" — LLM must parse the selection part and ignore the display hint).
- Frontend:
  - `nlp-api.ts` types mirror the new spec exactly. `postNlpQuery` unchanged. `applyClarifyPick` gains `filterIndex` for per-filter protocol pinning.
  - `NLPResults.tsx`: `TableView` replaced by `SelectionView` — shows `scope_sentence`, a large `Found N matching compounds` headline, a preview of the first 5 compound IDs (+ overflow count), and a primary `View as cards →` button that navigates to `redirect_url` via `window.location.assign`. Click-to-confirm UX per user decision (cards rendered server-side, not NLP page).
  - `NLPPanel.tsx`: unchanged orchestration, but the `handleClarifyPick` now threads `filterIndex` through from protocol-clarify responses.

- Tests: `test_executor.py` reworked with a `selectivity_world` fixture that demonstrates the intersection semantics (c1 passes both filters, c2 fails TM, c3 fails WT → correct single-compound selection). `test_view.py` asserts redirect-URL shape + protocol-names list + partial_selector round-trip. `test_llm.py` hits the new schema + new `_to_parse_result` selector branches. `test_golden.py` rewritten for the selector YAML shape.

Suite: 192 NLP tests passing (was 194 pre-pivot — the net delta reflects adjusted test structure, not capability loss). Frontend TypeScript compiles cleanly.

**Slice 7 — Frontend command-bar** (194 tests cumulative; one continuation-round-trip view test added for the slice-6 `partial_spec` retrofit)

- New page at `apps/compounds/frontend/app/nlp/page.tsx` — `/nlp` route. Discoverable via an "Ask (Natural-Language Query)" tile on both landing pages (`app/page.tsx` standalone + `app/app-selector/page.tsx` Docker overlay). Per-project placement deferred until per-project compounds pages exist — §18.3 slice-6 decision.
- New `apps/compounds/frontend/lib/compounds/nlp-api.ts` — TypeScript types mirroring `apps/compounds/nlp/spec.py` (QuerySpec, Threshold, TargetCandidate, ProtocolCandidate, TableRow, the NLPResponse discriminated union) + a thin `postNlpQuery` fetcher that wraps `authFetch` but (unlike the existing compounds `apiPost`) returns the JSON body on *any* HTTP status. Necessary because 4xx/5xx carry structured data the UI wants to render — 404 for disabled, 400 for spec/scope error, 429 for rate limit, 502 for parse error. Throws only on network / non-JSON failures. A convenience `applyClarifyPick(partial_spec, field, pickedId)` helper pins the chosen id on the correct `*_target_id` / `protocol_id` field for the continuation re-POST.
- New component `apps/compounds/frontend/components/compounds/nlp/NLPPanel.tsx` — owns the React state (prompt, response, loading, error). Multi-line `TextField` with Enter-to-submit (Shift+Enter for newline), prominent "Ask" button with spinner. Below the input: four example-prompt `Button`s that one-click into the bar and submit — the command-bar equivalent of the "not_a_query" examples panel (§19.6). Hidden once a response is present.
- New component `apps/compounds/frontend/components/compounds/nlp/NLPResults.tsx` — switches on `response.status`:
  - `table` → scope sentence + `MUI Table` (sticky header) with compound-id (linked to the registry compound page), value+unit, and one column per phys-chem field from `property_columns` using `PROPERTY_LABELS` from the aggregation-table shared module (keeps label convention consistent across the app). Footer `Paper` with `Chip`s per `footer_excluded` reason using human-readable labels (`unit_unknown` → *"no unit recorded"*, etc.).
  - `clarify` → non-customized MUI `Chip`s (per user decision), one per candidate, each showing `name` + gene symbols (targets) or `n_runs / n_compounds / last_run` (protocols). Click → `onClarifyPick(partial_spec, field, id)` from the panel, which re-POSTs with the pinning.
  - `miss` → `Alert severity="info"` with suggestion `Chip`s.
  - `not_a_query` → `Alert` with the LLM's reason + a static "what this bar can do" hint.
  - `error` → `Alert severity="warning"` for `kind="disabled"`, `severity="error"` otherwise, with the `kind` and optional `field` shown as sub-text for diagnostics.
- `lib/compounds/routes.ts` gains `routes.nlp.home()` returning `/nlp`.

**Slice-6 fix retrofit (bundled in slice 7):** `_serialize` now includes `partial_spec` on `clarify` responses — the QuerySpec that produced them, echoed back so continuations skip the LLM (decision 5). The slice-6 original missed this; it surfaced now because the frontend needs it to build a pinned continuation from a picker choice. New view test `test_clarify_response_round_trips_into_continuation` exercises the full pick → pin → continuation path and asserts the LLM was not re-called.

**Slice 13 — Substructure / scaffold queries** (240 tests cumulative; +21 for the slice)

The headline capability added by the slice: prompts like *"ARd compounds containing pyrimidine"*, *"CDK4 pyridines with HTRF IC50 < 100 nM"*, and *"CDK4 HTRF assays on pyrimidines"* now filter compounds / assays by on-the-fly RDKit substructure match. Architecture follows §19.3 exactly — the LLM never sees or emits SMARTS; it emits the user-typed name (*"pyrimidine"*, *"amide"*, *"trifluoromethyl"*) and the backend owns the catalog.

- New module `apps/compounds/nlp/substructures.py` — curated name→SMARTS catalog. ~30 entries across aromatic N-heterocycles, saturated N-heterocycles, functional groups, and carbocycles. Each `Scaffold` has `name` (canonical, also used as pinning ID), `aliases` (plurals, -yl suffix, common abbrevs), and `smarts`. Refine-and-redeploy: no migration, no retag; next query picks up the new semantics.
- `spec.py` additions: `FIELD_SCAFFOLD_HINT`; `ScaffoldCandidate`, `ResolvedScaffold`, `ScaffoldClarify`, `ScaffoldMiss`, `ScaffoldResolution` union. Both `CompoundSelector` and `AssaySelector` gain parallel `scaffold_hints: List[str]` + `scaffold_ids: List[str]` fields (pinning id for scaffolds is the canonical name — no UUID assignment). `ExecutionResult` grows to include Scaffold{Clarify,Miss}.
- `resolver.py` adds `resolve_scaffold(as_typed, *, scaffold_index, pinned_id) -> ScaffoldResolution`. Three-tier match: exact normalised → substring → fuzzy miss with `SequenceMatcher`-ranked suggestions. Matches the Target/Protocol/User resolver shape one-for-one, bringing the count of true-triad resolvers to four — which was the §19.4 trigger to consider a generic `Resolution[T, C]` protocol, not taken: candidate payload shapes remain genuinely different and the abstraction would save a few lines of types at the cost of design-lock-in.
- `executor.py` adds `_match_compounds_by_scaffold(compound_ids, resolved_scaffold)` — lazy RDKit import so the module only pulls in RDKit when a prompt actually names a scaffold. Scaffold narrowing is applied **before** measurement filters in both the compound path and the assay path (cheap compound-level predicate; intersect first). Multiple scaffold hints are ANDed. Scope sentence grows a *"containing pyridine AND amide"* phrase.
- `llm.py` adds `scaffold_hints` to the top-level schema and the `assay_selector` subschema. System prompt documents the substructure handling (LLM emits verbatim name; never SMARTS; single-multiple scaffolds each as its own list entry; plural/alias forms pass through and the backend resolver's catalog handles normalisation). Strict-mode JSON schema guards updated on both subschemas.
- `view.py`: `_selector_from_dict` / `_assay_selector_from_dict` deserialise `scaffold_hints` + `scaffold_ids`; `_serialize` includes ScaffoldClarify / ScaffoldMiss in the discriminated union dispatch.
- Frontend: `nlp-api.ts` adds `FIELD_SCAFFOLD_HINT`, `ScaffoldCandidate`, the new selector fields, `scaffold_index` on Clarify / Miss responses, `scaffoldIndex` on `applyClarifyPick` with a parallel-array `applyScaffoldPin` helper (hints + ids padded to index, write, preserve rest). `NLPResults.tsx` gains a scaffold branch in `CandidateLabel` (renders `name` + monospace SMARTS preview) and the clarify/miss label handling. `NLPPanel.tsx` threads `scaffoldIndex` through `handleClarifyPick`.
- Tests: new `tests/test_resolver_scaffold.py` (12 cases — three tiers, index/field tagging, pinning bypass). `test_executor.py` gains a `scaffold_world` fixture with pyridine / pyrimidine / amide / aliphatic compounds and 9 tests covering compound-path filtering, clarify / miss pass-through, pinning, two-scaffold AND (no compound hits both → empty selection), scope-sentence prose, composition with measurement filters, and the assay-selection path. `test_golden.py` grows a scaffold-aware diff and a scaffold-hints round-trip assertion on the assay path. `test_llm.py` schema-drift guard updated for the new `scaffold_hints` field on assay_selector.
- Golden set: 6 new entries covering single-scaffold, plural-verbatim, scaffold+threshold, two-scaffold-AND, functional-group (`trifluoromethyl`), and scaffold-on-assay-query.

The slice is deliberately RDKit-lazy at runtime: no scaffold mention in the prompt → RDKit stays un-imported and un-invoked. At DDU scale (few-thousand compounds per target), `HasSubstructMatch` is sub-second; the catalog's docstring notes pattern-fingerprint pre-screening as the escape valve if a real library hits the wall.

**Slice 14 — Compound-ID pinning** (262 tests cumulative; +19 for the slice)

The headline capability: prompts like *"compound 26007 and compounds made since January"*, *"NCL-00026007 alongside ARd pyrimidines"*, or *"just compound 12345"* let the user pin specific compound IDs into the selection alongside whatever the other predicates pick. Use case is "lead-compound-as-benchmark": a chemist scrolling through their recent set wants the known lead to stay visible for direct visual comparison. This is the **first additive operation** in the selector — every other predicate narrows; pins UNION on.

- Additions to `spec.py`:
  - `compound_refs_as_typed: List[str]` field on `CompoundSelector` (no parallel `_ids` field — there's no clarify path, so no pinning round-trip).
  - `FIELD_COMPOUND_REF` constant.
  - `ResolvedCompound`, `CompoundMiss`, `CompoundResolution` union — no clarify dataclass because compound IDs are deterministic.
  - `CompoundMiss` joins `ExecutionResult`.
- `resolver.py` adds `resolve_compound_ref(as_typed, *, ref_index) -> CompoundResolution`. Composes with the existing `compounds.formatting.extract_reg_number` helper, which already handles every prefix shape (`NCL-00026007` / `NCL26007` / `NCL 26007` / `26007` / `NCL000-26007`). One small function, no fuzzy tier — either the typed reference parses to an existing compound or it doesn't.
- `executor.py`:
  - New `_selector_has_non_pin_narrowing(selector)` helper distinguishes "selector has scaffold/user/date/filter narrowing" from "only pins". The pin-only case short-circuits the unscoped-base to the empty set so the final selection is exactly the pinned compounds (rather than every compound in the registry, which would be the naive "unscoped + UNION" result).
  - `_selector_has_no_narrowing_predicate` now treats pins as legitimate narrowing — a pin-only selector skips the friendly empty-selector ScopeError.
  - Pin resolution happens before measurement filters but its result is UNIONed in AFTER all narrowing — a pinned compound appears in the final selection even when it doesn't pass the threshold or sit in the target scope. This is the load-bearing semantics: a chemist pinning their lead expects to see it whether or not it survives the cutoff.
  - First unresolvable pin short-circuits to a `CompoundMiss` tagged with the offending ref's `ref_index` so the user sees which entry went wrong.
  - Scope sentence grows a *", plus NCL-00026007 (pinned)"* (or comma-joined list) clause.
- `llm.py` adds `compound_refs_as_typed: array of strings` to the top-level schema (required, strict-mode). System prompt teaches the pattern with examples spanning full prefix, bare integer, multi-pin, and mixed-with-other-predicates phrasings. Discipline is preserved: LLM emits the typed reference verbatim, never normalises.
- `view.py` deserialises `compound_refs_as_typed` from continuation selectors and adds `CompoundMiss` to the miss-branch dispatch. The redirect URL's `compound=` param naturally includes pinned IDs because the executor already UNIONed them into `CompoundSelection.compound_formatted_ids`.
- Frontend (`nlp-api.ts` + `NLPResults.tsx`):
  - New `FIELD_COMPOUND_REF` constant; `CompoundSelector.compound_refs_as_typed?: string[]` typing.
  - `CandidateKind` union grows a `'compound'` member; `MissView` renders compound-ref misses with a special helper line *"Check the format (e.g. NCL-00026007 or 26007) and that the compound has been registered."* — there are no fuzzy suggestions for compound IDs, so the message tells the user directly.
- Tests:
  - New `tests/test_resolver_compound.py` (8 tests — full ID, compact, lowercase, bare integer, unparseable, nonexistent, empty, ref_index tagging).
  - `test_executor.py` `pin_world` fixture spans two unrelated programmes; 9 new tests cover pin-only, pin+target, pin survives strict filter, bare-number, unknown ID, unparseable, mixed-with-bogus short-circuit, multi-pin selection, and scope-sentence rendering.
  - `test_view.py`: pinned compound appears in redirect URL + scope sentence; unparseable pin returns `compound_ref` miss.
  - `test_golden.py`: round-trip of `compound_refs_as_typed`, diff inclusion, the "target-less must have a narrowing predicate" guard now counts pins.
- Golden set: 5 new entries (pin-only full ID, pin-only bare number, pin + date range, pin + target + scaffold, two pins).

**Slice 15 — Diagnostic empty selections (MetricMiss + lenient KPI match)** (271 tests cumulative; +9 for the slice)

Problem: a chemist typing *"compounds registered to ARd where HTRF IC50 < 50"* against a protocol whose rows record `pIC50` (or `EC50`, or `IC-50` with a hyphen) saw `Found 0` — silently. The strict `results.get("KPI") != metric` check at the evaluator filtered every row with `FILTER_KPI_MISMATCH`; `FILTER_*` reasons aren't surfaced post-pivot, so the response carried no diagnostic at all. The slice fixes both halves:

- **Lenient KPI match** in [evaluator.py](apps/compounds/nlp/evaluator.py). User's typed metric and the row's stored `KPI` are both passed through `compounds.utils.normalize_ref` (lowercase + strip non-alphanumerics) before comparison, so *"IC50"* / *"ic50"* / *"IC-50"* / *"ic_50"* all unify. Once a match is established, the value-lookup key uses the row's stored KPI string verbatim (the row's own column name is the source of truth for the numeric column).
- **`metric=None` semantics fixed.** Pre-fix, `metric=None` triggered `if results.get("KPI") != None`, which filtered every row that had a KPI string — directly contradicting the system-prompt contract ("metric null → any numeric KPI counts"). Now `metric=None` means "any KPI counts; look up the value via the row's own KPI key".
- **MetricMiss surfacing** in [executor.py](apps/compounds/nlp/executor.py). Before evaluating each filter, a small pre-check walks the filter's scope queryset (target + protocol + date + assayed-by narrowed) and projects just the JSON `KPI` field; if no row's KPI matches the user's typed metric (lenient), the executor returns a `MetricMiss(query, available_metrics, filter_index, field)` with the distinct KPIs that ARE in scope. The frontend renders these as chips so the response answers itself ("did you mean `pIC50`?").
- **Shared scope helper.** `_filter_scope_qs` is factored out so the metric pre-check and the row-iterating filter loop see identical scope semantics — no risk of one path narrowing to a slightly different set than the other.

Wire-up:

- [spec.py](apps/compounds/nlp/spec.py): `FIELD_METRIC = "metric"`; `MetricMiss(query, available_metrics, filter_index, field)`; `MetricMiss` joins `ExecutionResult`.
- [view.py](apps/compounds/nlp/view.py): `MetricMiss` joins the miss-branch dispatch; `suggestions=[]` is added to the body so the frontend's uniform miss renderer doesn't trip on the extra-fields-only payload.
- Frontend: `FIELD_METRIC` constant, `'metric'` joins `CandidateKind`, `available_metrics?: string[]` on the miss-branch type, `MissView` renders the available KPIs as chips with the helper line *"Available KPIs in scope:"*.

Tests: 5 evaluator (lenient case + punctuation, `metric=None` two ways, plus the existing strict-EC50-vs-IC50 mismatch is preserved); 4 executor (`MetricMiss` with available KPIs, `filter_index` tagging on the second filter, lenient match end-to-end, empty-scope case where `available_metrics == []`); 1 view test asserting the response shape end-to-end. **271 passing** (was 262).

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
| 6 | Pinning IDs live on `QuerySpec` (`registration_target_id`, `assay_target_id`, `protocol_id`), not on a separate `ResolvedSpec` type | The LLM never emits these; only the view populates them on continuation. Putting them on `QuerySpec` keeps the continuation path through the *same* `execute()` entrypoint as fresh queries — no "resolved-entity second-stage" API to maintain. Cost is three nullable fields the LLM ignores. |
| 6 | Pinned protocol UUID is scope-validated — the resolver checks it against the resolved-target scope before accepting it | A stale clarify response from a prior resolved scope shouldn't be able to make the executor run a protocol from a different target. If the pinned id isn't in scope, return `ProtocolMiss` (same UX as a genuine no-match). |
| 6 | Per-project target fill (decision 18) is **not** implemented in slice 6 | No per-project compounds pages exist in the frontend today. Decision 18 presumes the request carries a project context signal. Implementing the backend half without a caller is dead code; the fill will land when slice 7 introduces per-project pages and has a concrete signal to send. |
| 6 | Daily cap uses Django's default `LocMemCache` (per-process) not a shared Redis/DB store | Decision 17 frames the cap as "catch retry storms / buggy integrations", not cost protection. At 500/day/user-per-process × N workers, a user hitting the absolute worst case still runs orders of magnitude below any plausible human workload. Upgrade to a shared counter when the cap purpose shifts to cost gating. |
| 6 | Response shape is a single flat object discriminated by `status` | Matches the "single JSON envelope with a kind tag" convention that keeps frontend fetchers simple — no union types on the client, just a `switch (resp.status)` at the render boundary. HTTP status code is a secondary signal: 200 covers successful dialogue (table/clarify/miss/not_a_query), 400 covers client-correctable errors, 429 rate limit, 502 LLM-glue failure, 404 feature-off. |
| 6 | `NotAQuery` maps to HTTP 200 (with `status: "not_a_query"`), not 4xx | It's the LLM behaving correctly, not an error. UI shows the LLM's reason string plus an example-prompts panel (decision 20). Treating it as 4xx would conflate a valid response with a transport failure. |
| 6 | `ParseError` maps to 502, not 500 | 502 ("bad gateway") is the HTTP idiom for "an upstream service we depend on returned something we couldn't process". `ParseError` means the Azure OpenAI response was malformed — it's the upstream's output that failed. Clients can distinguish "retry the same request" (502, transient) from "the code is broken" (500). |
| 8 | Golden set is YAML with entry-type discriminator, not separate files per type | One file keeps the curated corpus skimmable as a single list — reviewers can see "how many cross-project prompts do we cover?" by eye. Commented section headers within the file organise it without fragmenting the corpus. |
| 8 | The offline test round-trips each golden entry through `_to_parse_result` (using the entry's own expected shape as simulated LLM output) | Verifies the YAML is actually reachable from a well-formed LLM response — any drift between the YAML schema and the Python translator surfaces immediately as a test failure, not at first online eval. |
| 8 | Online eval is skipped by default, gated on `RUN_NLP_ONLINE_EVAL=true` | Every unskipped test in this repo is free to run. Online eval burns real LLM calls and cannot run in default CI. Default-off + explicit opt-in makes the behaviour legible. |
| 8 | Bar thresholds (≥90% spec, 100% not_a_query, 0 parse errors) are assertions, not warnings | Bar misses should fail loudly — the eval is the regression harness for future LLM / prompt changes. Soft warnings would let drift accumulate. |
| 8 | Diff rendering on mismatch is per-field `"; "`-joined rather than a generic dict diff | Every failure mode points at a specific field (e.g. `threshold.unit: expected 'uM', got 'µM'`) — far more actionable than "expected {…} got {…}" when debugging prompt-engineering regressions. |
| 7 | NLP-specific fetcher (`postNlpQuery`) returns the body on any HTTP status | The existing `apiPost` throws on non-2xx, which is right for CRUD endpoints but wrong here — 404 (feature-off), 400 (scope/spec error), 429 (cap), 502 (parse error) all carry structured `NLPResponse` bodies the UI wants to render. Cheaper to write a ~20-line wrapper than to add error-body propagation to the shared helper. |
| 7 | Discriminated union rendered via `switch (response.status)` in `NLPResults.tsx`, not via per-status sub-routes or context providers | Five branches, no shared state across them, each tested in isolation by shape. A state machine or reducer would be over-engineering for v1. The branches can be extracted into separate files later if any one of them grows. |
| 7 | Frontend reuses the aggregation module's `PROPERTY_LABELS` for column headers | Keeps "MW", "cLogP", "TPSA" etc. consistent across the aggregation table and the NLP table without a second source of truth. When the aggregation refactor added new properties, NLP picks them up automatically. |
| 7 | Example-prompt buttons below the command bar double as the example-set for the `not_a_query` branch | Decision 20 mandates showing example prompts when the LLM rejects a request. Serving the same examples as discoverability affordance on an empty page (and rejection recovery on a misfire) means the user learns one set and sees it in both contexts. |
| 7 | Per-project placement deferred until per-project pages exist | Matches the §18.3 slice-6 decision on backend target-context fill — building the UI half of a feature whose counterpart does nothing yet would rot. Reopens when per-project compounds pages land. |
| 7 | `partial_spec` on clarify responses was a slice-6 hole caught during slice-7 wiring; fixed alongside slice 7 with a continuation-round-trip view test | The slice-6 continuation story assumed the client somehow had the partial_spec to round-trip — but the clarify response didn't include it. Fixing it in a slice-6 maintenance pass bundled with slice 7 kept the fix reviewable alongside its first real caller (the `ClarifyView` chip `onClick`). |
| 9 | **Pivot shape**: NLP produces a compound selection, aggregation page produces the display | The aggregation page already has a rich URL contract (`targets` / `compound` / `protocols` / `format` / `aggregations` / …) and an automatic project-card / spider-plot mode. Duplicating any of that on the NLP side was an architectural mistake from v1 — it inevitably lagged behind what the chemist could already do via clicks. Pivot makes NLP do one job well (compound selection) and hand off to the specialist for display. |
| 9 | `measurement_filters` is a list with AND semantics, not a single top-level protocol/metric/threshold | Cross-protocol selectivity is not a fringe case — it's how kinase programmes, degrader selectivity, counter-screens are all phrased. A list of filters intersected makes the headline v2 capability fall out of the shape; cross-protocol is no longer a "feature to add" but a trivially-expressible query. |
| 9 | Default redirect `format=cards` (not `compact`) | `format=cards` renders project cards with auto-populated spider plots — the visual-triage view chemists actually use for "find the right compounds". A column-table default would have undersold the capability. The aggregation page's own format-switcher is one click away if the user wants otherwise. |
| 9 | Executor emits `protocol_names` from the filters (not from the output intent) | Pre-pivot, `protocol_hint` was both filter AND output. Post-pivot, only filters carry protocols; the output-axis concept goes away. `protocol_names` passed in the URL pre-populates the aggregation page's protocol-column picker so the user can see filtered-vs-other columns side by side after landing. |
| 9 | Click-to-confirm redirect (not auto-redirect) | User's call on the pivot's final UX. Auto-redirect when the LLM is confident, click when the user should verify. At v2, the LLM's confidence is unproven against real chemist prompts; click gives the user a chance to see "Found 0" or "Found 843" and decide before navigating. Can flip to auto for small-N results once eval data suggests it's safe. |
| 9 | `filter_index` tagged on `ProtocolClarify` / `ProtocolMiss` (not on the selector itself) | A clarify for filter #2 shouldn't look like a clarify for filter #0. The backend tags the filter index; the frontend (`applyClarifyPick(filterIndex)`) knows which `measurement_filters[i]` to pin `protocol_id` on. Avoids the alternative of requiring the client to scan filters for a matching hint. |
| 9 | Filter without a threshold means "has any valid measurement" (not "all compounds are eligible") | Chemists genuinely want *"compounds tested in HTRF"* — compound must have at least one valid HTRF row. A filter with null threshold naturally expresses this. The alternative (treat null-threshold filter as no-op) would silently drop filters from the intersection. |
| 9 | Window.location.assign for the redirect, not Next.js's router.push | The redirect target is on the same origin but a completely different page (aggregation, not nlp). A full navigation is cleaner than a client-side route change that'd re-mount layouts differently. Also side-steps any SSR/caching weirdness on the aggregation side. |
| 13 | Scaffold catalog lives as a curated Python module (`substructures.py`), not a database table | §19.3 flagged a `Scaffold` DB table as one option. The curator workflow is **"edit the tuple, commit, deploy"** — no migration, no retag, no re-indexing. SMARTS refinements take effect instantly for every compound on the next query. DB-backing would add a migration per edit and a service call per resolve for no scale benefit at DDU volume. Escape valve (FP pre-screening) can be a sidecar table later without disturbing the catalog. |
| 13 | Pinning ID for scaffolds is the canonical name, not a UUID | Canonical names are already unique and stable. Assigning UUIDs would require a registration/versioning story for no gain. `scaffold.name == scaffold.id` keeps the clarify round-trip trivial and the URL (when we add it to `/assays/aggregate`) human-readable. |
| 13 | Substructure filtering happens **before** measurement filters | Cheap compound-level predicate: matching ~thousands of SMILES against one SMARTS is sub-second, and the resulting compound-ID set can be intersected with the measurement-filter set without a second RDKit pass. Puts the narrowing early in the pipeline where it narrows the scope for everything downstream. |
| 13 | RDKit import is lazy (inside `_match_compounds_by_scaffold`) | Prompts without a `scaffold_hints` entry skip the import entirely. Keeps the executor's hot path free of RDKit cold-start cost and makes the rest of the `compounds.nlp` package safe to import in environments where RDKit is absent (test containers, minimal images). |
| 13 | Resolver three-tier match mirrors Target's (exact → substring → fuzzy miss) rather than Protocol's (strict/substring/overlap) | Scaffold names are single tokens (*"pyrimidine"*, *"sulfonamide"*) — the Protocol resolver's token-set containment doesn't apply. But typos + prefix abbreviations do, so `SequenceMatcher` ratio on name + aliases provides the fuzzy tier without needing a full multi-token ranking engine. |
| 13 | Multiple scaffolds ANDed as parallel `scaffold_hints` + `scaffold_ids` arrays, not as a single nested list of {name, id} pairs | Simpler JSON shape (two arrays the LLM can emit independently) and cheaper clarify-continuation — picker returns one scaffold, frontend writes into the indexed slot. The "parallel arrays indexed together" discipline is identical to the existing `measurement_filters[i]` pinning pattern. Frontend helper `applyScaffoldPin` handles the pad-to-length invariant. |
| 13 | Scope sentence renders multiple scaffolds as `"containing pyridine AND amide"` | Makes the AND semantics explicit in the echo-back. A user reading *"containing pyridine, amide"* might expect OR; the uppercase AND is disambiguating. Same convention as the existing multi-filter *"where HTRF IC50 < 100 nM AND TM IC50 > 1 uM"* rendering. |
| 13 | LLM emits aliases / plurals verbatim (*"pyrimidines"*, *"pyridinyl"*) — resolver normalises | Preserves the decision 2 discipline (LLM never normalises canonically; backend owns the catalog). The alias-pool in `substructures.py` knows that *"pyrimidines"* → *"pyrimidine"* and that *"pyridinyl"* → *"pyridine"*; pushing that knowledge into the LLM prompt would bloat it and re-centralise chemistry knowledge in the model. |
| 14 | Compound pinning is **additive (UNION)**, not narrowing | Every other predicate narrows the selection. Pins are different — the chemist's intent is "show this lead alongside whatever else I'm asking about", not "show only compounds that match this AND the other criteria". Implementing as UNION keeps the use case ("lead-as-benchmark for visual comparison") natural; an AND interpretation would just be "show me the lead, if it happens to match the filter" which is rarely what's wanted. The trade-off is the first non-narrowing operation in the executor, but the shape stays clean. |
| 14 | Pin reuses the existing `compounds.formatting.extract_reg_number` helper rather than rolling a new parser | That helper already handles every prefix variant (full / compact / spaced / underscored / malformed / bare integer) and is unit-tested across the codebase. Slice 14's resolver is ~10 lines on top of it. |
| 14 | No clarify path for compound refs — only Resolved / Miss | Compound IDs are unique integers — there's nothing fuzzy to disambiguate. A typed reference either parses to an existing compound or doesn't; the user revises the prompt on a miss. Saves the dataclass-and-frontend-chip dance that targets / protocols / scaffolds need. |
| 14 | `ref_index` tagged on `CompoundMiss` (mirroring `filter_index` on protocol clarify and `scaffold_index` on scaffold clarify) | When the LLM emits a list of pins, the user needs to know which entry went wrong. Tagging the index keeps the miss legible — *"the third pin you typed didn't parse"* rather than a generic "one of your pins missed". |
| 14 | Pin-only selector short-circuits the unscoped-base to empty set (rather than `Compound.objects.all()`) | Naive "unscoped + UNION pins" returns the entire registry plus the pin (a no-op union), which is wrong: the user pinned ONE compound and expects to see ONE compound. The `_selector_has_non_pin_narrowing` helper distinguishes "the selector has predicates that narrow the unscoped base" from "the selector is only pins" and skips the all-compounds default in the latter case. |
| 14 | Pin survives ALL narrowing — even a strict measurement filter or out-of-scope target | The lead-as-benchmark use case demands this: a chemist scrolling their recent ARd compounds with the lead pinned wants the lead visible REGARDLESS of whether the lead passes the strict cutoff or even sits in ARd. UNION-after-all-narrowing implements this naturally; UNION-before would let later narrowing accidentally drop the pinned compound, defeating the purpose. |
| 14 | LLM emits the typed reference verbatim (e.g. *"compound 26007"* → emit `"26007"`, not `"NCL-00026007"`) | The LLM-as-parser discipline (decision 2) — never canonicalise. Bare integer, full prefix, lower-case, malformed-with-extra-zeros: all valid LLM emissions; the backend's parser is one helper call away from the resolved Compound row. |
| 15 | KPI match is **lenient** (case + non-alphanumeric punctuation insensitive), not strict | The pre-fix strict match silently dropped chemist-equivalent KPIs that differed only in case/punctuation (*"IC50"* vs *"ic50"* vs *"IC-50"*). Re-using the existing `compounds.utils.normalize_ref` keeps the leniency consistent with how target / scaffold names are matched elsewhere — same discipline, same canonicalisation function. |
| 15 | The user's typed metric is the **search key**; the row's stored KPI is the **value-lookup key** | Once the lenient match establishes equivalence, the value column in `results` is keyed off whatever the row actually stored (`results["ic50"]`, not `results["IC50"]`). The chemist's typing doesn't have to match the column casing in the JSON blob — only the underlying KPI identity. |
| 15 | `metric=None` means *"any numeric KPI counts; look up via the row's own KPI key"* | Documented in the system prompt; pre-fix the evaluator implementation contradicted that ("filter every row with a KPI string"). Aligning the implementation with the prompt's contract is a behaviour fix, not a feature change. |
| 15 | Metric pre-check at the executor level, not the evaluator | The evaluator is per-row — it can't tell the difference between "no rows matched the threshold" and "no rows matched the metric". The executor sees the full scope and is the right place to surface "the typed metric isn't recorded anywhere in scope" as a Miss. |
| 15 | `MetricMiss` carries `available_metrics: List[str]`, not generic candidate suggestions | The "did you mean…" answer for KPIs is a list of the actual stored values, not a fuzzy-ranked suggestion. A flat string list is the simplest payload; the frontend renders chips off it. |
| 15 | First filter with a metric mismatch short-circuits (rather than collecting all misses) | Mirrors the existing per-filter clarify pattern (protocol clarify, scaffold clarify): one Miss at a time, `filter_index` tagged so the UI knows which filter the user is editing. Collecting and surfacing multiple misses simultaneously is doable but adds payload complexity for unclear UX gain. |
| 15 | Empty `available_metrics` (no rows in scope at all) is *also* a `MetricMiss`, not a separate empty-scope error | A user typing *"ARd FOO IC50 < 10 nM"* against an empty programme should see a Miss with empty `available_metrics`, which the frontend renders as "no measurements in scope under any KPI — try a different protocol". One code path, one component branch — and the right answer for the user is the same shape. |

### 18.4 Running the tests

From `server/`, with CCP4 sourced:

```bash
source ../../ccp4-20251105/bin/ccp4.setup-sh
PYTHONPATH="$PWD:$PWD/../apps" DJANGO_SETTINGS_MODULE=compounds.settings \
  ccp4-python -m pytest ../apps/compounds/nlp/tests/ -v
```

Expected: 193 passing + 2 skipped as of slice 8 (2026-04-23). The 2 skipped are the online eval tests in `test_golden.py` — they run only when `RUN_NLP_ONLINE_EVAL=true` and `AZURE_OPENAI_ENDPOINT` is set, because they burn real LLM calls. Offline tests (evaluator, LLM, golden-set validation/round-trip) do not touch the DB. LLM tests mock the Azure OpenAI factory so they don't require `openai` to be installed locally. Resolver / orchestrator / protocol / executor / view tests use `@pytest.mark.django_db` via the `db` fixture argument. View tests use DRF's `APIClient.force_authenticate` and patch `compounds.nlp.view.parse_prompt` to mock the LLM without touching `openai`. Fixtures are defined inline in each test module (no JSON fixture files, no separate conftest).

### 18.5 Pending work (post-v1)

All eight slices are shipped. Remaining non-code work to close out the v1 rollout:

- **§16.1 target-gene backfill on Demo and Kawamura** — before the feature is turned on for those instances, run `backfill_target_genes.py --dump` per instance and curate the plans. DDU is already fully backfilled.
- **Online-eval run against the live DDU endpoint** — `RUN_NLP_ONLINE_EVAL=true` with the golden set at `apps/compounds/nlp/tests/golden/prompts.yaml`. First run will surface which prompt families the live gpt-4o deployment handles cleanly vs. needs prompt-engineering tuning.
- **Per-instance enablement** — each instance flips to the feature by setting `AZURE_OPENAI_ENDPOINT=https://ddu-openai.openai.azure.com/` + `COMPOUNDS_NLP_ENABLED=true` in its `.env.*` file and redeploying server. Already done for DDU (2026-04-23).

Post-v1 feature directions live in §19 (chemistry-aware queries, scatter, compound-rooted queries, missing-data, explain-this-query). §19.7 names the top picks.

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

**Shipped in slice 13 (2026-04-23).** The design below was the pre-implementation sketch; the landed shape matches it closely, with two refinements noted after.

The most naturally-phrased chemistry query class — *"HTRF values of ARd compounds containing pyrimidine"* — fits the v1 architecture with **no change to the LLM's role**. The LLM still just emits a name (*"pyrimidine"*); the backend owns the authoritative chemistry catalog.

Implementation shape:

- **New resolver.** `resolve_scaffold(hint) -> ScaffoldResolution` follows the `Resolved / Clarify / Miss` triad established in §18.2. Matches against a curated `Scaffold` table (`name`, `aliases`, `smarts`, optional `svg`), edited by a curator the same way `Gene` rows are edited today. Clarify when a fragment name is ambiguous (*"pyrimidine"* vs *"pyrimidinone"* vs *"fused-pyrimidine core"*); miss with visual suggestions.
- **New executor predicate.** Substructure match is a per-row predicate in the same bucket as unit-threshold. `FILTER_SUBSTRUCTURE_ABSENT` slots next to `FILTER_THRESHOLD_NOT_MET`; `EXCLUDE_STRUCTURE_UNPARSEABLE` handles SMILES parse failures at the row level.
- **New `QuerySpec` field.** `scaffold_hint: Optional[str]` — the user's typed fragment name, verbatim.
- **No SMARTS in the prompt, ever.** The LLM never sees SMARTS and never emits SMARTS. Same discipline as "no catalogs in the prompt" for targets/protocols (§8, decision 2). This eliminates the class of failures where the LLM hallucinates structurally-plausible-but-wrong SMARTS.

**Scale path.** Naive per-row `rdkit.Chem.Mol.HasSubstructMatch` is fine at DDU-scale (low thousands of compounds) and should be the first slice. RDKit pattern fingerprints in a sidecar table add a screening layer for larger libraries (SSS via FP bit-subset match before exact substructure). The Postgres RDKit cartridge is the proper endgame for very large libraries but needs a prod/SQLite-tests split story. **Build naive first; add FP screening only when a real library hits the wall** — premature FP infrastructure is the canonical trap in this direction.

**Boundary with generative SAR prose.** *"Compounds containing pyrimidine with IC50 < 10 nM"* is tabular — fits this bucket. *"What's driving IC50 in the pyrimidine series?"* is generative SAR — belongs in bucket 5 (generative analysis) with a different LLM data-visibility model. Keeping the two in separate proposals matters.

**Refinements in the landed shape** (see §18.2 slice 13 entry for the full write-up):

- **Catalog is a Python module, not a DB table.** The above sketch mentioned a `Scaffold` DB table; the landed catalog is `apps/compounds/nlp/substructures.py`. Curator edits the module, commits, deploys — no migration. DDU scale doesn't justify the DB-backing overhead; FP sidecar (if it lands) would be an orthogonal addition.
- **Scope is both selectors, not a `MeasurementFilter` field.** `scaffold_hints` lives at the `CompoundSelector` / `AssaySelector` level (applies to *which compounds qualify overall*), not per-filter. Multiple scaffolds are ANDed. This matches how chemists actually phrase the predicate (*"pyridines that are potent in HTRF"*, not *"compounds potent in a pyridine-specific HTRF"*).

### 19.4 Chemistry queries collapse into another resolver

A useful framing: **chemistry-aware queries are structurally another resolver**, not a new subsystem. The substructure resolver carries the same `Resolved / Clarify / Miss` triad, the same field-tagging for the orchestrator, the same LLM-emits-name-backend-owns-catalog discipline. What differs is the matching engine (RDKit `HasSubstructMatch` instead of string normalization or token-set containment) and the candidate payload shape (SVG drawing + SMARTS in place of gene symbols or `n_runs`). The executor doesn't grow new surface; it just handles one more optional `QuerySpec` field.

**The family this seeds.** Compound-by-ID resolver (*"compound NCL-00026042"*) — **shipped in slice 14** as `resolve_compound_ref`; similarity resolver (*"structurally similar to X"* — Tanimoto k-NN with a tuneable threshold); scaffold-by-project resolver (project-local series curated as `Scaffold` rows linked to a `Target`). All cast into the same mould.

**When to abstract.** Four-plus concrete resolvers following the same shape would justify lifting a generic `Resolution[T, C]` protocol. We are at five today (target / protocol / user / scaffold / compound-by-id), but the compound-by-id resolver landed without a clarify branch — IDs are deterministic — which means the canonical `Resolved / Clarify / Miss` triad doesn't actually fit it. The resolvers genuinely differ in candidate payload shape and in whether clarify is sensible at all; the pressure to abstract that §19.4 anticipated has not actually materialised. Default position holds: **wait for the push, don't pre-empt the abstraction.** The cost of premature generality here is mostly design-lock-in around candidate payload shapes that are genuinely different per resolver, for the benefit of saving a few lines of type declarations.

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

### 19.7 First picks after v2 pivot

**Cross-protocol selectivity is now SHIPPED** in slice 9 — it's a natural consequence of the measurement_filters list shape. Two remaining top-priority post-pivot directions:

- **Compound-rooted queries** — **date range filters shipped in slice 10, user filters shipped in slice 11.** Selector gains `registered_date_range` + `registered_by_as_typed`; each measurement filter gains `assay_date_range` + `assayed_by_as_typed`. User resolver has its own clarify/miss triad (ambiguous "Alice" surfaces a chip picker with display name + email + compound/assay counts; pool scoped to users who actually have provenance data). LLM resolves relative dates via `[Today: YYYY-MM-DD]` on the user message. Remaining compound-rooted work: `missing_measurement_for` (anti-join — compound has NO matching row — one-liner in Django ORM, shape TBD).
- **Explain-this-query.** One-click reveal of the resolved `CompoundSelector`, scope sentence, and the redirect URL itself. Every chemist who uses an LLM over their data asks *"what did it actually do?"* The v2 pivot already does 80% of this — the redirect URL IS the explanation (bookmarkable, editable, shareable). Exposing the resolved selector + scope sentence alongside the click-to-redirect button is mostly a UI lift. Surfaces the LLM's interpretation before the user commits to navigating.

Post-pivot, two directions that WERE on the list naturally fold into either the shipped selector (selectivity) or are queued as field-additions (compound-rooted). **Chemistry-aware substructure is now SHIPPED in slice 13** (see §18.2; the landed shape diverged from the pre-implementation sketch by curating the catalog as a Python module rather than a DB table, and by attaching `scaffold_hints` at the selector level rather than per-filter). **Compound-ID pinning is now SHIPPED in slice 14** (the first additive / UNION operation in the selector — see §18.2 entry; landed as `compound_refs_as_typed: List[str]` on `CompoundSelector`, deterministic resolver via the existing `extract_reg_number` helper, no clarify path).

The remaining forward-looking buckets are:

- **NLP-driven scatter** (§19.5). Gated on the existing aggregation-page scatter component being easy to drive from a structured payload. If that integration is cheap (likely), this may well be the single most compelling demo — multi-axis data visualisation with narrative highlights — and is probably the cheapest high-impact win.
- **`missing_measurement_for`** — anti-join on the compound-selector path: *"ARd compounds with no HTRF measurement yet"*. One-liner in Django ORM once the shape lands. A new `MeasurementFilter.exists_mode` enum (`any` default, `none`) is the minimal shape; alternatively a top-level `missing_filters: List[MeasurementFilter]`. Sequence either ahead of or after scatter depending on chemist demand.
- **Explain-this-query.** One-click reveal of the resolved `CompoundSelector`, scope sentence, and the redirect URL itself. 80% already there via the redirect URL; the UI lift surfaces the resolved selector alongside the click-to-redirect button.

### 19.8 Deliberate non-goals for v2

To protect against scope creep in the v1 → v2 transition:

- **Generative SAR narrative** needs a separate proposal — different data-visibility model for the LLM, different prompt architecture, different evaluation harness. Do not graft it onto this one.
- **Chart code generation by the LLM** (§19.5 boundary) — the LLM emits a `ChartSpec`, never plotting code. This is an architectural non-negotiable for the scatter extension.
- **Write-capable commands** remain locked out by decision 11 until a trust/audit model is specified in its own proposal.
- **Cross-session conversational memory** is a distinct problem from extending the `QuerySpec` schema and should not be conflated with "more query shapes". Revisit when v1 usage data justifies a chat-style surface.
