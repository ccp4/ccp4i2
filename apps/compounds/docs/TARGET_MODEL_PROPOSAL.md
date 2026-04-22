# Target Model — Gene Symbols, Aliases, and Cross-References

**Status**: Draft for review. Design captured in dialogue during work on the NLP query proposal (§6.1 of `NLP_QUERY_PROPOSAL.md`); implementation not yet started.
**Author**: Martin Noble (with Claude)
**Date**: 2026-04-22
**Companion to**: `NLP_QUERY_PROPOSAL.md` (primary downstream consumer of this work)

## 1. Motivation

The `Target` model today carries a single free-text `name` field as its sole identifier:

```python
class Target(models.Model):
    name = models.CharField(max_length=255, unique=True)
    parent = models.ForeignKey('self', ...)          # hierarchy
    image = models.ImageField(...)                   # dashboard branding
    saved_aggregation_view = models.JSONField(...)   # UI state
    # + audit fields
```

That works for display but is the wrong shape for matching and lookup. Four concrete problems:

1. **Name is a poor matching key.** "AR degraders" and "AR" and "Androgen Receptor" and "NR3C4" are all things a chemist might type when they mean the same target. Today none of them equal the canonical `name` string, so lookup fails or relies on ad-hoc substring hacks.

2. **Multi-gene targets are invisible.** Protein–protein interaction targets (Skp2-Cks1, MDM2-p53, BCL-XL/BAX, cyclin-CDK complexes) span two or more genes. A compound flagged as a Skp2 binder in a counter-screen should resolve to the Skp2-Cks1 PPI target. Today there is no structured way to express that the PPI target *involves* SKP2 and CKS1 as genes.

3. **No external cross-references.** UniProt accession, Ensembl gene ID, HGNC ID — all absent. These are the keys that make every other cheminformatics tool interoperable (deposition, ELN integrations, external SAR databases, literature mining).

4. **Annotator burden is poorly managed.** A richer target model is only worth building if it doesn't turn target creation into a data-entry chore. The right shape is: annotator provides the **gene symbol(s)** (and not much else); the system hydrates aliases and cross-references automatically.

The immediate pressure comes from the NLP query work — matching a user's typed target against a deterministic rule needs a richer pool than a single canonical name — but the payoff is broader: every downstream feature that needs to reason about biological identity benefits (SAR search, cross-target annotations, external ELN sync, publication deposition).

## 2. Scope

| In scope | Out of scope (for now) |
|----------|------------------------|
| Extending `Target` with gene symbols, aliases, HGNC/UniProt cross-refs | Rewriting the aggregation UI or dashboards |
| A new `Gene` model holding HGNC-hydrated metadata | A bespoke gene/target picker component (reuse what's there) |
| HGNC as the canonical source of human gene metadata; UniProt as a second layer for protein cross-refs | Non-human targets (MGI for mouse, etc.) — revisit if cross-species targets land |
| A hydration pipeline (on save, or batch) that populates aliases from external APIs | Full biocuration (pathways, GO terms, disease associations) |
| A target-creation workflow that minimises annotator burden | Ontology-level reasoning (kinase family membership, etc.) |
| Data migration for existing `Target` rows (opt-in backfill, not enforced) | Variant / mutant representation — deliberately deferred to v1.1 (see §9) |
| A local-snapshot fallback for zero-external-call operation | Full offline-first design |

## 3. Current state

- `Target.name` is the only external-facing identifier. It's free-text, unique across the table, and currently carries everything an annotator might want to express (including, informally, variant/mutation info baked into the string).
- `Target.parent` exists as a self-FK for hierarchical organisation. Under-used today, but well-positioned to represent variant-of relationships later (§9).
- `Compound.target` is PROTECT/required — every compound belongs to exactly one target.
- `Assay.target` is SET_NULL/nullable — an assay may or may not record which target it was run against.
- There are no cross-references to HGNC, UniProt, Ensembl, or any other external system.

## 4. Proposed additions

A new `Gene` model and a many-to-many linkage to `Target`:

```python
# apps/compounds/registry/models.py

class Gene(models.Model):
    """
    HGNC-hydrated metadata for a single gene.

    Canonical by symbol. HGNC is the source of truth for human gene
    nomenclature; this table caches it so downstream matching does not
    depend on runtime external calls.
    """

    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    symbol = models.CharField(max_length=32, unique=True,
        help_text="HGNC-approved gene symbol, e.g. 'EGFR', 'SKP2'.")
    hgnc_id = models.CharField(max_length=32, blank=True, db_index=True,
        help_text="HGNC stable ID, e.g. 'HGNC:3236'.")
    name = models.CharField(max_length=255, blank=True,
        help_text="Full gene name from HGNC, e.g. 'epidermal growth factor receptor'.")
    aliases = models.JSONField(default=list, blank=True,
        help_text="Union of HGNC alias_symbol and prev_symbol lists.")
    uniprot_ids = models.JSONField(default=list, blank=True,
        help_text="UniProt accession(s); typically one primary + isoforms.")
    ensembl_gene_id = models.CharField(max_length=32, blank=True, db_index=True)

    hydration_status = models.CharField(max_length=16,
        choices=[('pending', 'Pending'),
                 ('ok', 'OK'),
                 ('partial', 'Partial'),
                 ('failed', 'Failed'),
                 ('manual', 'Manual override')],
        default='pending')
    hydrated_at = models.DateTimeField(null=True, blank=True)
    hydration_source = models.CharField(max_length=32, blank=True,
        help_text="'hgnc-api', 'hgnc-snapshot', 'manual', etc.")

    created_at = models.DateTimeField(auto_now_add=True)
    modified_at = models.DateTimeField(auto_now=True)


class Target(models.Model):
    # ... existing fields unchanged ...

    genes = models.ManyToManyField(
        Gene,
        related_name='targets',
        blank=True,
        help_text="Biological gene(s) underlying this target. "
                  "Multi-gene for PPI / complex / pathway targets."
    )
```

### 4.1 Why a through table and not an ArrayField

Three reasons to normalize `genes` as a proper M2M rather than stuff `["SKP2","CKS1"]` into an `ArrayField` or `JSONField` on `Target`:

- **Shared hydration.** Two targets (say "EGFR-WT" and "EGFR-L858R selective") both reference the EGFR gene. HGNC hydration for EGFR happens once, benefits both. With per-target arrays, you'd duplicate aliases everywhere.
- **Cross-target queries by gene.** "Find all targets involving SKP2" becomes `Gene.objects.get(symbol='SKP2').targets.all()` — no array containment tricks.
- **Portability.** `ArrayField` is Postgres-only; M2M + JSONField works on SQLite (used in tests) and Postgres (used in deployment) alike.

### 4.2 Why `Gene.aliases` as JSONField, not its own table

Aliases are a read-mostly flat list; they're identifiers that matter only when you're *looking up* a gene. They have no metadata worth modelling independently (when did this symbol become an alias? who says so? — that's HGNC's problem, not ours). JSONField is the right level of structure.

## 5. Hydration from HGNC / UniProt

### 5.1 HGNC API

HGNC (`rest.genenames.org`) is free, no-auth, and authoritative for human gene nomenclature. Single-gene lookup:

```
GET https://rest.genenames.org/fetch/symbol/EGFR
```

returns a JSON document with `symbol`, `name`, `alias_symbol` (list), `prev_symbol` (list), `uniprot_ids` (list), `ensembl_gene_id`, and a stable `hgnc_id`. That's every field `Gene` needs in one call.

### 5.2 UniProt (optional second layer)

If `Gene.uniprot_ids` is populated, an optional second call to UniProt:

```
GET https://rest.uniprot.org/uniprotkb/{accession}.json
```

can add protein-level aliases (recommended name, alternative names, submitted names). For v1 this is probably unnecessary — HGNC already gives us most of what we want — but the hook exists for later enrichment.

### 5.3 When hydration runs

Two triggers:

- **On `Gene` save** (via `post_save` signal) if `hydration_status == 'pending'` or `hydrated_at` is stale (> 90 days). HGNC responds in <500ms for single-gene lookups; synchronous is fine.
- **Nightly management command** `hydrate_genes` that re-fetches any `Gene` older than a configurable threshold. Idempotent; safe to re-run.

Errors are captured in `hydration_status = 'failed'` with the reason in a log line; the `Gene` row still exists and is usable with whatever metadata was typed in manually.

### 5.4 Local snapshot fallback

For deployments that want zero external runtime calls, HGNC publishes a **complete dataset** (~50k rows) as JSON/TSV on a regular cadence:

```
https://storage.googleapis.com/public-download-files/hgnc/json/json/hgnc_complete_set.json
```

A management command `import_hgnc_snapshot` can load the full dataset locally; the hydration code then reads from the snapshot table rather than making API calls. This is overkill for v1 but eliminates any concern about outbound traffic patterns.

## 6. External exposure analysis

HGNC / UniProt queries are pull-only, from our IP, with no authentication linking them to our internal data. The servers see queries for "EGFR" alongside millions of other queries per day; there is no realistic channel for them to correlate our query pattern with anything meaningful about our Target list. Other channels (commit logs, deployed config, user roles) leak more information than an HGNC query ever could.

If stricter isolation is required:
- Use the local-snapshot fallback (§5.4) — zero runtime external calls.
- Or route hydration through a dedicated outbound egress with no association to application traffic.

**Recommendation**: v1 uses live HGNC API; add snapshot fallback as an optional management command available to deployments that want it.

## 7. Annotator workflow

The goal is to make target creation feel lighter, not heavier, despite the richer backend model.

Target creation form (existing + proposed):

| Field | Current | Proposed |
|---|---|---|
| Name | Required, free text | Required, free text (unchanged) |
| Parent | Optional FK | Optional FK (unchanged) |
| Image | Optional upload | Optional upload (unchanged) |
| **Gene symbols** | *(does not exist)* | **Required for new targets; autocompletes against `Gene.symbol`; allows "create new" if symbol not in table** |

New-target save flow:
1. Annotator types target name and gene symbol(s).
2. Each gene symbol: if it already exists as a `Gene`, link it; if not, create a `Gene(symbol=..., hydration_status='pending')` and trigger hydration.
3. `Target.genes.set([...])`.
4. Hydration fires asynchronously (or synchronously if fast enough); the UI can show a spinner next to the gene chip and fill in the full name once hydrated.

**Nothing beyond the gene symbol is asked of the annotator.** Everything else (full name, aliases, UniProt IDs, HGNC ID) is filled in by the system.

### 7.1 Manual override

Some targets don't fit HGNC — internal project codes, consortium-defined entities, or genuinely non-gene targets (formulation-focused projects, delivery-focused work). For these:

- `Target.genes` can be empty (`blank=True`).
- A target with no genes is fine; NLP matching falls back to `Target.name` only for that row.
- An annotator can set `Gene.hydration_status = 'manual'` and edit the fields directly if they want to enter curated data for a non-HGNC entity.

## 8. Backfill of existing targets

Existing `Target` rows will have `genes=[]` after the migration. Three options:

- **(a) Leave them alone.** They keep working; NLP falls back to name-only matching for them. Acceptable if the existing target list is small and the maintainers will fill in genes when they next touch each target.
- **(b) Prompt in the UI.** Target detail page shows a banner when `genes` is empty: *"This target has no gene symbols recorded. Add them to enable richer search and external cross-refs."*
- **(c) Heuristic auto-suggest.** Run `Target.name` through a substring-match against `Gene.symbol` — if the name contains a known gene symbol, suggest that gene with a one-click confirm.

**Recommendation**: (a) + (b) for v1. (c) is tempting but risks bad auto-matches (target named "ARD" could be misidentified as ARG1 or AR or ARID1A); keep it human-in-the-loop.

No forced migration, no data loss, no user-visible regression.

## 9. Variants / mutants — deferred to v1.1

EGFR (L858R, T790M, C797S) is the canonical case. Other examples: KRAS-G12C, BRAF-V600E, IDH1-R132H, each with drug discovery programmes distinct from the wild-type target.

Options considered:

| Option | Shape | Trade-off |
|---|---|---|
| (a) Separate `Target` per variant | "EGFR-WT", "EGFR-L858R" as distinct rows | Simple, loses family semantics; noisy target list |
| (b) `Target.variant` CharField + `Target.mutation_codes` | "L858R", `["L858R"]` | Clean, small migration, composable |
| (c) Variant on `Assay` only | Compound is registered against "EGFR"; assay specifies variant | Loses registration intent (3rd-gen-compound was *designed* for T790M) |
| (d) Free-text in `Target.name` | "EGFR (T790M-selective)" | No modelling work |

**Recommendation for v1: (d), defer structure to v1.1.** Justification:
- You can add `variant` and `mutation_codes` later without data loss — backfill by parsing paren-expressions in existing `Target.name` strings.
- The NLP query feature can ship without variant structure; queries like *"EGFR L858R selective compounds"* can still route via name matching for now.
- Premature structure on variants risks modelling the wrong dimensions (selectivity-vs-activity, mutation-vs-isoform, single-vs-compound-mutant) before we've seen enough real cases.

The `Target.parent` field (already in the schema) gives us a natural path to v1.1 structure: a "EGFR-T790M" target can have `parent="EGFR"`, letting the wild-type target act as the family root without duplicating gene linkage.

## 10. Consumers

### 10.1 NLP query resolver (primary)

See `NLP_QUERY_PROPOSAL.md` §6.1. Once this proposal lands, the NLP resolver's target-matching pool expands from `{Target.name}` to:

```
normalize(Target.name)
  ∪ {normalize(g.symbol)    for g in Target.genes.all()}
  ∪ {normalize(a)            for g in Target.genes.all() for a in g.aliases}
  ∪ {normalize(g.name)       for g in Target.genes.all()}
```

Same Levenshtein-1-on-≥4-chars rule applies across the combined pool. Ambiguity resolution (UI picker) is unchanged.

### 10.2 Future consumers

- **External ELN / collaboration deposition**: UniProt accession and HGNC ID are the lingua franca for identifying targets in external systems.
- **Cross-target SAR search**: "compounds active against any BCL-2-family member" becomes `Gene.objects.filter(symbol__in=['BCL2','BCLXL','MCL1']).targets`.
- **Literature / database sync**: pulling bioactivity data from ChEMBL or PubChem requires UniProt accession mapping.
- **Dashboards**: gene-level rollups across multiple targets in the same family.

## 11. Locked decisions

| # | Decision | Rationale |
|---|----------|-----------|
| 1 | Gene is a separate model, linked M2M to Target | Shared hydration, cross-target queries by gene, portability across DBs |
| 2 | HGNC is the canonical source for human gene metadata | Authoritative, free, no auth, industry-standard |
| 3 | UniProt is an optional second enrichment layer, not required for v1 | HGNC already provides enough; avoid unnecessary complexity |
| 4 | Aliases live on `Gene` as JSONField, not a separate table | Read-mostly flat list with no metadata worth modelling |
| 5 | Annotator workflow: type gene symbol(s); system hydrates everything else | Minimises annotator burden; keeps target creation ergonomic |
| 6 | `Target.genes` may be empty; NLP matching falls back to name-only | Handles legacy rows and non-HGNC targets gracefully |
| 7 | Variant/mutant modelling deferred to v1.1; free-text in `Target.name` for now | Avoid premature structure; no data loss when upgraded |
| 8 | Live HGNC API in v1; local snapshot available as a management command for stricter deployments | Good default + escape hatch; no forced offline mode |
| 9 | Backfill is opt-in (UI banner), not forced | No data migration required; maintainers fill in as they touch each target |

## 12. Open questions

1. **Non-human targets.** If Akane's / another group ever runs compounds against a mouse or zebrafish target, HGNC is the wrong source (MGI / ZFIN respectively). Do we add a `Gene.species` field now, or wait? Lean: add now, default "human"; cheap.
2. **Gene isoforms.** AR has clinically important splice variants (AR-V7). Does an isoform-selective compound programme get its own target, or does isoform metadata attach to the Target? Probably isoform-as-variant (v1.1 work); flagging.
3. **Relationship to `Target.parent`.** If we add variant structure later, does `parent` become reserved for variant-of relationships, or does it continue to hold arbitrary hierarchies (project > sub-project > …)? Worth a one-time audit of current parent usage before making assumptions.
4. **Who curates `Gene.hydration_status = 'manual'`?** If an annotator overrides HGNC, should there be a UI audit trail? Probably yes, using the existing `modified_by` pattern — but explicit.
5. **Rehydration cadence.** 90 days is a guess. HGNC symbols do change (ABCA12 was once ICR2A, etc.) but slowly. Is there a user-visible consequence to a stale alias? If not, longer cadence is fine.
6. **Collision handling.** Two Target rows legitimately link to the same Gene (e.g. EGFR-WT programme and EGFR-T790M programme). NLP matching for "EGFR" would then return both, triggering the clarification picker. Expected, probably fine — flagging so it doesn't surprise.

## 13. Rollout / migration

- **Migration 1**: add `Gene` model; add `Target.genes` M2M. Pure additive, no data change.
- **Migration 2**: no data migration required — existing `Target` rows have `genes=[]` by default.
- **Management command 1**: `hydrate_genes` — nightly refresh of stale/pending Gene rows.
- **Management command 2**: `import_hgnc_snapshot` — optional bulk import from the HGNC JSON snapshot file.
- **Feature flag**: `TARGET_GENE_LINKAGE_ENABLED` — off by default while developing; on once annotator UI is ready.
- **Rollout order**:
  1. Ship models + migrations + `hydrate_genes` command.
  2. Add Target detail page "no genes recorded" banner and inline editor.
  3. Update Target creation form to require gene symbol(s) for new targets.
  4. Wire NLP resolver to use the expanded matching pool (gated on the feature flag).
  5. Encourage backfill of existing targets via the UI banner.

No database migration is destructive; no existing functionality breaks at any step.
