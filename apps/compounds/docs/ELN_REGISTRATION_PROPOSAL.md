# ELN-linked Compound Registration

**Status**: Phase 1 implemented (April 2026)
**Scope**: Phase 1 of supporting Kawamura-group-style registration. Independent of [HELM_REGISTRATION_PROPOSAL.md](HELM_REGISTRATION_PROPOSAL.md), which covers Phase 2 (formal monomer composition / structure derivation).

## Motivation

Some groups (currently the Kawamura group, generally any group using an Electronic Lab Notebook) do not record compound provenance as a paper-notebook number + page. Instead, they record a hyperlink to an ELN page (OneNote, in Akane's case). Their working registry is an Excel spreadsheet whose "ELN reference" cell carries a display label (e.g. `KF001`, `KF022`) and a hyperlink target (a OneNote sharepoint.com deep-link).

Concretely, a single ELN page describes a *synthesis campaign* that can produce many compounds. For example KF001 ("SA157 Cyclic peptide synthesis and characterisation") yielded 6 peptide scaffolds (ARD-CP1…CP6) in up to 3 cyclisation variants each (linear / DBX / disulfide), so up to ~18 distinct registered compounds reference the same ELN page. KF022 in turn takes products of KF001 and couples a FAM linker to them.

The current `Compound.labbook_number` + `Compound.page_number` fields cannot represent this:

1. There is no "page number" — the ELN entry *is* the unit, with a URL.
2. The relationship is many compounds → one ELN entry, so the reference is naturally a foreign key to a separate row, not a pair of integers per compound.
3. Other linked artefacts (Chemdraw files, QC PDFs) follow the same pattern in their workflow — an external URL rather than an upload.

## Data model

### New: `LabNotebookEntry`

Represents one ELN page or one paper notebook page. Many compounds may reference the same entry.

```python
class LabNotebookEntry(models.Model):
    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    supplier = models.ForeignKey(Supplier, on_delete=models.PROTECT, related_name='notebook_entries')
    sequence_number = models.IntegerField(
        help_text="Per-supplier sequential index (e.g. 1 for KF001, 22 for KF022)"
    )
    title = models.CharField(max_length=255, blank=True,
        help_text='Page title, e.g. "SA157 Cyclic peptide synthesis and characterisation"')
    date = models.DateField(null=True, blank=True)

    # ELN form
    url = models.URLField(max_length=2048, blank=True,
        help_text="Hyperlink to ELN page (OneNote, LabArchives, etc.)")

    # Paper-notebook form (mirrors current Compound.labbook_number/page_number)
    labbook_number = models.IntegerField(null=True, blank=True)
    page_number = models.IntegerField(null=True, blank=True)

    created_at = models.DateTimeField(auto_now_add=True)
    created_by = models.ForeignKey(User, on_delete=models.SET_NULL, null=True, blank=True,
        related_name='created_notebook_entries')

    class Meta:
        unique_together = [('supplier', 'sequence_number')]
        ordering = ['supplier', 'sequence_number']

    def __str__(self):
        return self.label

    @cached_property
    def label(self):
        """E.g. 'KF001', 'KF022'. Width is fixed at 3 digits."""
        if self.supplier and self.supplier.initials:
            return f'{self.supplier.initials}{self.sequence_number:03d}'
        return f'#{self.sequence_number:03d}'

    def clean(self):
        # Either ELN url or paper labbook+page should be set; allow neither (entry created
        # ahead of detail) but warn if both forms set.
        if self.url and (self.labbook_number or self.page_number):
            raise ValidationError(
                "Set either ELN url or paper labbook_number+page_number, not both")
```

### Changes to `Compound`

```python
notebook_entry = models.ForeignKey(
    'LabNotebookEntry', on_delete=models.PROTECT,
    null=True, blank=True, related_name='compounds')

helm_notation = models.CharField(
    max_length=4096, blank=True,
    help_text="Optional HELM string. Phase 1: opaque text, no validation.")

sequence_display = models.CharField(
    max_length=512, blank=True,
    help_text="Human-readable sequence with modifications shown inline, e.g. "
              "'FAM-linker-MCDWDIYRFPNHHC(1,4-Xylene)-NH2'. Free text.")

# Deprecated, kept for back-compat during migration:
# labbook_number, page_number, compound_number — see migration plan below
```

`Compound.barcode` becomes:

```python
@cached_property
def barcode(self):
    if self.notebook_entry:
        return self.notebook_entry.label  # e.g. 'KF001'
    if self.supplier and self.supplier.initials and self.labbook_number is not None:
        return f'{self.supplier.initials}-{self.labbook_number}-{self.page_number}'
    return None
```

### New: `CompoundDocument`

Per-compound document attachment for structure files (Chemdraw). Exactly one of `url` or `file` is set.

```python
class CompoundDocument(models.Model):
    KIND_CHOICES = [
        ('chemdraw', 'Chemdraw structure'),
        ('other',    'Other'),
    ]

    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    compound = models.ForeignKey(Compound, on_delete=models.CASCADE, related_name='documents')
    kind = models.CharField(max_length=16, choices=KIND_CHOICES)
    label = models.CharField(max_length=255, blank=True,
        help_text="Display text, e.g. 'CP1.cdxml'")

    url = models.URLField(max_length=2048, blank=True)
    file = models.FileField(upload_to=_compound_document_path, max_length=255, blank=True)

    created_at = models.DateTimeField(auto_now_add=True)
    created_by = models.ForeignKey(User, on_delete=models.SET_NULL, null=True, blank=True)

    def clean(self):
        if bool(self.url) == bool(self.file):
            raise ValidationError("Set exactly one of url or file")
```

### Extended: `BatchQCFile`

Extended to support URLs in addition to file uploads (for linking to QC documents in SharePoint/ELN without downloading them). QC documents (LCMS, HRMS, NMR, HPLC) are attached to batches, not compounds, because they pertain to physical samples.

```python
class BatchQCFile(models.Model):
    KIND_CHOICES = [
        ('lcms', 'LC-MS'),
        ('hrms', 'HR-MS'),
        ('nmr', 'NMR'),
        ('hplc', 'HPLC'),
        ('other', 'Other'),
    ]

    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    batch = models.ForeignKey(Batch, on_delete=models.CASCADE, related_name='qc_files')
    kind = models.CharField(max_length=16, choices=KIND_CHOICES, default='other')
    label = models.CharField(max_length=255, blank=True,
        help_text="Display text, e.g. 'LCMS-001'")

    # Either url or file should be set, not both
    url = models.URLField(max_length=2048, blank=True,
        help_text="External URL to QC document (SharePoint, ELN, cloud storage)")
    file = models.FileField(upload_to=_batch_qc_path, max_length=255, blank=True)
    original_filename = models.CharField(max_length=255, blank=True, null=True)
    comments = models.TextField(blank=True, null=True)
    uploaded_at = models.DateTimeField(auto_now_add=True, null=True, blank=True)
```

### Document Type Distinction

**Important design decision**: Structure files (Chemdraw) attach to **compounds** (the abstract entity), while QC/analysis files attach to **batches** (physical samples).

| Document Type | Model | Rationale |
|---------------|-------|-----------|
| Chemdraw structure | `CompoundDocument` | Structure is independent of physical batch |
| LCMS, HRMS | `BatchQCFile` | Analysis of a specific physical sample |
| NMR, HPLC | `BatchQCFile` | Analysis of a specific physical sample |
| Purity reports | `BatchQCFile` | Batch-specific quality assessment |

## Migration

1. **Add new tables and columns** (`LabNotebookEntry`, `CompoundDocument`, `Compound.notebook_entry`, `Compound.helm_notation`, `Compound.sequence_display`).
2. **Data migration**: for each `Compound` with `(supplier_id, labbook_number, page_number)` set, find-or-create a `LabNotebookEntry(supplier=..., sequence_number=labbook_number, labbook_number=labbook_number, page_number=page_number)` and set `compound.notebook_entry`. The `sequence_number` choice is a heuristic; for paper-notebook records with no obvious sequence, use `labbook_number * 1000 + page_number` to keep uniqueness.
3. **Keep legacy fields** (`Compound.labbook_number`, `Compound.page_number`, `Batch.labbook_number`, `Batch.page_number`, `Assay.labbook_number`, `Assay.page_number`) for one release cycle, marked deprecated in their `help_text`. Remove in a follow-up migration once UIs are clean.
4. **Backfill suppliers** for any Compounds whose labbook fields were set but supplier was null — flag those in the migration log; do not invent a supplier.

## Spreadsheet importer

Akane's spreadsheet (see [Sample: ARd selection_Peptide]) has these columns: NCL_ID, Batch, Peptide ID, Target, **ELN reference (hyperlinked)**, **Chemdraw file (hyperlinked)**, Chemical Formula, Exact mass, Sequence (and downstream columns).

`apps/compounds/registry/management/commands/import_legacy_compounds.py` is the existing importer for the old format. Add a sibling command `import_eln_spreadsheet.py` (or extend with a `--format=eln` flag) that:

1. Reads the workbook with `openpyxl` in `data_only=False` so hyperlinks are accessible via `cell.hyperlink.target`.
2. For each row, extracts:
   - `eln_label` from the cell display value of column **ELN reference** (e.g. `KF022`)
   - `eln_url` from `cell.hyperlink.target`
   - Same pair for the **Chemdraw file** column.
3. Parses `eln_label` against a regex `^([A-Z]{1,4})(\d+)$` → `(initials, sequence_number)`. Looks up `Supplier` by `initials`; creates if unknown (with a warning in the import log — do not silently auto-create in production runs, gate behind `--create-suppliers`).
4. Upserts `LabNotebookEntry` keyed on `(supplier, sequence_number)`. First-seen-wins for `url`; if a later row has a different URL for the same label, log a warning.
5. Creates the `Compound` with `notebook_entry` set, `sequence_display` from the sequence column, `formula`, `exact_mass`, etc.
6. If a Chemdraw URL is present, creates `CompoundDocument(kind='chemdraw', url=..., label=display_text)`.
7. Leaves `helm_notation` blank.

Importer should be idempotent: re-running with the same workbook updates existing rows by `(supplier, supplier_ref)` or `(notebook_entry, peptide_id)` rather than duplicating.

## Frontend

### Bulk Import UI

The bulk import page (`app/registry/import/page.tsx`) includes a **Document Categorization** section that:

1. **Auto-detects hyperlinked columns** from the spreadsheet using ExcelJS hyperlink extraction.
2. **Heuristic-based pre-assignment**:
   - Columns containing "chemdraw", "cdx", "structure" → Compound Documents
   - Columns containing "lcms", "hrms", "nmr", "hplc", "qc", "purity" → Batch QC Documents
   - Other hyperlinked columns → Default to Compound Documents (user can reclassify)
3. **Two multi-select autocompletes** allow users to move columns between categories.
4. **Auto-enable batch creation** when any Batch QC columns are assigned (batch is required to attach QC files).

### Components affected

| File | Change |
|---|---|
| [components/compounds/SpreadsheetUpload.tsx](../frontend/components/compounds/SpreadsheetUpload.tsx) | Extract hyperlinks from Excel cells using ExcelJS `cell.hyperlink` property. |
| [app/registry/new/page.tsx](../frontend/app/registry/new/page.tsx) | Single compound registration form: ELN Reference fields (label + URL), HELM notation, sequence display, legacy paper notebook fields. |
| [app/registry/import/page.tsx](../frontend/app/registry/import/page.tsx) | Bulk import: Document categorization UI with two autocomplete multi-selects. Auto-enable batch creation when QC columns present. |
| [app/registry/compounds/[id]/page.tsx](../frontend/app/registry/compounds/[id]/page.tsx) | Compound detail view: ELN/Sequence accordion section showing notebook entry link, sequence, HELM, and documents. |

### TypeScript types

New types in [types/compounds/models.ts](../frontend/types/compounds/models.ts):
- `LabNotebookEntry`, `LabNotebookEntryCompact`
- `CompoundDocument`, `CompoundDocumentKind`
- `BatchQCFileKind` (lcms, hrms, nmr, hplc, other)

### REST endpoints

- `LabNotebookEntryViewSet` (`/api/compounds/notebook-entries/`) - list, retrieve, create, update; filter by supplier, search by label
- `CompoundDocumentViewSet` (`/api/compounds/compound-documents/`) - CRUD for compound documents; filter by compound, kind

## Out of scope (Phase 2 / future)

- **HELM parsing or validation**. `helm_notation` is opaque text in Phase 1.
- **Monomer library**. Comes with HELM in [HELM_REGISTRATION_PROPOSAL.md](HELM_REGISTRATION_PROPOSAL.md).
- **OneNote scraping via Microsoft Graph**. The OneNote page has rich per-compound data (theoretical mass, observed mass, yield, success/fail) that the spreadsheet flattens. A future `CompoundMeasurement` model could capture this; out of scope here.
- **Provenance graph between ELN entries** (e.g. KF022 modifies products of KF001). Out of scope.
- **Removal of legacy labbook_number/page_number fields**. Deferred to a follow-up migration.

## Resolved Questions

1. **Should `LabNotebookEntry.title` be auto-extracted from the OneNote page on first import?**
   → Hand-entered, optional. Phase 2 OneNote scraper can backfill.

2. **Width of the supplier sequence — `KF001` (3 digits) is the convention but the existing Akane spreadsheet has a typo (`KF0022` should be `KF022`).**
   → Importer normalises to 3 digits and warns on mismatch.

3. **Should `CompoundDocument` be per-batch instead of per-compound for QC?**
   → **Resolved**: Structure files (Chemdraw) go to `CompoundDocument` (per-compound). QC files (LCMS, HRMS, NMR, HPLC) go to `BatchQCFile` (per-batch). The frontend bulk import UI allows users to categorize hyperlinked columns into either category, with heuristic-based defaults. Selecting any batch QC columns forces batch creation.
