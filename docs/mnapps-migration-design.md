# Compounds App Migration Design Document

**Status:** Planning → Implementation
**Last Updated:** 2025-01-05
**Authors:** Martin Noble, Claude

## Overview

Migration of AssayCompounds and RegisterCompounds Django apps from the legacy MNApps codebase into CCP4i2, with:
- Clean data model naming
- Proper audit trail via django-reversion
- New React/Next.js frontend (replacing Django templates)
- Flexible deployment (Azure, desktop, or standalone)

## Architecture Decisions

### 1. Directory Structure

**Decision:** Place in `apps/compounds/` at repository root for deployment flexibility.

```
ccp4i2/
├── apps/                           # Optional/pluggable applications
│   └── compounds/                  # Compound registration & assay management
│       ├── __init__.py
│       ├── settings.py             # Extends core or azure settings
│       ├── urls.py                 # Combined API routes
│       ├── registry/               # Compound registration (was RegisterCompounds)
│       │   ├── __init__.py
│       │   ├── apps.py
│       │   ├── models.py           # Target, Compound, Batch, Supplier
│       │   ├── serializers.py
│       │   ├── views.py            # DRF ViewSets
│       │   └── migrations/
│       └── assays/                 # Assay data (was AssayCompounds)
│           ├── __init__.py
│           ├── apps.py
│           ├── models.py           # Protocol, Assay, DataSeries, Analysis
│           ├── serializers.py
│           ├── views.py
│           └── migrations/
│
├── server/                         # Core CCP4i2 Django (unchanged)
├── client/                         # Core CCP4i2 React (unchanged)
└── Docker/azure/
    └── azure_extensions/           # Azure-specific features only
```

**Why `apps/compounds/` instead of `Docker/azure/`:**
- Not tied to Docker or Azure in naming
- Supports desktop/Electron deployment
- Clear "optional apps" semantic
- Easy to add other optional apps later

**Settings chain (Azure deployment):**
```
compounds.settings
  → imports azure_extensions.settings
    → imports ccp4i2.config.settings
```

**Settings chain (Desktop/standalone):**
```
compounds.settings
  → imports ccp4i2.config.settings
```

**Local development:**
```bash
# With Compounds app
export PYTHONPATH="$PWD:$PWD/apps"
export DJANGO_SETTINGS_MODULE=compounds.settings
ccp4-python server/manage.py runserver

# Without Compounds app (core CCP4i2 only)
export DJANGO_SETTINGS_MODULE=ccp4i2.config.settings
ccp4-python server/manage.py runserver
```

**Azure Container deployment:**
```dockerfile
COPY apps/ /app/apps/
ENV PYTHONPATH="/app:/app/apps"
ENV DJANGO_SETTINGS_MODULE="compounds.settings"
```

**Desktop/Electron deployment:**
```javascript
process.env.PYTHONPATH = `${appPath}:${appPath}/apps`;
process.env.DJANGO_SETTINGS_MODULE = "compounds.settings";
```

### 2. The "Project" Naming Collision

**Problem:** "Project" means different things:
- CCP4i2: A crystallographic experiment folder
- RegisterCompounds: A drug discovery target/campaign

**Decision:** Rename in UI and code (not just UI):

| Old Model | New Model | Django App | UI Term | Route |
|-----------|-----------|------------|---------|-------|
| `RegProjects` | `Target` | `compounds.registry` | Target/Campaign | `/api/targets/` |
| `RegData` | `Compound` | `compounds.registry` | Compound | `/api/compounds/` |
| `RegBatch` | `Batch` | `compounds.registry` | Batch | `/api/batches/` |
| `RegSuppliers` | `Supplier` | `compounds.registry` | Supplier | `/api/suppliers/` |
| `RegBatchQCFile` | `BatchQCFile` | `compounds.registry` | QC File | nested under batch |
| `RegDataTemplate` | `CompoundTemplate` | `compounds.registry` | Template | `/api/templates/` |
| `Protocol` | `Protocol` | `compounds.assays` | Protocol | `/api/protocols/` |
| `Experiment` | `Assay` | `compounds.assays` | Assay | `/api/assays/` |
| `DilutionSeries` | `DilutionSeries` | `compounds.assays` | Dilution Series | `/api/dilutions/` |
| `DataSeries` | `DataSeries` | `compounds.assays` | Data Series | nested under assay |
| `Analysis` | `AnalysisResult` | `compounds.assays` | Result | nested under series |
| `Hypothesis` | `Hypothesis` | `compounds.assays` | Hypothesis | `/api/hypotheses/` |
| CCP4i2 `Project` | (unchanged) | `ccp4i2.db` | Project | `/api/projects/` |

### 3. User Authentication & JIT Provisioning

**Decision:** Use Django's built-in `auth.User` with Just-In-Time provisioning.

**Flow:**
```
Azure AD / AWS Cognito / Other
            │
            ▼
    Provider Auth Middleware
    (validates token, extracts claims)
    Sets: request.auth_email, request.auth_claims
            │
            ▼
    JIT User Provisioning Middleware
    User.objects.get_or_create(email=auth_email)
    Sets: request.user = django_user
            │
            ▼
    Django Views / DRF ViewSets
    request.user → Django User instance
    reversion.set_user(request.user) ✓
```

**Key:** Email address is the stable identifier across providers.

**Middleware order in settings.py:**
```python
MIDDLEWARE = [
    # ... security middleware ...
    "ccp4i2.middleware.azure_auth.AzureADAuthMiddleware",
    "ccp4i2.middleware.jit_user.JITUserProvisioningMiddleware",
    # ... other middleware ...
]
```

**Future providers:** Each creates a middleware that validates its tokens and sets `request.auth_email` and `request.auth_claims`. The JIT middleware handles user creation uniformly.

### 4. Audit Trail with django-reversion

**Decision:** Use django-reversion with explicit context managers in DRF ViewSets.

```python
import reversion

class CompoundViewSet(viewsets.ModelViewSet):
    queryset = Compound.objects.all()
    serializer_class = CompoundSerializer

    def perform_create(self, serializer):
        with reversion.create_revision():
            instance = serializer.save()
            reversion.set_user(self.request.user)
            reversion.set_comment("Created via API")

    def perform_update(self, serializer):
        with reversion.create_revision():
            instance = serializer.save()
            reversion.set_user(self.request.user)
            reversion.set_comment("Updated via API")

    @action(detail=True, methods=['get'])
    def history(self, request, pk=None):
        """Expose version history via API."""
        obj = self.get_object()
        versions = Version.objects.get_for_object(obj)
        return Response([{
            'id': v.id,
            'date': v.revision.date_created,
            'user': v.revision.user.email if v.revision.user else None,
            'comment': v.revision.comment,
        } for v in versions])
```

---

## Data Model Cleanup

### Naming Inconsistencies to Fix

#### 1. ForeignKey `_id` suffix problem
Django automatically adds `_id` to FK field names in the database. Current code uses `project_id = ForeignKey(...)` which becomes `project_id_id` in DB.

| Current | Fixed |
|---------|-------|
| `project_id = ForeignKey(RegProjects)` | `target = ForeignKey(Target)` |
| `regdata_id = ForeignKey(RegData)` | `compound = ForeignKey(Compound)` |
| `batch_id = ForeignKey(RegBatch)` | `batch = ForeignKey(CompoundBatch)` |

#### 2. Remove "Reg" prefix
| Current | Fixed |
|---------|-------|
| `RegProjects` | `Target` |
| `RegData` | `Compound` |
| `RegBatch` | `CompoundBatch` |
| `RegSuppliers` | `Supplier` |
| `RegBatchQCFile` | `BatchQCFile` |
| `RegDataTemplate` | `CompoundTemplate` |

#### 3. Consistent UUID primary keys
Use `id` as the field name for UUID PKs (Django convention):

| Current | Fixed |
|---------|-------|
| `project_uuid = UUIDField(primary_key=True)` | `id = UUIDField(primary_key=True)` |
| `reg_uuid = UUIDField(primary_key=True)` | `id = UUIDField(primary_key=True)` |
| etc. | etc. |

#### 4. snake_case consistency
| Current (camelCase) | Fixed (snake_case) |
|---------------------|-------------------|
| `compoundName` | `compound_name` |
| `startColumn` | `start_column` |
| `endColumn` | `end_column` |
| `skipPoints` | `skip_points` |
| `extractedData` | `extracted_data` |
| `dilutionSeries` | `dilution_series` |
| `analysisMethod` | `analysis_method` |
| `preferredDilutions` | `preferred_dilutions` |
| `svgFile` | `svg_file` |
| `dataFile` | `data_file` |
| `rdkitSmiles` | `rdkit_smiles` |

#### 5. Remove legacy db_column mappings
The `RegProjects` model has `db_column="PROJECT_NAME"` etc. from Oracle migration. Drop these in new schema.

---

## Proposed New Models

### register_compounds/models.py

```python
import uuid
from django.db import models
from django.contrib.auth import get_user_model

User = get_user_model()


class Supplier(models.Model):
    """Chemical supplier / synthesis source."""
    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    name = models.CharField(max_length=64, unique=True)
    initials = models.CharField(max_length=16, unique=True, blank=True)

    class Meta:
        ordering = ['name']

    def __str__(self):
        return self.name


class Target(models.Model):
    """
    Drug discovery target/campaign.

    Note: This was previously called "RegProjects". Renamed to avoid
    confusion with CCP4i2's crystallographic "Project" concept.
    """
    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    name = models.CharField(max_length=255, unique=True)
    parent = models.ForeignKey(
        'self', on_delete=models.SET_NULL, null=True, blank=True,
        related_name='children'
    )
    created_at = models.DateTimeField(auto_now_add=True)

    class Meta:
        ordering = ['name']

    def __str__(self):
        return self.name


class Compound(models.Model):
    """
    Registered compound with unique NCL-XXXXXXXX identifier.
    """
    STEREO_CHOICES = [
        ('unset', 'Unset'),
        ('achiral', 'Achiral'),
        ('racemic', 'Racemic mixture'),
        ('single_unknown', 'Single enantiomer, configuration unknown'),
        ('r_enantiomer', 'R enantiomer'),
        ('s_enantiomer', 'S enantiomer'),
        # ... etc
    ]

    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    reg_number = models.IntegerField(unique=True, editable=False)  # Auto-assigned
    target = models.ForeignKey(Target, on_delete=models.PROTECT, related_name='compounds')

    # Chemistry
    smiles = models.TextField()
    rdkit_smiles = models.TextField(blank=True)  # Canonical form
    inchi = models.TextField(blank=True)
    molecular_weight = models.FloatField(null=True, blank=True)
    stereo_comment = models.CharField(max_length=70, choices=STEREO_CHOICES, default='unset')

    # Provenance
    supplier = models.ForeignKey(Supplier, on_delete=models.SET_NULL, null=True)
    supplier_ref = models.CharField(max_length=255, blank=True)
    labbook_number = models.IntegerField()
    page_number = models.IntegerField()
    compound_number = models.IntegerField(default=1)

    # Metadata
    registered_by = models.ForeignKey(User, on_delete=models.SET_NULL, null=True)
    registered_at = models.DateTimeField(auto_now_add=True)
    modified_at = models.DateTimeField(auto_now=True)
    comments = models.TextField(blank=True)

    # Generated files
    svg_file = models.ImageField(upload_to='compounds/svg/', blank=True)

    class Meta:
        ordering = ['-reg_number']

    def __str__(self):
        return self.formatted_id

    @property
    def formatted_id(self):
        return f'NCL-{self.reg_number:08d}'

    def save(self, *args, **kwargs):
        if not self.reg_number:
            from django.db.models import Max
            max_reg = Compound.objects.aggregate(Max('reg_number'))['reg_number__max']
            self.reg_number = (max_reg or 25999) + 1
        super().save(*args, **kwargs)


class CompoundBatch(models.Model):
    """A synthesis batch of a compound."""
    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    compound = models.ForeignKey(Compound, on_delete=models.CASCADE, related_name='batches')
    batch_number = models.IntegerField()

    # Provenance
    supplier = models.ForeignKey(Supplier, on_delete=models.SET_NULL, null=True)
    supplier_ref = models.CharField(max_length=255, blank=True)
    labbook_number = models.IntegerField()
    page_number = models.IntegerField()

    # Properties
    amount = models.DecimalField(max_digits=10, decimal_places=4, null=True)
    salt_code = models.CharField(max_length=128, blank=True)
    molecular_weight = models.FloatField(null=True, blank=True)

    # Metadata
    registered_at = models.DateTimeField(auto_now_add=True)
    comments = models.TextField(blank=True)

    class Meta:
        ordering = ['compound', 'batch_number']
        unique_together = [['compound', 'batch_number']]

    def __str__(self):
        return f'{self.compound.formatted_id}/{self.batch_number}'


class BatchQCFile(models.Model):
    """QC documentation for a batch (NMR, LCMS, etc.)."""
    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    batch = models.ForeignKey(CompoundBatch, on_delete=models.CASCADE, related_name='qc_files')
    file = models.FileField(upload_to='batches/qc/')
    comments = models.TextField(blank=True)
    uploaded_at = models.DateTimeField(auto_now_add=True)

    def __str__(self):
        return f'QC for {self.batch}'
```

### assay_compounds/models.py

```python
import uuid
from django.db import models
from django.contrib.auth import get_user_model
from mnapps.register_compounds.models import Target, Compound

User = get_user_model()


class DilutionSeries(models.Model):
    """Standard concentration series for dose-response experiments."""
    UNIT_CHOICES = [
        ('nM', 'nanomolar'),
        ('uM', 'micromolar'),
        ('mM', 'millimolar'),
    ]

    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    concentrations = models.JSONField(default=list)  # [100, 30, 10, 3, 1, 0.3, ...]
    unit = models.CharField(max_length=10, choices=UNIT_CHOICES, default='nM')

    def __str__(self):
        return f"{','.join(str(c) for c in self.concentrations)} {self.unit}"


class AssayProtocol(models.Model):
    """
    Assay protocol definition.

    Renamed from 'Protocol' to avoid ambiguity.
    """
    ANALYSIS_METHOD_CHOICES = [
        ('hill_langmuir', 'Hill-Langmuir'),
        ('hill_langmuir_fix_hill', 'Hill-Langmuir (fixed Hill)'),
        ('hill_langmuir_fix_minmax', 'Hill-Langmuir (fixed min/max)'),
        ('ms_intact', 'MS-Intact'),
        ('table_of_values', 'Table of values'),
    ]

    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    name = models.CharField(max_length=256)
    analysis_method = models.CharField(max_length=50, choices=ANALYSIS_METHOD_CHOICES)
    preferred_dilutions = models.ForeignKey(
        DilutionSeries, on_delete=models.SET_NULL, null=True, blank=True
    )
    pherastar_table = models.CharField(max_length=32, blank=True)

    created_by = models.ForeignKey(User, on_delete=models.SET_NULL, null=True)
    created_at = models.DateTimeField(auto_now_add=True)
    comments = models.TextField(blank=True)

    class Meta:
        ordering = ['name']

    def __str__(self):
        return self.name


class Assay(models.Model):
    """
    An assay experiment (plate read, etc.).

    Renamed from 'Experiment' for clarity.
    """
    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    protocol = models.ForeignKey(AssayProtocol, on_delete=models.PROTECT, related_name='assays')
    target = models.ForeignKey(
        Target, on_delete=models.SET_NULL, null=True, blank=True,
        related_name='assays',
        help_text="Target this assay is testing against (may differ from compound's registered target)"
    )

    data_file = models.FileField(upload_to='assays/data/')
    labbook_number = models.IntegerField(null=True, blank=True)
    page_number = models.IntegerField(null=True, blank=True)

    created_by = models.ForeignKey(User, on_delete=models.SET_NULL, null=True)
    created_at = models.DateTimeField(auto_now_add=True)
    comments = models.TextField(blank=True)

    class Meta:
        ordering = ['-created_at']
        verbose_name_plural = 'assays'

    def __str__(self):
        return f'{self.protocol.name} - {self.created_at:%Y-%m-%d}'


class DataSeries(models.Model):
    """A single compound's dose-response data within an assay."""
    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    assay = models.ForeignKey(Assay, on_delete=models.CASCADE, related_name='data_series')
    compound = models.ForeignKey(
        Compound, on_delete=models.SET_NULL, null=True, blank=True,
        related_name='assay_data'
    )
    compound_name = models.CharField(max_length=256, blank=True)  # Original name from file

    # Position in source file
    row = models.IntegerField()
    start_column = models.IntegerField()
    end_column = models.IntegerField()

    # Data
    dilution_series = models.ForeignKey(DilutionSeries, on_delete=models.SET_NULL, null=True)
    extracted_data = models.JSONField(default=dict)
    skip_points = models.JSONField(default=list)  # Indices to exclude from fit

    # Results
    analysis = models.OneToOneField(
        'AnalysisResult', on_delete=models.SET_NULL, null=True,
        related_name='data_series'
    )
    svg_file = models.ImageField(upload_to='assays/plots/', blank=True)

    class Meta:
        ordering = ['compound_name']
        verbose_name_plural = 'data series'

    def __str__(self):
        return f'{self.compound_name} in {self.assay}'


class AnalysisResult(models.Model):
    """Fitted results for a data series."""
    STATUS_CHOICES = [
        ('valid', 'Valid'),
        ('invalid', 'Invalid'),
        ('unassigned', 'Unassigned'),
    ]

    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    status = models.CharField(max_length=20, choices=STATUS_CHOICES, default='unassigned')
    results = models.JSONField(default=dict)  # {EC50: ..., Hill: ..., KPI: ...}

    def __str__(self):
        if self.status != 'valid':
            return self.status.upper()
        kpi = self.results.get('KPI')
        if kpi and kpi in self.results:
            return f"{self.results[kpi]:.2f}"
        return 'Valid'


class Hypothesis(models.Model):
    """Compound design hypothesis for a target."""
    STATUS_CHOICES = [
        ('pending', 'Pending'),
        ('rejected', 'Rejected'),
        ('chemistry', 'In Chemistry'),
        ('shelved', 'Shelved'),
        ('made', 'Made'),
    ]

    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    target = models.ForeignKey(Target, on_delete=models.CASCADE, related_name='hypotheses')
    parent = models.ForeignKey('self', on_delete=models.SET_NULL, null=True, blank=True)

    smiles = models.CharField(max_length=1024)
    rationale = models.TextField()
    model_url = models.URLField(blank=True)

    status = models.CharField(max_length=20, choices=STATUS_CHOICES, default='pending')
    product_compound = models.ForeignKey(
        Compound, on_delete=models.SET_NULL, null=True, blank=True,
        help_text="Link to registered compound if hypothesis was synthesized"
    )
    completion_notes = models.TextField(blank=True)

    svg_file = models.ImageField(upload_to='hypotheses/svg/', blank=True)
    created_at = models.DateTimeField(auto_now_add=True)
    updated_at = models.DateTimeField(auto_now=True)

    class Meta:
        verbose_name_plural = 'hypotheses'
        ordering = ['-created_at']
```

---

## Data Migration Strategy

### Outstanding Questions

**Q1: Same database or separate?**
Will MNApps use the same PostgreSQL instance as CCP4i2?

- **Same DB** (recommended): Simpler deployment, can have FK relationships to auth.User
- **Separate DB**: More isolation, but complicates user references

**Q2: Table naming approach?**

- **(a) New tables, migrate data**: Create `mnapps_target`, `mnapps_compound`, etc. Write migration scripts to copy data from old tables.
- **(b) Keep old table names**: Use `db_table = 'RegisterCompounds_regprojects'` to map new models to existing tables. Less disruption but perpetuates old naming.

**Q3: Reversion history migration?**

The `reversion_version` and `reversion_revision` tables reference content types by old model names.

- **(a) Migrate history**: Update content_type references. Complex but preserves full audit trail.
- **(b) Fresh start**: New models get new reversion tracking. Old history remains in old tables (read-only reference).
- **(c) Hybrid**: Keep old tables, create views to unify history queries.

**Q4: User references in historical data?**

Old records have `user_name = CharField(...)` strings. New records will have `ForeignKey(User)`.

- Keep both fields during transition?
- Write migration to match usernames to User records where possible?

---

## Frontend Architecture

### Route Structure

```
/                           # CCP4i2 home (unchanged)
/project/                   # CCP4i2 projects (unchanged)
/project/{id}/              # CCP4i2 project detail (unchanged)

# MNApps routes (feature-flagged)
/targets/                   # Target list
/targets/{id}/              # Target detail + compounds
/compounds/                 # Compound search/list
/compounds/{id}/            # Compound detail + batches
/compounds/{id}/batches/    # Batch management
/assays/                    # Assay list
/assays/{id}/               # Assay detail + data series
/protocols/                 # Protocol management
/hypotheses/                # Hypothesis tracker
```

### Feature Flag

```typescript
// client/renderer/lib/features.ts
export const FEATURES = {
  MNAPPS_ENABLED: process.env.NEXT_PUBLIC_MNAPPS_ENABLED === 'true',
}
```

```tsx
// In navigation
{FEATURES.MNAPPS_ENABLED && (
  <>
    <NavItem href="/targets">Targets</NavItem>
    <NavItem href="/compounds">Compounds</NavItem>
    <NavItem href="/assays">Assays</NavItem>
  </>
)}
```

### API Hooks Pattern (useSWR)

```typescript
// client/renderer/hooks/useCompounds.ts
import useSWR from 'swr';

interface Compound {
  id: string;
  formatted_id: string;
  smiles: string;
  target: { id: string; name: string };
  // ...
}

export function useCompound(id: string) {
  return useSWR<Compound>(`/api/compounds/${id}/`);
}

export function useCompounds(targetId?: string) {
  const params = targetId ? `?target=${targetId}` : '';
  return useSWR<Compound[]>(`/api/compounds/${params}`);
}

export function useCompoundHistory(id: string) {
  return useSWR(`/api/compounds/${id}/history/`);
}
```

---

## Implementation Phases

### Phase 1: Infrastructure
- [ ] Create `Docker/azure/mnapps/` directory structure
- [ ] Create `mnapps/settings.py` extending azure_extensions
- [ ] Create JIT user provisioning middleware
- [ ] Add django-reversion to requirements

### Phase 2: Models & API
- [ ] Create cleaned-up models (register_compounds, assay_compounds)
- [ ] Create DRF serializers
- [ ] Create ViewSets with reversion integration
- [ ] Create URL routing
- [ ] Write and test migrations

### Phase 3: Data Migration
- [ ] Decide on migration strategy (Q1-Q4 above)
- [ ] Write data migration scripts
- [ ] Test with copy of production data
- [ ] Plan rollback strategy

### Phase 4: Frontend
- [ ] Add feature flag infrastructure
- [ ] Create TypeScript types from API
- [ ] Build Target/Compound/Assay pages
- [ ] Add molecule visualization (RDKit or similar)
- [ ] History/audit trail UI

### Phase 5: Deployment
- [ ] Update Dockerfile to include mnapps
- [ ] Update Azure deployment scripts
- [ ] Test end-to-end in staging
- [ ] Production deployment

---

## Dependencies

### Python (add to requirements.txt for Azure)
```
django-reversion>=5.0
rdkit  # For molecule handling - check Azure container compatibility
```

### JavaScript (client)
```
rdkit-js  # or similar for client-side molecule rendering
swr       # Already in use
```

---

## Resolved Questions

| # | Question | Decision |
|---|----------|----------|
| Q1 | Same database or separate? | **Same DB** - MNApps uses same PostgreSQL instance as CCP4i2 |
| Q2 | Table naming | **(a) New clean names** - create fresh tables, migrate data |
| Q3 | Reversion history | **(c) Hybrid** - keep old tables read-only, new tracking for new models |
| Q5 | Local dev auth | **Auto dev-user** - auto-create test user in DEBUG mode |
| Q6 | RDKit availability | **Yes** - available via ccp4-python, no container changes needed |

---

## Q4: Historical `user_name` String Migration

### Current State (from code analysis)

**RegisterCompounds:**
- `RegData.user_name` = CharField, but it's a **user-editable form field**, not auto-populated from `request.user`
- Users typed their own names (inconsistently, or may have left blank)
- The admin has a bug: `obj.username = request.user.username` but the field is `user_name`
- This field is **NOT a reliable auth link** - just a text annotation

**AssayCompounds:**
- Uses `created_by = ForeignKey(User)` - properly links to Django auth
- `RevisionMixin` is used on views, so reversion captured the user

**Key insight:** The reversion tables (`reversion_revision`) likely have more reliable user data than the model's `user_name` field, since `RevisionMixin` captures `request.user` automatically.

### The Challenge

Historical `user_name` strings are unreliable. Going forward, JIT provisioning creates Users keyed by **email**. We need to bridge these where possible, but accept that some historical data may not link cleanly.

### Options

**(a) Direct FK lookup during migration**
```python
for old_compound in OldRegData.objects.all():
    user = User.objects.filter(username=old_compound.user_name).first()
    new_compound.registered_by = user  # May be None
```
- **Risk:** Users may not exist yet (pre-JIT), losing attribution

**(b) Pre-migration: Update user_name to emails**
```sql
-- Before migration, replace usernames with emails where possible
UPDATE RegisterCompounds_regdata r
SET user_name = u.email
FROM auth_user u
WHERE r.user_name = u.username;
```
Then migrate with email-based lookup.

**(c) Keep legacy field**
```python
class Compound(models.Model):
    registered_by = models.ForeignKey(User, null=True, ...)
    legacy_registered_by = models.CharField(blank=True)  # Preserve original
```

### Recommendation: Hybrid (b) + (c)

1. **Pre-migration phase:**
   - Export mapping of `username → email` from existing auth_user table
   - Update `user_name` field to email where matches exist
   - Log unmatched usernames for manual review

2. **During migration:**
   - Create `registered_by` FK via email lookup
   - Store original `user_name` value in `legacy_registered_by` field

3. **Post-migration:**
   - For unmatched records: display "Registered by: martin (legacy)" in UI
   - When those users first login via JIT, optionally prompt to claim historical records

4. **Going forward:**
   - JIT provisioning ensures User always exists before any API writes
   - `registered_by` FK is always populated for new records

### Migration Script Sketch

```python
def migrate_compounds():
    # Load email mapping from old auth_user
    email_map = {u.username: u.email for u in OldUser.objects.all()}

    for old in OldRegData.objects.all():
        # Try to find/create user by email
        email = email_map.get(old.user_name) or old.user_name  # May already be email
        user = None
        if email and '@' in email:
            user, _ = User.objects.get_or_create(
                email=email.lower(),
                defaults={'username': email.lower(), 'is_active': True}
            )

        Compound.objects.create(
            # ... other fields ...
            registered_by=user,
            legacy_registered_by=old.user_name or '',
        )
```

This approach:
- Preserves all historical attribution (nothing lost)
- Creates User records for known emails (enabling reversion tracking)
- Allows graceful handling of unmatched legacy usernames

### Simplified Recommendation (given code analysis)

Since `user_name` is unreliable, the simplest approach is:

1. **Keep `legacy_registered_by = CharField`** - copy the original `user_name` value as-is
2. **Set `registered_by = None`** for all migrated records (don't try to link)
3. **Going forward**: JIT provisioning ensures proper `registered_by` FK for new records

**Why this works:**
- The old `user_name` was just a text field users typed into
- Reversion already has the real user info via `RevisionMixin`
- Trying to match unreliable strings to emails adds complexity with little benefit
- Clean break: new system does it right, old data preserved for reference

**If you want to attempt matching:**
- Join via reversion: Look up the user from `reversion_revision.user_id` for each object's first revision
- This gives the actual authenticated user, not the form-field string
