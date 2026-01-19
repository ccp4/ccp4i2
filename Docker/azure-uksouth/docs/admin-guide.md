# Platform Administrator Guide

Reference documentation for platform administrators.

## Data Management

The platform manages data across multiple apps with automated backups and migration tools.

### Data Apps

| App | Description |
|-----|-------------|
| **CCP4i2** | Projects, jobs, and files for crystallographic computing |
| **Compound Registry** | Compound registration, batches, suppliers, and targets |
| **Assays** | Assay experiments, protocols, and analysis results |
| **Constructs** | Plasmid and construct database with GenBank files |

### Fixture Naming Convention

All fixtures use Django's JSON format with timestamp prefixes:

```
YYYYMMDD-HH-MM-{app}.json
```

| App | Fixture Name |
|-----|--------------|
| CCP4i2 | `*-CCP4i2.json` |
| Registry | `*-RegisterCompounds.json` |
| Assays | `*-AssayCompounds.json` |
| Constructs | `*-ConstructDatabase.json` |
| Users | `*-auth.json` |

---

## Automated Backups

Database backups run automatically via Azure Container Apps Job.

| Setting | Value |
|---------|-------|
| **Schedule** | Daily at 2:00 AM UTC |
| **Retention** | 30 days (older backups automatically deleted) |
| **Location** | Azure Files: `/mnt/azure-files/fixtures/` |
| **Format** | Django JSON fixtures |

### Manual Backup (CLI)

```bash
# Trigger backup manually
./scripts/deploy-backup-job.sh trigger

# View backup job status
./scripts/deploy-backup-job.sh logs

# Check job status
./scripts/deploy-backup-job.sh status
```

---

## Data Migration

Tools for migrating data from legacy systems to the new platform.

### Fixture Import (Web UI)

Use the import sections in the Admin panel to upload legacy fixture files directly.
For files larger than 100MB, use the CLI tools below.

### Fixture Migration (CLI)

```bash
# List fixtures in legacy storage
./scripts/migrate-fixtures.sh list

# Copy latest fixtures to new storage
./scripts/migrate-fixtures.sh latest

# Copy all historical fixtures
./scripts/migrate-fixtures.sh copy
```

### Media Files Migration

```bash
# Check what will be migrated
./scripts/migrate-media.sh check

# Perform the migration
./scripts/migrate-media.sh copy

# Verify migrated files
./scripts/migrate-media.sh list
```

**Path transformations:**
- Most files: `media/*` → `media/*` (unchanged)
- Batch QC: `media/RegBatchQCFile_NCL-*` → `media/RegisterCompounds/BatchQCFiles/NCL-*`

---

## Import Commands Reference

Management commands for importing legacy data via CLI.

### Accessing the Container App CLI

To run management commands on the Azure deployment, open a shell session on the server container:

```bash
# Using Azure CLI
az containerapp exec \
  --name ccp4i2-bicep-server \
  --resource-group ccp4i2-bicep-rg-uksouth \
  --command /bin/bash

# Once connected, activate the environment and navigate to the app
source /mnt/ccp4data/ccp4-20251105/bin/ccp4.setup-sh
cd /app
```

Alternatively, use the Azure Portal:
1. Navigate to **Container Apps** → **ccp4i2-bicep-server**
2. Select **Console** from the left menu
3. Choose `/bin/bash` as the startup command
4. Run `source /mnt/ccp4data/ccp4-20251105/bin/ccp4.setup-sh && cd /app`

### Import CCP4i2 Data

```bash
python manage.py import_legacy_ccp4i2 \
  --fixture ccp4i2.json \
  --remap-dirs /old/path:/new/path
```

### Import Compound Registry and Assays

```bash
python manage.py import_legacy_compounds \
  --auth-fixture auth.json \
  --registry-fixture RegisterCompounds.json \
  --assays-fixture AssayCompounds.json
```

### Import Construct Database

```bash
python manage.py import_legacy_constructs \
  --auth-fixture auth.json \
  --constructs-fixture ConstructDatabase.json
```

### Import Users

```bash
python manage.py import_legacy_users auth.json
```

### Common Options

| Option | Description |
|--------|-------------|
| `--dry-run` | Validate without saving changes |
| `--verbose` | Show detailed progress |
| `--output-dir DIR` | Save transformed fixtures for inspection |

---

## Troubleshooting

### Fixture Import Errors

**Problem**: Import fails with "Invalid JSON" error

**Solution**: Legacy fixtures may contain Django debug output before the JSON array.
The import commands handle this automatically by finding the first `[` character.
If manual inspection is needed:

```bash
# Check first few characters of fixture
head -c 100 fixture.json

# If there's text before the JSON, extract it
python -c "
import json
with open('fixture.json') as f:
    content = f.read()
    start = content.find('[')
    data = json.loads(content[start:])
    print(f'Found {len(data)} records')
"
```

### Large File Imports

For fixtures larger than 100MB, the web UI upload may time out.
Use the CLI tools instead:

1. Copy fixtures to Azure Files storage:
   ```bash
   ./scripts/migrate-fixtures.sh latest
   ```

2. Connect to the container and run import from mounted storage:
   ```bash
   # Connect to the container
   az containerapp exec \
     --name ccp4i2-bicep-server \
     --resource-group ccp4i2-bicep-rg-uksouth \
     --command /bin/bash

   # Set up environment and run import
   source /mnt/ccp4data/ccp4-20251105/bin/ccp4.setup-sh
   cd /app
   python manage.py import_legacy_compounds \
     --registry-fixture /mnt/azure-files/fixtures/RegisterCompounds.json
   ```
