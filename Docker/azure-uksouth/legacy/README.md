# Legacy Scripts

These scripts are from earlier deployment approaches and are retained for reference only.

**Do not use these scripts for deployment.** Use the scripts in `../scripts/` instead.

## Contents

| File | Original Purpose | Replaced By |
|------|------------------|-------------|
| `deploy.sh` | All-in-one build/push/deploy script | `scripts/build-and-push.sh`, `scripts/deploy-applications.sh` |
| `azcopy.sh` | Generic Azure blob copy utility | `scripts/upload-ccp4data.sh` |
| `azcopy-files.sh` | Azure Files upload utility | `scripts/upload-ccp4data.sh` |
| `container-app.bicep` | Single-file Bicep template | `infrastructure/` directory (modular Bicep) |
| `parameters.example.json` | ARM parameters template | `.env.deployment` (environment-based config) |

## Current Deployment Approach

The current deployment uses:

1. **Modular Bicep templates** in `infrastructure/`
2. **Specialized scripts** in `scripts/` for each operation
3. **Environment-based configuration** via `.env.deployment`
4. **Layered Docker images** for efficient builds

See `../OPERATIONS.md` for current deployment procedures.
