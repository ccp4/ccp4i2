# i2remote - Remote CLI for CCP4i2

`i2remote` is a command-line interface for interacting with a remote CCP4i2 server via HTTP API. It mirrors the syntax of the local `i2` CLI, making it easy to switch between local and remote workflows.

## Features

- **Familiar syntax**: Same commands as local `i2` CLI
- **CCP4 data awareness**: Uses `ccp4-python` for file validation (MTZ, PDB, mmCIF)
- **Azure AD authentication**: Secure login with browser or device code flow
- **Token management**: Automatic token storage and refresh

## Requirements

### Python Dependencies

```bash
# Required
pip install requests

# Optional (for login command)
pip install msal
```

### CCP4 Environment

For full functionality including file validation, run with `ccp4-python`:

```bash
source /path/to/ccp4/bin/ccp4.setup-sh
```

## Installation

The CLI is located at `cli/i2remote.py`. You can run it directly or create an alias:

```bash
# Direct execution
ccp4-python cli/i2remote.py projects

# Or create an alias
alias i2remote='ccp4-python /path/to/ccp4i2/cli/i2remote.py'
```

## Configuration

### Setting the Server URL

```bash
# Via environment variable
export CCP4I2_API_URL=https://your-server.azurecontainerapps.io/api/proxy

# Or via config command
i2remote config set api_url https://your-server.azurecontainerapps.io/api/proxy
```

### Configuration File

Settings are stored in `~/.ccp4i2remote.json` with secure permissions (0600).

View current configuration:
```bash
i2remote config
```

## Authentication

The remote server may require Azure AD authentication. There are several ways to authenticate:

### Method 1: Browser Login (Recommended)

Opens your default browser for interactive authentication:

```bash
# First, configure Azure AD settings
i2remote config set azure_client_id <your-app-client-id>
i2remote config set azure_tenant_id <your-tenant-id>

# Then login
i2remote login
```

### Method 2: Device Code Flow (SSH/Headless)

For environments without a browser (SSH sessions, servers):

```bash
i2remote login --device-code
```

This displays a code and URL:
```
============================================================
To sign in, use a web browser to open the page
https://microsoft.com/devicelogin and enter the code ABC123XYZ
============================================================
```

### Method 3: Azure CLI

If you have Azure CLI installed:

```bash
# Login to Azure
az login

# Get token and set it
export CCP4I2_API_TOKEN=$(az account get-access-token \
  --resource <client-id> --query accessToken -o tsv)

# Use i2remote
i2remote projects
```

### Method 4: Manual Token

If you have a token from another source:

```bash
i2remote config set api_token <your-jwt-token>
```

### Checking Authentication Status

```bash
i2remote whoami
```

Output:
```
Current authentication:
  Name:     John Smith
  Email:    john.smith@example.com
  Subject:  abc123def456...
  Expires:  2024-12-29T15:30:00+00:00 (in 0:45:23)

Config file: /home/user/.ccp4i2remote.json
```

### Logging Out

```bash
i2remote logout
```

## Commands

### Projects

```bash
# List all projects
i2remote projects

# Create a new project
i2remote projects create myproject

# Show project details
i2remote projects show myproject

# Show project job tree
i2remote projects tree myproject
```

### Jobs

```bash
# List jobs in a project
i2remote jobs myproject

# Create a new job
i2remote jobs create myproject refmac

# Run a job
i2remote jobs run myproject 1

# Check job status
i2remote jobs status myproject 1

# Wait for job completion
i2remote jobs wait myproject 1

# Show job KPIs (R-factors, etc.)
i2remote jobs kpi myproject 1

# Clone a job
i2remote jobs clone myproject 1

# Validate output files (requires ccp4-python)
i2remote jobs validate myproject 1

# Set/get job parameters
i2remote jobs set-param myproject 1 XYZIN fullPath=/path/to/file.pdb
i2remote jobs get-param myproject 1 XYZIN

# Upload files as job parameters
i2remote jobs upload-param myproject 1 XYZIN /path/to/model.pdb
i2remote jobs upload-param myproject 1 F_SIGF /path/to/data.mtz '/*/*/[FP,SIGFP]'
```

### Setting Job Parameters with FileUse References

The `set-param` command supports referencing output files from previous jobs using fileUse syntax. This allows you to chain jobs together by passing outputs from one job as inputs to another.

#### FileUse Syntax

```
[jobIndex].jobParamName[paramIndex]         # Reference any job
task_name[jobIndex].jobParamName[paramIndex] # Reference specific task type
```

- **jobIndex**: Job index (negative for counting from end, e.g., `-1` = most recent)
- **jobParamName**: Output parameter name (e.g., `XYZOUT`, `HKLOUT`, `F_SIGF_OUT`)
- **paramIndex**: (Optional) Index when multiple files have the same param name
- **task_name**: (Optional) Filter to jobs of a specific task type

#### Examples

```bash
# Use output from most recent job
i2remote jobs set-param myproject 5 XYZIN '[-1].XYZOUT[0]'

# Use output from most recent refmac job
i2remote jobs set-param myproject 5 XYZIN 'refmac[-1].XYZOUT'

# Use output from second-to-last aimless job
i2remote jobs set-param myproject 5 HKLIN 'aimless[-2].HKLOUT'

# Explicit fileUse= prefix (equivalent to auto-detected patterns)
i2remote jobs set-param myproject 5 XYZIN 'fileUse=[-1].XYZOUT[0]'
```

#### How it Works

1. The CLI detects fileUse patterns automatically (or via explicit `fileUse=` prefix)
2. Queries the server to resolve the reference to actual file metadata
3. Sets the job parameter with the resolved file path

This enables scripting complex crystallographic workflows:

```bash
# Run a refinement pipeline
i2remote jobs create myproject refmac
i2remote jobs set-param myproject 2 XYZIN 'phaser[-1].XYZOUT'
i2remote jobs set-param myproject 2 HKLIN 'phaser[-1].HKLOUT'
i2remote jobs run myproject 2
```

### Uploading Files as Job Parameters

The `upload-param` command uploads a local file to the server and sets it as a job parameter. This is useful when you need to provide input files that don't already exist on the server.

```bash
# Upload a PDB file as the XYZIN parameter
i2remote jobs upload-param myproject 5 XYZIN /local/path/to/model.pdb

# Upload an MTZ file with column selection
i2remote jobs upload-param myproject 5 F_SIGF /local/path/to/data.mtz '/*/*/[FP,SIGFP]'
```

#### MTZ Column Selectors

For reflection data (MTZ files), you can optionally specify which columns to extract:

- `/*/*/[FP,SIGFP]` - Extract F and SIGF columns from any crystal/dataset
- `/*/*/[I,SIGI]` - Extract I and SIGI columns (intensity data)
- `/crystal1/dataset1/[F,SIGF]` - Specify exact crystal and dataset names

If no column selector is provided for general MTZ containers, the file is stored as-is. For specific reflection types (like `F_SIGF`), appropriate columns are auto-detected.

### Files

```bash
# List files in a job
i2remote files myproject 1

# Display file contents
i2remote files cat myproject 1 output.log
```

### Reports

```bash
# Get job report XML
i2remote report myproject 1
```

### Export/Import

```bash
# Export a job
i2remote export job myproject 1

# Export a project
i2remote export project myproject

# Import a project
i2remote import project_backup.zip
```

## Environment Variables

| Variable | Description |
|----------|-------------|
| `CCP4I2_API_URL` | Base URL for the CCP4i2 API |
| `CCP4I2_API_TOKEN` | Authentication token (overrides saved token) |
| `AZURE_CLIENT_ID` | Azure AD application client ID |
| `AZURE_TENANT_ID` | Azure AD tenant ID |

## Examples

### Complete Workflow

```bash
# Configure server
i2remote config set api_url https://ccp4i2.azurecontainerapps.io/api/proxy
i2remote config set azure_client_id 12345678-abcd-1234-abcd-123456789012
i2remote config set azure_tenant_id 87654321-dcba-4321-dcba-987654321098

# Authenticate
i2remote login

# Create project and job
i2remote projects create myexperiment
i2remote jobs create myexperiment aimless

# Run and monitor
i2remote jobs run myexperiment 1
i2remote jobs wait myexperiment 1

# Check results
i2remote jobs status myexperiment 1
i2remote jobs kpi myexperiment 1
i2remote jobs validate myexperiment 1
```

### Scripted Usage

```bash
#!/bin/bash
# Example: Monitor all running jobs

PROJECT="myproject"

# List all jobs and their status
i2remote jobs $PROJECT | while read line; do
    echo "$line"
done

# Wait for specific job
i2remote jobs wait $PROJECT 5
if [ $? -eq 0 ]; then
    echo "Job completed successfully"
    i2remote jobs kpi $PROJECT 5
fi
```

### Using with Azure CLI in CI/CD

```bash
#!/bin/bash
# Example: GitHub Actions or Azure DevOps script

# Login using service principal (set in CI environment)
az login --service-principal \
  -u $AZURE_CLIENT_ID \
  -p $AZURE_CLIENT_SECRET \
  --tenant $AZURE_TENANT_ID

# Get token
export CCP4I2_API_TOKEN=$(az account get-access-token \
  --resource $AZURE_CLIENT_ID --query accessToken -o tsv)

export CCP4I2_API_URL=https://ccp4i2.azurecontainerapps.io/api/proxy

# Run commands
i2remote projects
i2remote jobs myproject
```

## Troubleshooting

### "API URL not configured"

Set the server URL:
```bash
i2remote config set api_url https://your-server/api/proxy
```

### "Azure Client ID not configured"

For login functionality, configure Azure AD:
```bash
i2remote config set azure_client_id <your-client-id>
i2remote config set azure_tenant_id <your-tenant-id>
```

### "MSAL library required for login"

Install the MSAL package:
```bash
pip install msal
```

### "Authentication failed" or 401 errors

1. Check if token is expired: `i2remote whoami`
2. Re-authenticate: `i2remote login`
3. Verify Azure AD configuration matches the server

### "requests package required"

Install the requests package:
```bash
pip install requests
```

### File validation shows "(gemmi not available)"

Run with ccp4-python for full CCP4 library access:
```bash
source /path/to/ccp4/bin/ccp4.setup-sh
ccp4-python cli/i2remote.py jobs validate myproject 1
```

## Security Notes

- Tokens are stored in `~/.ccp4i2remote.json` with 0600 permissions
- Tokens are not displayed in `config` output (shown as `****`)
- Use `logout` to clear credentials when done
- For CI/CD, prefer environment variables over stored tokens
