# CCP4i2 Command Line Interface

Modern CLI for the cdata-codegen project (CCP4i2 replacement).

## Table of Contents

- [Installation](#installation)
- [Quick Start](#quick-start)
- [Command Structure](#command-structure)
- [Commands Reference](#commands-reference)
  - [Projects](#projects)
  - [Jobs](#jobs)
  - [Files](#files)
  - [Plugins](#plugins)
  - [Shortcuts](#shortcuts)
- [Global Flags](#global-flags)
- [Examples](#examples)
- [Tab Completion](#tab-completion)

## Installation

### Unix/Linux/macOS

1. Make the script executable:
```bash
chmod +x ccp4i2
```

2. Add to PATH (choose one):

**Option A: Symlink to /usr/local/bin**
```bash
sudo ln -s /path/to/cdata-codegen/ccp4i2 /usr/local/bin/ccp4i2
```

**Option B: Add to your PATH in ~/.bashrc or ~/.zshrc**
```bash
export PATH="/path/to/cdata-codegen:$PATH"
```

### Windows

**Option A: Using .bat file**
Add the directory containing `ccp4i2.bat` to your PATH environment variable.

**Option B: Using PowerShell**
1. Ensure PowerShell execution policy allows scripts:
```powershell
Set-ExecutionPolicy -ExecutionPolicy RemoteSigned -Scope CurrentUser
```

2. Add to PATH or create an alias in your PowerShell profile:
```powershell
Set-Alias ccp4i2 "C:\path\to\cdata-codegen\ccp4i2.ps1"
```

## Quick Start

```bash
# List all projects
ccp4i2 projects list

# Create a new project
ccp4i2 projects create my_project --description "My crystallography project"

# List jobs in a project
ccp4i2 jobs list my_project

# Create and run a job
ccp4i2 jobs create my_project import_merged_mtz --title "Import MTZ file"
ccp4i2 jobs run my_project 1

# Or create and run in one step
ccp4i2 run my_project import_merged_mtz
```

## Command Structure

The CLI follows a resource-oriented design:

```
ccp4i2 <resource> <action> [<identifier>] [options]
```

- **Resource**: The entity type (projects, jobs, files, plugins)
- **Action**: What to do (list, create, show, run, etc.)
- **Identifier**: Name, UUID, or number to identify a specific entity
- **Options**: Flags and arguments to modify behavior

### Smart Identifier Resolution

The CLI automatically detects identifier types:

- **Names**: `toxd`, `my_project`
- **UUIDs**: `a1b2c3d4-e5f6-7890-abcd-ef1234567890` (full or short prefix like `a1b2c3d4`)
- **Numbers**: `1`, `42` (for jobs)

You don't need to specify which type - the CLI figures it out!

## Commands Reference

### Projects

Manage CCP4i2 projects.

#### `projects list` (aliases: `ls`)

List all projects with filtering and sorting.

```bash
ccp4i2 projects list [options]

Options:
  --filter PATTERN      Filter by name pattern (supports wildcards like 'toxd*')
  --sort-by FIELD       Sort by: name, created, accessed, jobs (default: created)
  --reverse             Reverse sort order
  --limit N             Limit number of results
  --json                Output as JSON
  -q, --quiet           Output only project names
```

Examples:
```bash
# List all projects
ccp4i2 projects list

# Find projects matching a pattern
ccp4i2 projects list --filter "toxd*"

# Show recently accessed projects
ccp4i2 projects list --sort-by accessed --reverse --limit 5

# Get JSON output for scripting
ccp4i2 projects list --json
```

#### `projects create` (aliases: `new`)

Create a new project with directory structure.

```bash
ccp4i2 projects create <name> [options]

Arguments:
  name                  Project name (alphanumeric, underscore, hyphen only)

Options:
  -d, --description TEXT   Project description
  --directory PATH         Custom project directory (optional)
  --json                   Output as JSON
  -q, --quiet              Output only project UUID
```

Examples:
```bash
# Simple project
ccp4i2 projects create my_project

# With description
ccp4i2 projects create toxd --description "Toxd structure determination"

# Custom directory
ccp4i2 projects create test --directory /custom/path/test
```

#### `projects show` (aliases: `info`, `get`)

Show detailed information about a project.

```bash
ccp4i2 projects show <project> [options]

Arguments:
  project               Project name or UUID

Options:
  --json                Output as JSON
```

Examples:
```bash
ccp4i2 projects show toxd
ccp4i2 projects show a1b2c3d4  # Using UUID prefix
```

#### `projects export`

Export a project to a ZIP archive.

```bash
ccp4i2 projects export <project> -o <output> [options]

Arguments:
  project               Project name or UUID

Options:
  -o, --output PATH     Output ZIP file path (required)
  --jobs NUMBERS        Comma-separated job numbers to export (e.g., "1,3,5")
  --detach              Run export in background
```

Examples:
```bash
# Export entire project
ccp4i2 projects export toxd -o toxd_backup.zip

# Export specific jobs only
ccp4i2 projects export toxd -o toxd_subset.zip --jobs "1,2,5"

# Run in background
ccp4i2 projects export large_project -o export.zip --detach
```

#### `projects import`

Import a project from a ZIP archive.

```bash
ccp4i2 projects import <zip_file>

Arguments:
  zip_file              ZIP file to import
```

Examples:
```bash
ccp4i2 projects import toxd_backup.zip
```

#### `projects tree`

Show a beautiful directory tree structure for a project.

```bash
ccp4i2 projects tree <project> [options]

Arguments:
  project               Project name or UUID

Options:
  --depth N             Maximum depth to display (default: 3, 0 = unlimited)
  --no-sizes            Do not show file sizes
  --show-hidden         Show hidden files and directories
```

Examples:
```bash
# Show project directory structure
ccp4i2 projects tree toxd

# Show deeper tree
ccp4i2 projects tree toxd --depth 5

# Show all files (unlimited depth)
ccp4i2 projects tree toxd --depth 0

# Hide file sizes
ccp4i2 projects tree toxd --no-sizes
```

#### `projects cat` (aliases: `view`, `show-file`)

View file contents within a project directory.

```bash
ccp4i2 projects cat <project> <file_path> [options]

Arguments:
  project               Project name or UUID
  file_path             Relative path to file within project directory

Options:
  --lines N             Show first N lines (like head -n)
  --tail N              Show last N lines (like tail -n)
  --binary              Allow viewing binary files (shows hex dump)
  --no-line-numbers     Do not show line numbers
```

Examples:
```bash
# View a file in the project
ccp4i2 projects cat toxd CCP4_JOBS/job_1/input_params.xml

# Show first 20 lines
ccp4i2 projects cat toxd CCP4_JOBS/job_1/program.log --lines 20

# Show last 50 lines (like tail)
ccp4i2 projects cat toxd CCP4_JOBS/job_1/program.log --tail 50

# View without line numbers
ccp4i2 projects cat toxd CCP4_JOBS/job_1/output.txt --no-line-numbers

# View binary file (hex dump)
ccp4i2 projects cat toxd CCP4_JOBS/job_1/output.mtz --binary
```

### Jobs

Manage jobs within projects.

#### `jobs list` (aliases: `ls`)

List jobs in a project with filtering and sorting.

```bash
ccp4i2 jobs list <project> [options]

Arguments:
  project               Project name or UUID

Options:
  --status STATUS       Filter by status (PENDING, RUNNING, FINISHED, FAILED, etc.)
  --task TASK           Filter by task name
  --filter PATTERN      Filter by title pattern
  --sort-by FIELD       Sort by: number, title, task, status, created, finished
  --reverse             Reverse sort order
  --limit N             Limit number of results
  --json                Output as JSON
  -q, --quiet           Output only job numbers
```

Examples:
```bash
# List all jobs
ccp4i2 jobs list toxd

# Show only finished jobs
ccp4i2 jobs list toxd --status FINISHED

# Find jobs by title pattern
ccp4i2 jobs list toxd --filter "refine*"

# Show most recent jobs
ccp4i2 jobs list toxd --sort-by created --reverse --limit 10
```

#### `jobs create` (aliases: `new`)

Create a new job in a project.

```bash
ccp4i2 jobs create <project> <task> [options]

Arguments:
  project               Project name or UUID
  task                  Task/plugin name (e.g., import_merged_mtz, refmac5)

Options:
  --from-job JOB        Copy inputs from this job (number or UUID)
  --title TEXT          Job title
```

Examples:
```bash
# Create a simple job
ccp4i2 jobs create toxd import_merged_mtz

# With custom title
ccp4i2 jobs create toxd refmac5 --title "Refinement round 1"

# Copy inputs from another job
ccp4i2 jobs create toxd refmac5 --from-job 3
```

#### `jobs show` (aliases: `info`, `get`)

Show detailed information about a job.

```bash
ccp4i2 jobs show <project> <job>

Arguments:
  project               Project name or UUID
  job                   Job number or UUID
```

Examples:
```bash
ccp4i2 jobs show toxd 1
ccp4i2 jobs show toxd a1b2c3d4  # Using job UUID
```

#### `jobs run` (aliases: `start`, `execute`)

Run a job.

```bash
ccp4i2 jobs run <project> <job> [options]

Arguments:
  project               Project name or UUID
  job                   Job number or UUID

Options:
  --detach              Run in background (detached process)
```

Examples:
```bash
# Run in foreground
ccp4i2 jobs run toxd 1

# Run in background
ccp4i2 jobs run toxd 1 --detach
```

#### `jobs validate` (aliases: `check`)

Validate job parameters before running.

```bash
ccp4i2 jobs validate <project> <job> [options]

Arguments:
  project               Project name or UUID
  job                   Job number or UUID

Options:
  -o, --output PATH     Save validation report to file (XML format)
  --json                Output as JSON
```

Examples:
```bash
# Quick validation check
ccp4i2 jobs validate toxd 1

# Save detailed report
ccp4i2 jobs validate toxd 1 -o validation_report.xml

# Get JSON output
ccp4i2 jobs validate toxd 1 --json
```

#### `jobs export`

Export a job to a ZIP archive.

```bash
ccp4i2 jobs export <project> <job> -o <output>

Arguments:
  project               Project name or UUID
  job                   Job number or UUID

Options:
  -o, --output PATH     Output ZIP file path (required)
  --json                Output as JSON
```

Examples:
```bash
ccp4i2 jobs export toxd 1 -o job_1_backup.zip
```

#### `jobs tree`

Show a beautiful directory tree structure for a job.

```bash
ccp4i2 jobs tree <project> <job> [options]

Arguments:
  project               Project name or UUID
  job                   Job number or UUID

Options:
  --depth N             Maximum depth to display (default: 2, 0 = unlimited)
  --no-sizes            Do not show file sizes
  --show-hidden         Show hidden files and directories
```

Examples:
```bash
# Show job directory structure
ccp4i2 jobs tree toxd 1

# Show deeper tree
ccp4i2 jobs tree toxd 1 --depth 4

# Show all files (unlimited depth)
ccp4i2 jobs tree toxd 1 --depth 0

# Hide file sizes
ccp4i2 jobs tree toxd 1 --no-sizes
```

#### `jobs cat` (aliases: `view`, `show-file`)

View file contents within a job directory.

```bash
ccp4i2 jobs cat <project> <job> <file_path> [options]

Arguments:
  project               Project name or UUID
  job                   Job number or UUID
  file_path             Relative path to file within job directory

Options:
  --lines N             Show first N lines (like head -n)
  --tail N              Show last N lines (like tail -n)
  --binary              Allow viewing binary files (shows hex dump)
  --no-line-numbers     Do not show line numbers
```

Examples:
```bash
# View a job file (relative to job directory)
ccp4i2 jobs cat toxd 1 input_params.xml

# Show first 30 lines of a log
ccp4i2 jobs cat toxd 1 program.log --lines 30

# Show last 100 lines (like tail -n 100)
ccp4i2 jobs cat toxd 1 program.log --tail 100

# View without line numbers
ccp4i2 jobs cat toxd 1 output.txt --no-line-numbers

# View binary file (hex dump)
ccp4i2 jobs cat toxd 1 output.mtz --binary
```

#### `jobs set-status`

Set the status of a job.

```bash
ccp4i2 jobs set-status <project> <job> <status>

Arguments:
  project               Project name or UUID
  job                   Job number or UUID
  status                New status value
```

Examples:
```bash
ccp4i2 jobs set-status toxd 1 FINISHED
```

#### `jobs set-param`

Set a job parameter value.

```bash
ccp4i2 jobs set-param <project> <job> <param> <value>

Arguments:
  project               Project name or UUID
  job                   Job number or UUID
  param                 Parameter path
  value                 Parameter value
```

Examples:
```bash
ccp4i2 jobs set-param toxd 1 "cycles" "10"
```

#### `jobs report`

Get a job report.

```bash
ccp4i2 jobs report <project> <job> [options]

Arguments:
  project               Project name or UUID
  job                   Job number or UUID

Options:
  --format FORMAT       Output format: text, json, xml (default: text)
```

Examples:
```bash
ccp4i2 jobs report toxd 1
ccp4i2 jobs report toxd 1 --format json
```

### Files

File operations (preview, list).

#### `files preview` (aliases: `show`, `cat`)

Preview file contents.

```bash
ccp4i2 files preview <path> [options]

Arguments:
  path                  File path

Options:
  --lines N             Number of lines to show (default: 50)
```

Examples:
```bash
ccp4i2 files preview /path/to/logfile.log
ccp4i2 files preview output.mtz --lines 100
```

### Plugins

Plugin/task discovery and information.

#### `plugins list` (aliases: `ls`)

List available plugins/tasks.

```bash
ccp4i2 plugins list [options]

Options:
  --filter PATTERN      Filter by name pattern
```

#### `plugins show` (aliases: `info`)

Show detailed information about a plugin.

```bash
ccp4i2 plugins show <plugin>

Arguments:
  plugin                Plugin name
```

### Shortcuts

Convenience commands for common workflows.

#### `run`

Create and run a job in one step.

```bash
ccp4i2 run <project> <task> [options]

Arguments:
  project               Project name or UUID
  task                  Task/plugin name

Options:
  --from-job JOB        Copy inputs from this job
  --detach              Run in background
  --title TEXT          Job title
```

Examples:
```bash
# Quick create and run
ccp4i2 run toxd import_merged_mtz

# With context job
ccp4i2 run toxd refmac5 --from-job 1 --detach
```

#### `status` (aliases: `st`)

Quick status overview of projects or a specific project.

```bash
ccp4i2 status [project]

Arguments:
  project               Project name or UUID (optional)
```

Examples:
```bash
# Show all projects
ccp4i2 status

# Show specific project
ccp4i2 status toxd
```

## Global Flags

These flags work with all commands:

- `--json` - Output as JSON (machine-readable)
- `-q, --quiet` - Minimal output (IDs/names only)
- `-v, --verbose` - Extra debug information
- `--no-color` - Disable colored output
- `-h, --help` - Show help for any command

## Examples

### Complete Workflow Example

```bash
# 1. Create a new project
ccp4i2 projects create toxd_structure --description "Toxd refinement project"

# 2. List all projects to confirm
ccp4i2 projects list

# 3. Create an import job
ccp4i2 jobs create toxd_structure import_merged_mtz --title "Import initial MTZ"

# 4. Validate the job
ccp4i2 jobs validate toxd_structure 1

# 5. Run the job
ccp4i2 jobs run toxd_structure 1 --detach

# 6. Check job status
ccp4i2 jobs list toxd_structure --status RUNNING

# 7. Create a refinement job using the import as context
ccp4i2 jobs create toxd_structure refmac5 --from-job 1 --title "Initial refinement"

# 8. Export the project
ccp4i2 projects export toxd_structure -o backup.zip
```

### Scripting with JSON Output

```bash
# Get all finished jobs as JSON
ccp4i2 jobs list toxd --status FINISHED --json > finished_jobs.json

# Get project list for processing
projects=$(ccp4i2 projects list --json)
echo "$projects" | jq -r '.[].name'

# Create project and capture UUID
project_uuid=$(ccp4i2 projects create test --quiet)
echo "Created project: $project_uuid"
```

### Using Filters and Patterns

```bash
# Find all projects starting with "toxd"
ccp4i2 projects list --filter "toxd*"

# Find refinement jobs
ccp4i2 jobs list my_project --filter "*refine*"

# Show failed jobs
ccp4i2 jobs list my_project --status FAILED
```

## Tab Completion

### Bash

Add to your `~/.bashrc`:

```bash
# ccp4i2 completion (to be implemented)
# eval "$(ccp4i2 --completion bash)"
```

### Zsh

Add to your `~/.zshrc`:

```bash
# ccp4i2 completion (to be implemented)
# eval "$(ccp4i2 --completion zsh)"
```

### PowerShell

Add to your PowerShell profile:

```powershell
# ccp4i2 completion (to be implemented)
# ccp4i2 --completion powershell | Out-String | Invoke-Expression
```

## Command Aliases

Many commands have short aliases for faster typing:

- `projects` → `project`, `proj`, `p`
- `jobs` → `job`, `j`
- `files` → `file`, `f`
- `plugins` → `plugin`, `tasks`
- `list` → `ls`
- `create` → `new`
- `show` → `info`, `get`
- `run` → `start`, `execute`
- `validate` → `check`
- `delete` → `rm`, `remove`
- `status` → `st`

Examples:
```bash
ccp4i2 p ls          # Same as: ccp4i2 projects list
ccp4i2 j ls toxd     # Same as: ccp4i2 jobs list toxd
ccp4i2 p new test    # Same as: ccp4i2 projects create test
```

## Migration from Old Commands

If you're familiar with the old management commands, here's the mapping:

| Old Command | New Command |
|-------------|-------------|
| `python manage.py list_projects` | `ccp4i2 projects list` |
| `python manage.py create_project toxd` | `ccp4i2 projects create toxd` |
| `python manage.py list_jobs toxd` | `ccp4i2 jobs list toxd` |
| `python manage.py create_job --projectname toxd --taskname refmac5` | `ccp4i2 jobs create toxd refmac5` |
| `python manage.py run_job --projectname toxd --jobnumber 1` | `ccp4i2 jobs run toxd 1` |
| `python manage.py validate_job --projectname toxd --jobnumber 1` | `ccp4i2 jobs validate toxd 1` |
| `python manage.py export_project --projectname toxd -o out.zip` | `ccp4i2 projects export toxd -o out.zip` |

The new CLI is:
- **Shorter**: Fewer flags, more intuitive arguments
- **Smarter**: Auto-detects UUIDs vs names
- **Consistent**: Same patterns across all commands
- **Discoverable**: Built-in help at every level

## Troubleshooting

### Command not found

Make sure the `ccp4i2` script is in your PATH or create a symlink:
```bash
sudo ln -s /path/to/cdata-codegen/ccp4i2 /usr/local/bin/ccp4i2
```

### Virtual environment not activated

The wrapper script should handle this automatically. If you see warnings, check:
```bash
ls -la /path/to/cdata-codegen/.venv
```

### CCP4 environment not found

Update the CCP4 path in the wrapper script (`ccp4i2`, `ccp4i2.bat`, or `ccp4i2.ps1`).

### Permission denied (Unix)

Make the script executable:
```bash
chmod +x ccp4i2
```

### Getting Help

For any command, add `--help`:
```bash
ccp4i2 --help
ccp4i2 projects --help
ccp4i2 jobs create --help
```
