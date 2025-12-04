# CCP4i2 CLI - Quick Reference

One-page cheat sheet for the modern CCP4i2 command line interface.

## Installation

```bash
# Unix/Linux/macOS
chmod +x ccp4i2
sudo ln -s $(pwd)/ccp4i2 /usr/local/bin/ccp4i2

# Windows - Add directory to PATH or use:
# ccp4i2.bat (cmd) or ccp4i2.ps1 (PowerShell)
```

## Command Structure

```
ccp4i2 <resource> <action> [<identifier>] [options]
```

## Common Commands

### Projects

```bash
# List all projects
ccp4i2 projects list
ccp4i2 p ls                    # Short alias

# Create a project
ccp4i2 projects create my_project --description "My work"
ccp4i2 p new my_project -d "My work"

# Show project details
ccp4i2 projects show toxd
ccp4i2 p info toxd

# Export/import
ccp4i2 projects export toxd -o backup.zip
ccp4i2 projects import backup.zip

# Show directory tree
ccp4i2 projects tree toxd
ccp4i2 p tree toxd --depth 5

# View files in project
ccp4i2 projects cat toxd CCP4_JOBS/job_1/input_params.xml
ccp4i2 p cat toxd CCP4_JOBS/job_1/program.log --tail 50
```

### Jobs

```bash
# List jobs in a project
ccp4i2 jobs list toxd
ccp4i2 j ls toxd               # Short alias

# Filter jobs
ccp4i2 jobs list toxd --status FINISHED
ccp4i2 jobs list toxd --filter "refine*"

# Create a job
ccp4i2 jobs create toxd refmac5
ccp4i2 j new toxd refmac5 --title "Refinement round 1"

# Create from another job's inputs
ccp4i2 jobs create toxd refmac5 --from-job 1

# Run a job
ccp4i2 jobs run toxd 1
ccp4i2 j run toxd 1 --detach   # Run in background

# Validate before running
ccp4i2 jobs validate toxd 1

# Export a job
ccp4i2 jobs export toxd 1 -o job1.zip

# Show directory tree
ccp4i2 jobs tree toxd 1
ccp4i2 j tree toxd 1 --depth 4

# View files in job
ccp4i2 jobs cat toxd 1 input_params.xml
ccp4i2 j cat toxd 1 program.log --tail 100
```

### Shortcuts

```bash
# Create and run in one step
ccp4i2 run toxd refmac5

# Quick status overview
ccp4i2 status
ccp4i2 status toxd
ccp4i2 st                      # Short alias
```

### Files

```bash
# Preview file contents
ccp4i2 files preview /path/to/logfile.log
ccp4i2 f preview output.mtz --lines 100
```

## Filters and Sorting

### Projects

```bash
# Filter by pattern
ccp4i2 projects list --filter "toxd*"

# Sort by different fields
ccp4i2 projects list --sort-by created
ccp4i2 projects list --sort-by jobs --reverse

# Limit results
ccp4i2 projects list --limit 10
```

### Jobs

```bash
# Filter by status
ccp4i2 jobs list toxd --status RUNNING
ccp4i2 jobs list toxd --status FINISHED

# Filter by task
ccp4i2 jobs list toxd --task refmac5

# Filter by title pattern
ccp4i2 jobs list toxd --filter "refine*"

# Sort jobs
ccp4i2 jobs list toxd --sort-by created --reverse
```

## Output Formats

```bash
# JSON output (machine-readable)
ccp4i2 --json projects list
ccp4i2 --json jobs list toxd

# Quiet mode (names/IDs only)
ccp4i2 -q projects list        # Just names
ccp4i2 -q jobs list toxd       # Just job numbers

# Verbose mode (debug info)
ccp4i2 -v projects create test
```

## Identifiers

The CLI auto-detects identifier types:

```bash
# By name
ccp4i2 projects show toxd
ccp4i2 jobs run toxd 1

# By UUID (full or short prefix)
ccp4i2 projects show a1b2c3d4-e5f6-7890-abcd-ef1234567890
ccp4i2 projects show a1b2c3d4

# Jobs can use number or UUID
ccp4i2 jobs run toxd 1
ccp4i2 jobs run toxd a1b2c3d4
```

## Aliases

### Resource Aliases
- `projects` → `project`, `proj`, `p`
- `jobs` → `job`, `j`
- `files` → `file`, `f`
- `plugins` → `plugin`, `tasks`

### Action Aliases
- `list` → `ls`
- `create` → `new`
- `show` → `info`, `get`
- `run` → `start`, `execute`
- `validate` → `check`
- `delete` → `rm`, `remove`
- `status` → `st`

### Examples
```bash
ccp4i2 p ls                    # projects list
ccp4i2 p new test              # projects create test
ccp4i2 j ls toxd               # jobs list toxd
ccp4i2 j new toxd refmac5      # jobs create toxd refmac5
ccp4i2 st                      # status
```

## Global Flags

Work with any command:

```bash
--json          # JSON output
-q, --quiet     # Minimal output
-v, --verbose   # Debug info
--no-color      # Disable colors
-h, --help      # Show help
```

## Help System

```bash
# Top-level help
ccp4i2 --help
ccp4i2 -h

# Resource-level help
ccp4i2 projects --help
ccp4i2 jobs --help

# Action-level help
ccp4i2 projects list --help
ccp4i2 jobs create --help
```

## Complete Workflow Example

```bash
# 1. Create project
ccp4i2 projects create my_structure -d "Structure determination"

# 2. Import data
ccp4i2 jobs create my_structure import_merged_mtz --title "Import MTZ"

# 3. Validate and run
ccp4i2 jobs validate my_structure 1
ccp4i2 jobs run my_structure 1 --detach

# 4. Check status
ccp4i2 jobs list my_structure --status RUNNING

# 5. Create next job using previous output
ccp4i2 jobs create my_structure refmac5 --from-job 1

# 6. Backup project
ccp4i2 projects export my_structure -o backup.zip
```

## Scripting Examples

```bash
# Get all finished jobs as JSON
ccp4i2 --json jobs list toxd --status FINISHED > finished.json

# Get project UUID
project_uuid=$(ccp4i2 -q projects create test)
echo "Created: $project_uuid"

# Process all projects
for project in $(ccp4i2 -q projects list); do
  echo "Processing $project"
  ccp4i2 jobs list "$project"
done

# Count jobs by status
ccp4i2 --json jobs list toxd | jq '.[] | .status' | sort | uniq -c
```

## Troubleshooting

```bash
# Check if command is accessible
which ccp4i2

# Test with verbose mode
ccp4i2 -v projects list

# Get help on any command
ccp4i2 <command> --help
```

## Migration from Old Commands

| Old | New |
|-----|-----|
| `python manage.py list_projects` | `ccp4i2 projects list` |
| `python manage.py create_project toxd` | `ccp4i2 projects create toxd` |
| `python manage.py list_jobs toxd` | `ccp4i2 jobs list toxd` |
| `python manage.py run_job --projectname toxd --jobnumber 1` | `ccp4i2 jobs run toxd 1` |

## Tips

1. **Use aliases** for faster typing: `ccp4i2 p ls` vs `ccp4i2 projects list`
2. **Tab completion** (coming soon) will make discovery easier
3. **JSON output** is perfect for scripting: `ccp4i2 --json projects list`
4. **Help is everywhere**: Try `--help` on any command level
5. **UUIDs are flexible**: Use full UUID or just the first 8 characters
