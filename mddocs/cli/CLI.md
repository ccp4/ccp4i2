# CCP4i2 Command Line Interface

Modern CLI for the ccp4i2 project.

The `i2` command is installed as a console script entry point (defined in `pyproject.toml`).
It wraps Django management commands with a clean, resource-oriented syntax.

## Table of Contents

- [Installation](#installation)
- [Quick Start](#quick-start)
- [Command Structure](#command-structure)
- [Commands Reference](#commands-reference)
  - [Projects](#projects)
  - [Jobs](#jobs)
  - [Files](#files)
  - [File Types](#file-types)
  - [Report](#report)
  - [Export / Import](#export--import)
  - [Run (i2run shortcut)](#run-i2run-shortcut)
- [Examples](#examples)

## Installation

The `i2` command is registered as a console script in `server/pyproject.toml`:

```toml
[project.scripts]
i2 = "ccp4i2.cli.i2:main"
```

After installing the package (e.g. `pip install -e server/`), the `i2` command is available.
Alternatively, run directly with:

```bash
ccp4-python -m ccp4i2.cli.i2 <command>
```

## Quick Start

```bash
# List all projects
i2 projects

# Create a new project
i2 projects create my_project

# List jobs in a project
i2 jobs toxd

# Create and run a job
i2 jobs create toxd import_merged_mtz
i2 jobs run toxd 1

# Or create and run in one step via i2run
i2 run refmac --project_name toxd --hklin data.mtz --xyzin model.pdb
```

## Command Structure

```
i2 <resource> [<action>] [<args>...] [options]
```

- **Resource**: `projects`, `jobs`, `files`, `filetypes`, `report`, `export`, `import`, `run`
- **Action**: `list`, `create`, `show`, `run`, `tree`, `clone`, `cat`, etc.
- **Args**: Positional arguments (project name, job number, file path)
- **Options**: Flags passed through to the underlying management command

## Commands Reference

### Projects

```bash
i2 projects                         # List all projects (default action)
i2 projects list                    # List all projects
i2 projects create <name> [opts]    # Create a new project
i2 projects show <project>          # Show project details
i2 projects tree <project> [opts]   # Show project directory tree
i2 projects cat <project> <path> [opts]  # View file in project directory
```

Examples:
```bash
i2 projects create toxd --description "Toxd structure"
i2 projects show toxd
i2 projects tree toxd --depth 5
i2 projects cat toxd CCP4_JOBS/job_1/program.log --lines 20
```

### Jobs

```bash
i2 jobs <project>                   # List jobs in project (default action)
i2 jobs list <project>              # List jobs in project
i2 jobs create <project> <task> [opts]  # Create a new job
i2 jobs run <project> <job> [opts]  # Run a job
i2 jobs tree <project> <job> [opts] # Show job directory tree
i2 jobs clone <project> <job>       # Clone a job
i2 jobs kpi <project> <job>         # Show job KPIs (float/char values)
i2 jobs validate <project> <job> [opts]  # Validate job parameters
i2 jobs set-status <project> <job> <status>  # Change job status
i2 jobs cat <project> <job> <path> [opts]    # View file in job directory
i2 jobs set-param <project> <job> <param> <value>  # Set job parameter
i2 jobs upload-param <project> <job> <param> <file> [col]  # Upload file param
```

#### set-param FileUse Syntax

When setting file parameters, you can reference outputs from previous jobs:

```bash
i2 jobs set-param myproject 5 XYZIN '[-1].XYZOUT[0]'       # Last job's XYZOUT
i2 jobs set-param myproject 5 HKLIN 'refmac[-1].HKLOUT'     # Last refmac's HKLOUT
i2 jobs set-param myproject 5 XYZIN fullPath=/path/to/file.pdb
```

Parameter names without a dot are auto-prefixed with `inputData.` (e.g. `XYZIN` becomes `inputData.XYZIN`).

#### upload-param

Upload a file and set it as a job parameter:

```bash
i2 jobs upload-param myproject 5 XYZIN /path/to/model.pdb
i2 jobs upload-param myproject 5 F_SIGF /path/to/data.mtz '/*/*/[FP,SIGFP]'
```

Examples:
```bash
i2 jobs toxd
i2 jobs create toxd refmac5
i2 jobs run toxd 5 -y
i2 jobs validate toxd 5
i2 jobs set-status toxd 5 FINISHED
i2 jobs cat toxd 5 program.log --tail 50
i2 jobs kpi toxd 5 --json
i2 jobs clone toxd 3
```

### Files

```bash
i2 files <project> <job>            # List files in a job (default action)
i2 files list <project> <job>       # List files in a job
i2 files cat <project> <job> <name> # Display file contents
i2 files uses [opts]                # List file-job relationships
i2 files imports [opts]             # List file imports
```

The `uses` and `imports` subcommands accept filter options passed through to the
management command (e.g. `--project`, `--job`, `--format`).

Examples:
```bash
i2 files toxd 5
i2 files cat toxd 5 output.mtz
i2 files uses --project toxd --format json
i2 files imports --job <uuid>
```

### File Types

```bash
i2 filetypes [opts]                 # List registered file types
```

Options are passed through to `list_filetypes` (e.g. `--format json`, `--filter mtz`).

### Report

```bash
i2 report <project> <job> [opts]    # Generate job report
```

Examples:
```bash
i2 report toxd 5
i2 report toxd 5 --json
```

### Export / Import

```bash
i2 export job <project> <job> [opts]     # Export a job to zip
i2 export project <project> [opts]       # Export a project to zip
i2 import <zipfile> [opts]               # Import a project from zip
```

Examples:
```bash
i2 export job toxd 5 -o job5.zip
i2 export project toxd -o toxd_backup.zip --jobs "1,2,5"
i2 import toxd_backup.zip
```

### Run (i2run shortcut)

```bash
i2 run <task> [options]             # Create and run a task via i2run
```

This delegates to the `i2run` management command with all arguments passed through.

Examples:
```bash
i2 run refmac --project_name toxd --hklin data.mtz --xyzin model.pdb
```

## Examples

### Complete Workflow

```bash
# 1. Create a new project
i2 projects create toxd_structure

# 2. Create an import job
i2 jobs create toxd_structure import_merged_mtz

# 3. Validate the job
i2 jobs validate toxd_structure 1

# 4. Run the job
i2 jobs run toxd_structure 1

# 5. Check results
i2 jobs kpi toxd_structure 1
i2 report toxd_structure 1

# 6. Create a refinement job using output from import
i2 jobs create toxd_structure refmac5
i2 jobs set-param toxd_structure 2 XYZIN '[-1].XYZOUT[0]'
i2 jobs set-param toxd_structure 2 HKLIN '[-1].HKLOUT'
i2 jobs run toxd_structure 2

# 7. Export the project
i2 export project toxd_structure -o backup.zip
```

### Migration from Management Commands

| Management Command | i2 Command |
|-------------------|------------|
| `manage.py list_projects` | `i2 projects` |
| `manage.py create_project toxd` | `i2 projects create toxd` |
| `manage.py list_jobs toxd` | `i2 jobs toxd` |
| `manage.py create_job --projectname toxd --taskname refmac5` | `i2 jobs create toxd refmac5` |
| `manage.py run_job --projectname toxd --jobnumber 1` | `i2 jobs run toxd 1` |
| `manage.py validate_job --projectname toxd --jobnumber 1` | `i2 jobs validate toxd 1` |
| `manage.py export_project --projectname toxd -o out.zip` | `i2 export project toxd -o out.zip` |
| `manage.py set_job_status --projectname toxd --jobnumber 1 --status FINISHED` | `i2 jobs set-status toxd 1 FINISHED` |
| `manage.py list_filetypes` | `i2 filetypes` |
