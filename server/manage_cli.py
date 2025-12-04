#!/usr/bin/env python3
"""
Modern CLI for cdata-codegen (CCP4i2 replacement).

This is the main entry point for the user-facing CLI. It provides a clean,
resource-oriented interface following modern CLI design patterns.

Architecture:
    - Resource-first: ccp4i2 <resource> <action> [args]
    - Smart identifier resolution (name, UUID, ID)
    - Consistent flags across all commands
    - Rich human output + JSON for scripting
    - Contextual help at every level

Usage:
    ccp4i2 projects list
    ccp4i2 jobs create toxd import_merged_mtz
    ccp4i2 jobs run toxd 1 --detach
"""

import sys
import os
import argparse
from pathlib import Path
from typing import List, Optional

# Ensure Django settings are configured
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'ccp4x.config.settings')

# Add server directory to path
SERVER_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SERVER_DIR))

# Ensure CCP4I2_ROOT is set for plugin discovery
# The ccp4i2 wrapper script should set this, but if not, detect it
if 'CCP4I2_ROOT' not in os.environ:
    # Server is in cdata-codegen/server, so project root is parent
    project_root = SERVER_DIR.parent
    os.environ['CCP4I2_ROOT'] = str(project_root)
    print(f"DEBUG: Auto-detected CCP4I2_ROOT={os.environ['CCP4I2_ROOT']}", file=sys.stderr)

import django
django.setup()

from django.core.management import call_command


class CLI:
    """Main CLI dispatcher."""

    def __init__(self):
        self.parser = self._build_parser()

    def _build_parser(self) -> argparse.ArgumentParser:
        """Build the main argument parser with all subcommands."""
        parser = argparse.ArgumentParser(
            prog='ccp4i2',
            description='Modern CLI for CCP4i2 crystallographic workflow management',
            formatter_class=argparse.RawDescriptionHelpFormatter,
            epilog="""
Examples:
  # Projects
  ccp4i2 projects list
  ccp4i2 projects create my_project --description "Crystal structure"
  ccp4i2 projects show toxd

  # Jobs
  ccp4i2 jobs list toxd
  ccp4i2 jobs create toxd import_merged_mtz --title "Import MTZ"
  ccp4i2 jobs run toxd 1 --detach

  # Shortcuts
  ccp4i2 run toxd import_merged_mtz  # create + run in one step

For help on any command:
  ccp4i2 <resource> --help
  ccp4i2 <resource> <action> --help
            """
        )

        # Global flags
        parser.add_argument('--json', action='store_true',
                           help='Output as JSON (machine-readable)')
        parser.add_argument('-q', '--quiet', action='store_true',
                           help='Minimal output')
        parser.add_argument('-v', '--verbose', action='store_true',
                           help='Verbose output with debug info')
        parser.add_argument('--no-color', action='store_true',
                           help='Disable colored output')

        # Subcommands
        subparsers = parser.add_subparsers(dest='resource', help='Resource type')

        # Projects
        self._add_projects_commands(subparsers)

        # Jobs
        self._add_jobs_commands(subparsers)

        # Files
        self._add_files_commands(subparsers)

        # Plugins
        self._add_plugins_commands(subparsers)

        # File types
        self._add_filetypes_commands(subparsers)

        # Shortcuts
        self._add_shortcut_commands(subparsers)

        # Legacy i2run compatibility
        self._add_i2run_command(subparsers)

        return parser

    def _add_projects_commands(self, subparsers):
        """Add all project-related commands."""
        projects = subparsers.add_parser('projects', aliases=['project', 'proj', 'p'],
                                         help='Project management')
        projects_sub = projects.add_subparsers(dest='action', help='Project action')

        # projects list
        list_cmd = projects_sub.add_parser('list', aliases=['ls'],
                                           help='List all projects')
        list_cmd.add_argument('--filter', type=str,
                             help='Filter by name pattern (wildcards supported)')
        list_cmd.add_argument('--sort-by', choices=['name', 'created', 'accessed', 'jobs'],
                             default='created', help='Sort field')
        list_cmd.add_argument('--reverse', action='store_true',
                             help='Reverse sort order')
        list_cmd.add_argument('--limit', type=int,
                             help='Limit number of results')

        # projects create
        create_cmd = projects_sub.add_parser('create', aliases=['new'],
                                             help='Create a new project')
        create_cmd.add_argument('name', help='Project name')
        create_cmd.add_argument('-d', '--description', default='',
                               help='Project description')
        create_cmd.add_argument('--directory', default='__default__',
                               help='Custom project directory')

        # projects show
        show_cmd = projects_sub.add_parser('show', aliases=['info', 'get'],
                                           help='Show project details')
        show_cmd.add_argument('project', help='Project name or UUID')

        # projects export
        export_cmd = projects_sub.add_parser('export',
                                             help='Export project to ZIP')
        export_cmd.add_argument('project', help='Project name or UUID')
        export_cmd.add_argument('-o', '--output', required=True,
                               help='Output ZIP file path')
        export_cmd.add_argument('--jobs', type=str,
                               help='Comma-separated job numbers to export')
        export_cmd.add_argument('--detach', action='store_true',
                               help='Run export in background')

        # projects import
        import_cmd = projects_sub.add_parser('import',
                                             help='Import project from ZIP')
        import_cmd.add_argument('zip_file', help='ZIP file to import')

        # projects tree
        tree_cmd = projects_sub.add_parser('tree',
                                           help='Show directory tree structure')
        tree_cmd.add_argument('project', help='Project name or UUID')
        tree_cmd.add_argument('--depth', type=int, default=3,
                             help='Maximum depth to display (default: 3, 0 = unlimited)')
        tree_cmd.add_argument('--no-sizes', action='store_true',
                             help='Do not show file sizes')
        tree_cmd.add_argument('--show-hidden', action='store_true',
                             help='Show hidden files and directories')

        # projects cat
        cat_cmd = projects_sub.add_parser('cat', aliases=['view', 'show-file'],
                                          help='View file contents within project')
        cat_cmd.add_argument('project', help='Project name or UUID')
        cat_cmd.add_argument('file_path', help='Relative path to file within project')
        cat_cmd.add_argument('--lines', type=int,
                            help='Number of lines to show from start')
        cat_cmd.add_argument('--tail', type=int,
                            help='Number of lines to show from end')
        cat_cmd.add_argument('--binary', action='store_true',
                            help='Allow viewing binary files (hex dump)')
        cat_cmd.add_argument('--no-line-numbers', action='store_true',
                            help='Do not show line numbers')

        # projects delete
        delete_cmd = projects_sub.add_parser('delete', aliases=['rm', 'remove'],
                                             help='Delete a project')
        delete_cmd.add_argument('project', help='Project name or UUID')
        delete_cmd.add_argument('--force', action='store_true',
                               help='Skip confirmation')

    def _add_jobs_commands(self, subparsers):
        """Add all job-related commands."""
        jobs = subparsers.add_parser('jobs', aliases=['job', 'j'],
                                     help='Job management')
        jobs_sub = jobs.add_subparsers(dest='action', help='Job action')

        # jobs list
        list_cmd = jobs_sub.add_parser('list', aliases=['ls'],
                                       help='List jobs in a project')
        list_cmd.add_argument('project', help='Project name or UUID')
        list_cmd.add_argument('--status', choices=[
            'UNKNOWN', 'PENDING', 'QUEUED', 'RUNNING', 'INTERRUPTED',
            'FAILED', 'FINISHED', 'RUNNING_REMOTELY', 'FILE_HOLDER',
            'TO_DELETE', 'UNSATISFACTORY'
        ], help='Filter by status')
        list_cmd.add_argument('--task', type=str,
                             help='Filter by task name')
        list_cmd.add_argument('--filter', type=str,
                             help='Filter by title pattern')
        list_cmd.add_argument('--sort-by',
                             choices=['number', 'title', 'task', 'status', 'created', 'finished'],
                             default='number', help='Sort field')
        list_cmd.add_argument('--reverse', action='store_true',
                             help='Reverse sort order')
        list_cmd.add_argument('--limit', type=int,
                             help='Limit number of results')

        # jobs create
        create_cmd = jobs_sub.add_parser('create', aliases=['new'],
                                         help='Create a new job')
        create_cmd.add_argument('project', help='Project name or UUID')
        create_cmd.add_argument('task', help='Task/plugin name')
        create_cmd.add_argument('--from-job', type=str, dest='from_job',
                               help='Context job (name or UUID) to copy inputs from')
        create_cmd.add_argument('--title', type=str,
                               help='Job title')

        # jobs show
        show_cmd = jobs_sub.add_parser('show', aliases=['info', 'get'],
                                       help='Show job details')
        show_cmd.add_argument('project', help='Project name or UUID')
        show_cmd.add_argument('job', help='Job number or UUID')

        # jobs run
        run_cmd = jobs_sub.add_parser('run', aliases=['start', 'execute'],
                                      help='Run a job')
        run_cmd.add_argument('project', help='Project name or UUID')
        run_cmd.add_argument('job', help='Job number or UUID')
        run_cmd.add_argument('--detach', action='store_true',
                            help='Run in background')

        # jobs validate
        validate_cmd = jobs_sub.add_parser('validate', aliases=['check'],
                                           help='Validate job parameters')
        validate_cmd.add_argument('project', help='Project name or UUID')
        validate_cmd.add_argument('job', help='Job number or UUID')
        validate_cmd.add_argument('-o', '--output', type=str,
                                 help='Save validation report to file')

        # jobs export
        export_cmd = jobs_sub.add_parser('export',
                                         help='Export job to ZIP')
        export_cmd.add_argument('project', help='Project name or UUID')
        export_cmd.add_argument('job', help='Job number or UUID')
        export_cmd.add_argument('-o', '--output', required=True,
                               help='Output ZIP file path')

        # jobs set-status
        status_cmd = jobs_sub.add_parser('set-status',
                                         help='Set job status')
        status_cmd.add_argument('project', help='Project name or UUID')
        status_cmd.add_argument('job', help='Job number or UUID')
        status_cmd.add_argument('status', help='New status')

        # jobs set-param
        param_cmd = jobs_sub.add_parser('set-param',
                                        help='Set job parameter')
        param_cmd.add_argument('project', help='Project name or UUID')
        param_cmd.add_argument('job', help='Job number or UUID')
        param_cmd.add_argument('param', help='Parameter path')
        param_cmd.add_argument('value', help='Parameter value')

        # jobs report
        report_cmd = jobs_sub.add_parser('report',
                                         help='Get job report')
        report_cmd.add_argument('project', help='Project name or UUID')
        report_cmd.add_argument('job', help='Job number or UUID')
        report_cmd.add_argument('--format', choices=['text', 'json', 'xml'],
                               default='text', help='Report format')

        # jobs tree
        tree_cmd = jobs_sub.add_parser('tree',
                                       help='Show directory tree structure')
        tree_cmd.add_argument('project', help='Project name or UUID')
        tree_cmd.add_argument('job', help='Job number or UUID')
        tree_cmd.add_argument('--depth', type=int, default=2,
                             help='Maximum depth to display (default: 2, 0 = unlimited)')
        tree_cmd.add_argument('--no-sizes', action='store_true',
                             help='Do not show file sizes')
        tree_cmd.add_argument('--show-hidden', action='store_true',
                             help='Show hidden files and directories')

        # jobs cat
        cat_cmd = jobs_sub.add_parser('cat', aliases=['view', 'show-file'],
                                      help='View file contents within job')
        cat_cmd.add_argument('project', help='Project name or UUID')
        cat_cmd.add_argument('job', help='Job number or UUID')
        cat_cmd.add_argument('file_path', help='Relative path to file within job directory')
        cat_cmd.add_argument('--lines', type=int,
                            help='Number of lines to show from start')
        cat_cmd.add_argument('--tail', type=int,
                            help='Number of lines to show from end')
        cat_cmd.add_argument('--binary', action='store_true',
                            help='Allow viewing binary files (hex dump)')
        cat_cmd.add_argument('--no-line-numbers', action='store_true',
                            help='Do not show line numbers')

        # jobs delete
        delete_cmd = jobs_sub.add_parser('delete', aliases=['rm', 'remove'],
                                         help='Delete a job')
        delete_cmd.add_argument('project', help='Project name or UUID')
        delete_cmd.add_argument('job', help='Job number or UUID')
        delete_cmd.add_argument('--force', action='store_true',
                               help='Skip confirmation')

        # jobs clone
        clone_cmd = jobs_sub.add_parser('clone', aliases=['copy', 'duplicate'],
                                        help='Clone a job with identical parameters')
        clone_cmd.add_argument('project', help='Project name or UUID')
        clone_cmd.add_argument('job', help='Job number or UUID to clone')
        clone_cmd.add_argument('--json', action='store_true',
                              help='Output result as JSON')

        # jobs files
        files_cmd = jobs_sub.add_parser('files', aliases=['list-files'],
                                        help='List files used by a job')
        files_cmd.add_argument('project', help='Project name or UUID')
        files_cmd.add_argument('job', help='Job number or UUID')
        files_cmd.add_argument('--format', choices=['table', 'csv', 'json'],
                              default='table', help='Output format')

        # jobs file-uses
        uses_cmd = jobs_sub.add_parser('file-uses', aliases=['uses'],
                                       help='List file uses for a job')
        uses_cmd.add_argument('project', help='Project name or UUID')
        uses_cmd.add_argument('job', help='Job number or UUID')
        uses_cmd.add_argument('--format', choices=['table', 'csv', 'json'],
                             default='table', help='Output format')

        # jobs file-imports
        imports_cmd = jobs_sub.add_parser('file-imports', aliases=['imports'],
                                          help='List file imports for a job')
        imports_cmd.add_argument('project', help='Project name or UUID')
        imports_cmd.add_argument('job', help='Job number or UUID')
        imports_cmd.add_argument('--format', choices=['table', 'csv', 'json'],
                                 default='table', help='Output format')

    def _add_files_commands(self, subparsers):
        """Add file-related commands."""
        files = subparsers.add_parser('files', aliases=['file', 'f'],
                                      help='File operations')
        files_sub = files.add_subparsers(dest='action', help='File action')

        # files list
        list_cmd = files_sub.add_parser('list', aliases=['ls'],
                                        help='List files in database')
        list_cmd.add_argument('--project', type=str,
                             help='Project name or UUID')
        list_cmd.add_argument('--job', type=str,
                             help='Job UUID')
        list_cmd.add_argument('--all', action='store_true',
                             help='List all files across all projects')
        list_cmd.add_argument('--format', choices=['table', 'csv', 'json'],
                             default='table', help='Output format')
        list_cmd.add_argument('--filter', type=str,
                             help='Filter by filename')

        # files uses
        uses_cmd = files_sub.add_parser('uses', aliases=['file-uses'],
                                        help='List file uses (file-job relationships)')
        uses_cmd.add_argument('--job', type=str, help='Job UUID')
        uses_cmd.add_argument('--file', type=str, help='File UUID')
        uses_cmd.add_argument('--project', type=str, help='Project name or UUID')
        uses_cmd.add_argument('--format', choices=['table', 'csv', 'json'],
                             default='table', help='Output format')

        # files imports
        imports_cmd = files_sub.add_parser('imports', aliases=['file-imports'],
                                           help='List file imports')
        imports_cmd.add_argument('--job', type=str, help='Job UUID')
        imports_cmd.add_argument('--file', type=str, help='File UUID')
        imports_cmd.add_argument('--project', type=str, help='Project name or UUID')
        imports_cmd.add_argument('--format', choices=['table', 'csv', 'json'],
                                 default='table', help='Output format')

        # files preview
        preview_cmd = files_sub.add_parser('preview', aliases=['show', 'cat'],
                                           help='Preview file contents')
        preview_cmd.add_argument('path', help='File path')
        preview_cmd.add_argument('--lines', type=int, default=50,
                                help='Number of lines to show')

    def _add_plugins_commands(self, subparsers):
        """Add plugin-related commands."""
        plugins = subparsers.add_parser('plugins', aliases=['plugin', 'tasks'],
                                        help='Plugin/task management')
        plugins_sub = plugins.add_subparsers(dest='action', help='Plugin action')

        # plugins list
        list_cmd = plugins_sub.add_parser('list', aliases=['ls'],
                                          help='List available plugins')
        list_cmd.add_argument('--filter', type=str,
                             help='Filter by name pattern')

        # plugins show
        show_cmd = plugins_sub.add_parser('show', aliases=['info'],
                                          help='Show plugin details')
        show_cmd.add_argument('plugin', help='Plugin name')

    def _add_filetypes_commands(self, subparsers):
        """Add file type commands."""
        filetypes = subparsers.add_parser('filetypes', aliases=['filetype', 'ft'],
                                          help='File type information')
        filetypes_sub = filetypes.add_subparsers(dest='action', help='File type action')

        # filetypes list
        list_cmd = filetypes_sub.add_parser('list', aliases=['ls'],
                                            help='List registered file types')
        list_cmd.add_argument('--format', type=str, choices=['table', 'csv', 'json'],
                             default='table', help='Output format')
        list_cmd.add_argument('--filter', type=str,
                             help='Filter by MIME type name')

    def _add_shortcut_commands(self, subparsers):
        """Add convenience shortcut commands."""
        # run (create + execute in one)
        run_cmd = subparsers.add_parser('run',
                                        help='Create and run a job in one step')
        run_cmd.add_argument('project', help='Project name or UUID')
        run_cmd.add_argument('task', help='Task/plugin name')
        run_cmd.add_argument('--from-job', type=str, dest='from_job',
                            help='Context job to copy inputs from')
        run_cmd.add_argument('--detach', action='store_true',
                            help='Run in background')
        run_cmd.add_argument('--title', type=str,
                            help='Job title')

        # status (quick overview)
        status_cmd = subparsers.add_parser('status', aliases=['st'],
                                           help='Show project/job status overview')
        status_cmd.add_argument('project', nargs='?',
                               help='Project name or UUID (optional, shows all if omitted)')

    def _add_i2run_command(self, subparsers):
        """Add legacy i2run compatibility command."""
        i2run_cmd = subparsers.add_parser('i2run',
                                          help='Legacy i2run command-line runner')
        i2run_cmd.add_argument('task_name', type=str,
                              help='Task name (e.g., prosmart_refmac)')
        i2run_cmd.add_argument('--project_name', type=str, required=True,
                              help='Project name or UUID')
        i2run_cmd.add_argument('--project_path', type=str,
                              help='Project path (optional)')
        # All other arguments are passed through to the plugin

    def run(self, args: Optional[List[str]] = None):
        """Run the CLI with the given arguments."""
        # For i2run, use parse_known_args to allow plugin-specific parameters
        # First, peek at the command to see if it's i2run
        import sys
        argv_to_check = args if args is not None else sys.argv[1:]

        # Check if 'i2run' is in the arguments (may be after global flags like -v)
        is_i2run = 'i2run' in argv_to_check

        if is_i2run:
            # Use parse_known_args to preserve unknown arguments for plugin
            parsed_args, unknown_args = self.parser.parse_known_args(args)
            # Store unknown args for i2run handler
            parsed_args.unknown_args = unknown_args
        else:
            # Normal parsing for other commands
            parsed_args = self.parser.parse_args(args)
            parsed_args.unknown_args = []

        # Handle no command
        if not parsed_args.resource:
            self.parser.print_help()
            return 1

        # Dispatch to appropriate handler
        try:
            return self._dispatch(parsed_args)
        except KeyboardInterrupt:
            print("\n\nInterrupted by user")
            return 130
        except Exception as e:
            if parsed_args.verbose:
                import traceback
                traceback.print_exc()
            else:
                print(f"Error: {e}", file=sys.stderr)
            return 1

    def _dispatch(self, args):
        """Dispatch to the appropriate command handler."""
        resource = args.resource
        action = getattr(args, 'action', None)

        # Map to Django management commands
        # This is a temporary bridge - eventually we'll refactor management commands
        # to be called directly from here

        if resource in ('projects', 'project', 'proj', 'p'):
            return self._handle_projects(args)
        elif resource in ('jobs', 'job', 'j'):
            return self._handle_jobs(args)
        elif resource in ('files', 'file', 'f'):
            return self._handle_files(args)
        elif resource in ('plugins', 'plugin', 'tasks'):
            return self._handle_plugins(args)
        elif resource in ('filetypes', 'filetype', 'ft'):
            return self._handle_filetypes(args)
        elif resource == 'run':
            return self._handle_run_shortcut(args)
        elif resource in ('status', 'st'):
            return self._handle_status_shortcut(args)
        elif resource == 'i2run':
            return self._handle_i2run(args)
        else:
            print(f"Unknown resource: {resource}", file=sys.stderr)
            return 1

    def _handle_projects(self, args):
        """Handle project commands."""
        action = args.action

        if not action:
            print("Error: No action specified for projects", file=sys.stderr)
            print("Try: ccp4i2 projects --help", file=sys.stderr)
            return 1

        # Map to existing management commands
        if action in ('list', 'ls'):
            cmd_args = ['list_projects']
            if args.json:
                cmd_args.append('--json')
            if args.quiet:
                cmd_args.append('--names-only')
            if hasattr(args, 'filter') and args.filter:
                cmd_args.extend(['--filter', args.filter])
            if hasattr(args, 'sort_by'):
                # Map friendly names to actual fields
                sort_map = {'created': 'creation_time', 'accessed': 'last_access'}
                sort_field = sort_map.get(args.sort_by, args.sort_by)
                cmd_args.extend(['--sort-by', sort_field])
            if hasattr(args, 'reverse') and args.reverse:
                cmd_args.append('--reverse')
            if hasattr(args, 'limit') and args.limit:
                cmd_args.extend(['--limit', str(args.limit)])

            call_command(*cmd_args)
            return 0

        elif action in ('create', 'new'):
            cmd_args = ['create_project', args.name]
            if args.description:
                cmd_args.extend(['--description', args.description])
            if args.directory != '__default__':
                cmd_args.extend(['--directory', args.directory])
            if args.json:
                cmd_args.append('--json')
            elif args.quiet:
                cmd_args.append('--quiet')

            call_command(*cmd_args)
            return 0

        elif action in ('show', 'info', 'get'):
            cmd_args = ['list_project', '--projectname', args.project]
            if args.json:
                cmd_args.append('--json')

            call_command(*cmd_args)
            return 0

        elif action == 'export':
            cmd_args = ['export_project']
            # Smart project identifier
            if self._is_uuid(args.project):
                cmd_args.extend(['--projectuuid', args.project])
            else:
                cmd_args.extend(['--projectname', args.project])
            cmd_args.extend(['--output', args.output])
            if hasattr(args, 'jobs') and args.jobs:
                cmd_args.extend(['--jobs', args.jobs])
            if hasattr(args, 'detach') and args.detach:
                cmd_args.append('--detach')

            call_command(*cmd_args)
            return 0

        elif action == 'import':
            cmd_args = ['import_ccp4_project_zip', args.zip_file]
            call_command(*cmd_args)
            return 0

        elif action == 'tree':
            cmd_args = ['tree_project', args.project]
            if hasattr(args, 'depth') and args.depth != 3:
                cmd_args.extend(['--depth', str(args.depth)])
            if hasattr(args, 'no_sizes') and args.no_sizes:
                cmd_args.append('--no-sizes')
            if hasattr(args, 'show_hidden') and args.show_hidden:
                cmd_args.append('--show-hidden')

            call_command(*cmd_args)
            return 0

        elif action in ('cat', 'view', 'show-file'):
            cmd_args = ['cat_project_file', args.project, args.file_path]
            if hasattr(args, 'lines') and args.lines:
                cmd_args.extend(['--lines', str(args.lines)])
            if hasattr(args, 'tail') and args.tail:
                cmd_args.extend(['--tail', str(args.tail)])
            if hasattr(args, 'binary') and args.binary:
                cmd_args.append('--binary')
            if hasattr(args, 'no_line_numbers') and args.no_line_numbers:
                cmd_args.append('--no-line-numbers')

            call_command(*cmd_args)
            return 0

        elif action in ('delete', 'rm', 'remove'):
            # This command doesn't exist yet - would need to be implemented
            print("Error: Project deletion not yet implemented", file=sys.stderr)
            return 1

        else:
            print(f"Unknown action: {action}", file=sys.stderr)
            return 1

    def _handle_jobs(self, args):
        """Handle job commands."""
        action = args.action

        if not action:
            print("Error: No action specified for jobs", file=sys.stderr)
            print("Try: ccp4i2 jobs --help", file=sys.stderr)
            return 1

        if action in ('list', 'ls'):
            cmd_args = ['list_jobs', args.project]
            if args.json:
                cmd_args.append('--json')
            if args.quiet:
                cmd_args.append('--numbers-only')
            if hasattr(args, 'status') and args.status:
                cmd_args.extend(['--status', args.status])
            if hasattr(args, 'task') and args.task:
                cmd_args.extend(['--task', args.task])
            if hasattr(args, 'filter') and args.filter:
                cmd_args.extend(['--filter', args.filter])
            if hasattr(args, 'sort_by'):
                # Map friendly names
                sort_map = {'created': 'creation_time', 'finished': 'finish_time', 'task': 'task_name'}
                sort_field = sort_map.get(args.sort_by, args.sort_by)
                cmd_args.extend(['--sort-by', sort_field])
            if hasattr(args, 'reverse') and args.reverse:
                cmd_args.append('--reverse')
            if hasattr(args, 'limit') and args.limit:
                cmd_args.extend(['--limit', str(args.limit)])

            call_command(*cmd_args)
            return 0

        elif action in ('create', 'new'):
            cmd_args = ['create_job']
            if self._is_uuid(args.project):
                cmd_args.extend(['--projectuuid', args.project])
            else:
                cmd_args.extend(['--projectname', args.project])
            cmd_args.extend(['--taskname', args.task])
            if hasattr(args, 'from_job') and args.from_job:
                if self._is_uuid(args.from_job):
                    cmd_args.extend(['--contextjobuuid', args.from_job])
                else:
                    cmd_args.extend(['--contextjobnumber', args.from_job])

            call_command(*cmd_args)
            return 0

        elif action in ('show', 'info', 'get'):
            # Not implemented yet - would show detailed job info
            print("Error: Job show not yet implemented", file=sys.stderr)
            return 1

        elif action in ('run', 'start', 'execute'):
            cmd_args = ['run_job']
            if self._is_uuid(args.job):
                cmd_args.extend(['--jobuuid', args.job])
            else:
                cmd_args.extend(['--projectname', args.project, '--jobnumber', args.job])
            if hasattr(args, 'detach') and args.detach:
                cmd_args.append('--detach')

            call_command(*cmd_args)
            return 0

        elif action in ('validate', 'check'):
            cmd_args = ['validate_job']
            if self._is_uuid(args.job):
                cmd_args.extend(['--jobuuid', args.job])
            else:
                cmd_args.extend(['--projectname', args.project, '--jobnumber', args.job])
            if hasattr(args, 'output') and args.output:
                cmd_args.extend(['--output', args.output])
            if args.json:
                cmd_args.append('--json')

            call_command(*cmd_args)
            return 0

        elif action == 'export':
            cmd_args = ['export_job']
            if self._is_uuid(args.job):
                cmd_args.extend(['--jobuuid', args.job])
            else:
                cmd_args.extend(['--projectname', args.project, '--jobnumber', args.job])
            cmd_args.extend(['--output', args.output])
            if args.json:
                cmd_args.append('--json')

            call_command(*cmd_args)
            return 0

        elif action == 'set-status':
            cmd_args = ['set_job_status']
            if self._is_uuid(args.job):
                cmd_args.extend(['--jobuuid', args.job])
            else:
                cmd_args.extend(['--projectname', args.project, '--jobnumber', args.job])
            cmd_args.extend(['--status', args.status])

            call_command(*cmd_args)
            return 0

        elif action == 'set-param':
            cmd_args = ['set_job_parameter']
            if self._is_uuid(args.job):
                cmd_args.extend(['--jobuuid', args.job])
            else:
                cmd_args.extend(['--projectname', args.project, '--jobnumber', args.job])
            cmd_args.extend(['--path', args.param, '--value', args.value])

            call_command(*cmd_args)
            return 0

        elif action == 'report':
            cmd_args = ['get_job_report']
            if self._is_uuid(args.job):
                cmd_args.extend(['--jobuuid', args.job])
            else:
                cmd_args.extend(['--projectname', args.project, '--jobnumber', args.job])
            if hasattr(args, 'format') and args.format == 'json':
                cmd_args.append('--json')

            call_command(*cmd_args)
            return 0

        elif action == 'tree':
            cmd_args = ['tree_job', args.project, args.job]
            if hasattr(args, 'depth') and args.depth != 2:
                cmd_args.extend(['--depth', str(args.depth)])
            if hasattr(args, 'no_sizes') and args.no_sizes:
                cmd_args.append('--no-sizes')
            if hasattr(args, 'show_hidden') and args.show_hidden:
                cmd_args.append('--show-hidden')

            call_command(*cmd_args)
            return 0

        elif action in ('cat', 'view', 'show-file'):
            cmd_args = ['cat_job_file', args.project, args.job, args.file_path]
            if hasattr(args, 'lines') and args.lines:
                cmd_args.extend(['--lines', str(args.lines)])
            if hasattr(args, 'tail') and args.tail:
                cmd_args.extend(['--tail', str(args.tail)])
            if hasattr(args, 'binary') and args.binary:
                cmd_args.append('--binary')
            if hasattr(args, 'no_line_numbers') and args.no_line_numbers:
                cmd_args.append('--no-line-numbers')

            call_command(*cmd_args)
            return 0

        elif action in ('clone', 'copy', 'duplicate'):
            # Build command arguments for clone_job
            cmd_args = ['clone_job']

            # Determine if project is UUID or name
            if self._is_uuid(args.project):
                cmd_args.extend(['--projectuuid', args.project])
            else:
                cmd_args.extend(['--projectname', args.project])

            # Determine if job is UUID or number
            if self._is_uuid(args.job):
                cmd_args.extend(['--jobuuid', args.job])
            else:
                cmd_args.extend(['--jobnumber', args.job])

            # Add JSON flag if requested
            if hasattr(args, 'json') and args.json:
                cmd_args.append('--json')

            call_command(*cmd_args)
            return 0

        elif action in ('files', 'list-files'):
            # Get job UUID to pass to list_files
            job_uuid = self._resolve_job_uuid(args.project, args.job)
            if not job_uuid:
                return 1

            cmd_args = ['list_files', '--job', str(job_uuid)]
            if hasattr(args, 'format') and args.format:
                cmd_args.extend(['--format', args.format])

            call_command(*cmd_args)
            return 0

        elif action in ('file-uses', 'uses'):
            # Get job UUID
            job_uuid = self._resolve_job_uuid(args.project, args.job)
            if not job_uuid:
                return 1

            cmd_args = ['list_file_uses', '--job', str(job_uuid)]
            if hasattr(args, 'format') and args.format:
                cmd_args.extend(['--format', args.format])

            call_command(*cmd_args)
            return 0

        elif action in ('file-imports', 'imports'):
            # Get job UUID
            job_uuid = self._resolve_job_uuid(args.project, args.job)
            if not job_uuid:
                return 1

            cmd_args = ['list_file_imports', '--job', str(job_uuid)]
            if hasattr(args, 'format') and args.format:
                cmd_args.extend(['--format', args.format])

            call_command(*cmd_args)
            return 0

        elif action in ('delete', 'rm', 'remove'):
            print("Error: Job deletion not yet implemented", file=sys.stderr)
            return 1

        else:
            print(f"Unknown action: {action}", file=sys.stderr)
            return 1

    def _handle_files(self, args):
        """Handle file commands."""
        action = args.action

        if not action:
            print("Error: No action specified for files", file=sys.stderr)
            return 1

        if action in ('list', 'ls'):
            cmd_args = ['list_files']
            if hasattr(args, 'project') and args.project:
                cmd_args.extend(['--project', args.project])
            if hasattr(args, 'job') and args.job:
                cmd_args.extend(['--job', args.job])
            if hasattr(args, 'all') and args.all:
                cmd_args.append('--all')
            if hasattr(args, 'format') and args.format:
                cmd_args.extend(['--format', args.format])
            if hasattr(args, 'filter') and args.filter:
                cmd_args.extend(['--filter', args.filter])

            call_command(*cmd_args)
            return 0

        elif action in ('uses', 'file-uses'):
            cmd_args = ['list_file_uses']
            if hasattr(args, 'job') and args.job:
                cmd_args.extend(['--job', args.job])
            if hasattr(args, 'file') and args.file:
                cmd_args.extend(['--file', args.file])
            if hasattr(args, 'project') and args.project:
                cmd_args.extend(['--project', args.project])
            if hasattr(args, 'format') and args.format:
                cmd_args.extend(['--format', args.format])

            call_command(*cmd_args)
            return 0

        elif action in ('imports', 'file-imports'):
            cmd_args = ['list_file_imports']
            if hasattr(args, 'job') and args.job:
                cmd_args.extend(['--job', args.job])
            if hasattr(args, 'file') and args.file:
                cmd_args.extend(['--file', args.file])
            if hasattr(args, 'project') and args.project:
                cmd_args.extend(['--project', args.project])
            if hasattr(args, 'format') and args.format:
                cmd_args.extend(['--format', args.format])

            call_command(*cmd_args)
            return 0

        elif action in ('preview', 'show', 'cat'):
            cmd_args = ['preview_file', args.path]
            if hasattr(args, 'lines'):
                cmd_args.extend(['--lines', str(args.lines)])

            call_command(*cmd_args)
            return 0

        else:
            print(f"Unknown action: {action}", file=sys.stderr)
            return 1

    def _handle_plugins(self, args):
        """Handle plugin commands."""
        action = args.action

        if not action:
            print("Error: No action specified for plugins", file=sys.stderr)
            return 1

        # These would need new management commands
        print(f"Error: Plugin commands not yet implemented", file=sys.stderr)
        return 1

    def _handle_filetypes(self, args):
        """Handle file type commands."""
        action = args.action

        if not action:
            print("Error: No action specified for filetypes", file=sys.stderr)
            print("Try: ccp4i2 filetypes --help", file=sys.stderr)
            return 1

        if action in ('list', 'ls'):
            cmd_args = ['list_filetypes']
            if hasattr(args, 'format') and args.format:
                cmd_args.extend(['--format', args.format])
            if hasattr(args, 'filter') and args.filter:
                cmd_args.extend(['--filter', args.filter])

            call_command(*cmd_args)
            return 0

        else:
            print(f"Unknown action: {action}", file=sys.stderr)
            return 1

    def _handle_run_shortcut(self, args):
        """Handle the 'run' shortcut (create + run)."""
        # First create the job
        print(f"Creating job in project '{args.project}' with task '{args.task}'...")
        create_args = ['create_job']
        if self._is_uuid(args.project):
            create_args.extend(['--projectuuid', args.project])
        else:
            create_args.extend(['--projectname', args.project])
        create_args.extend(['--taskname', args.task])
        if hasattr(args, 'from_job') and args.from_job:
            if self._is_uuid(args.from_job):
                create_args.extend(['--contextjobuuid', args.from_job])
            else:
                create_args.extend(['--contextjobnumber', args.from_job])

        # TODO: Capture the job UUID/number from create_job output
        # For now, this is a placeholder
        call_command(*create_args)

        print("Note: Automatic run after create not yet fully implemented")
        print("Please use: ccp4i2 jobs run <project> <job>")
        return 0

    def _handle_status_shortcut(self, args):
        """Handle the 'status' shortcut."""
        if hasattr(args, 'project') and args.project:
            # Show single project status
            call_command('list_project', '--projectname', args.project)
        else:
            # Show all projects
            call_command('list_projects')
        return 0

    def _handle_i2run(self, args):
        """
        Handle legacy i2run command.

        This dispatches to CCP4i2RunnerDjango which:
        1. Creates or resolves project and job in database
        2. Parses task-specific parameters from command line
        3. Creates and configures plugin instance
        4. Saves parameters to database
        5. Executes the job via run_job_context_aware
        """
        import sys

        # Reconstruct sys.argv for i2run parser
        # Format: manage.py i2run task_name --project_name foo --param1 val1 --param2 val2
        i2run_argv = ['manage.py', 'i2run', args.task_name]

        # Add recognized i2run arguments
        i2run_argv.extend(['--project_name', args.project_name])

        if hasattr(args, 'project_path') and args.project_path:
            i2run_argv.extend(['--project_path', args.project_path])

        # Add plugin-specific parameters (unknown args from parse_known_args)
        if hasattr(args, 'unknown_args') and args.unknown_args:
            i2run_argv.extend(args.unknown_args)

        # Update sys.argv for i2run's parser (it reads sys.argv directly)
        original_sys_argv = sys.argv
        try:
            sys.argv = i2run_argv

            # Import and run the i2run command directly
            from ccp4x.db.management.commands.i2run import Command as I2RunCommand
            cmd = I2RunCommand()

            # Execute the command - handle() reads sys.argv directly
            cmd.handle()
            return 0

        except Exception as e:
            print(f"Error running i2run: {e}", file=sys.stderr)
            if hasattr(args, 'verbose') and args.verbose:
                import traceback
                traceback.print_exc()
            return 1
        finally:
            sys.argv = original_sys_argv

    def _is_uuid(self, value: str) -> bool:
        """Check if a string looks like a UUID."""
        import uuid
        try:
            uuid.UUID(value)
            return True
        except (ValueError, AttributeError):
            return False

    def _resolve_job_uuid(self, project_identifier: str, job_identifier: str):
        """Resolve project + job identifiers to job UUID."""
        from ccp4x.db.models import Job, Project
        import uuid

        try:
            # If job is already a UUID, use it directly
            if self._is_uuid(job_identifier):
                return uuid.UUID(job_identifier)

            # Otherwise, need to get project first, then job by number
            if self._is_uuid(project_identifier):
                project = Project.objects.get(uuid=uuid.UUID(project_identifier))
            else:
                project = Project.objects.get(name=project_identifier)

            job = Job.objects.get(number=job_identifier, project=project)
            return job.uuid

        except (Job.DoesNotExist, Project.DoesNotExist) as e:
            print(f"Error: Job or project not found: {e}", file=sys.stderr)
            return None


def main():
    """Main entry point."""
    cli = CLI()
    sys.exit(cli.run())


if __name__ == '__main__':
    main()
