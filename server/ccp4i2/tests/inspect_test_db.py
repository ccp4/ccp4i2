#!/usr/bin/env python
"""
Inspect database contents from a test run using raw SQLite queries.
"""
import sqlite3
import sys
from pathlib import Path

def inspect_database(db_path):
    """Inspect and display database contents."""
    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row
    cursor = conn.cursor()

    print("\n" + "="*120)
    print("DATABASE INSPECTION REPORT")
    print("="*120)

    # 1. PROJECT
    print("\n📁 PROJECT")
    print("-" * 120)
    cursor.execute("SELECT * FROM ccp4i2_project ORDER BY creation_time DESC LIMIT 1")
    project = cursor.fetchone()

    if not project:
        print("No projects found")
        return

    print(f"  UUID            : {project['uuid']}")
    print(f"  Name            : {project['name']}")
    print(f"  Directory       : {project['directory']}")
    print(f"  Created         : {project['creation_time']}")

    # 2. JOBS
    print("\n⚙️  JOBS")
    print("-" * 120)
    print(f"{'Job #':<8} {'Task':<18} {'Title':<30} {'Status':<12} {'Started':<20} {'Finished':<20}")
    print("-" * 120)
    cursor.execute("""
        SELECT * FROM ccp4i2_job WHERE project_id = ? ORDER BY start_time
    """, (project['id'],))

    for job in cursor.fetchall():
        print(f"{job['number']:<8} {job['task_name']:<18} {(job['title'] or 'N/A')[:29]:<30} "
              f"{job['status']:<12} {job['start_time'] or 'N/A':<20} {job['finish_time'] or 'N/A':<20}")

    # 3. FILES
    print("\n📄 FILES")
    print("-" * 120)
    print(f"{'Job':<6} {'Parameter':<18} {'Filename':<25} {'Type':<35} {'Size':<15}")
    print("-" * 120)

    cursor.execute("""
        SELECT f.*, j.number as job_number, ft.name as type_name
        FROM ccp4i2_file f
        JOIN ccp4i2_job j ON f.job_id = j.id
        LEFT JOIN ccp4i2_filetype ft ON f.type_id = ft.id
        WHERE j.project_id = ?
    """, (project['id'],))

    for file in cursor.fetchall():
        size = "N/A"
        if file['path']:
            fpath = Path(file['path'])
            if fpath.exists():
                fsize = fpath.stat().st_size
                if fsize > 1024 * 1024:
                    size = f"{fsize/(1024*1024):.2f} MB"
                elif fsize > 1024:
                    size = f"{fsize/1024:.2f} KB"
                else:
                    size = f"{fsize} bytes"

        print(f"{file['job_number']:<6} {(file['job_param_name'] or 'N/A')[:17]:<18} {file['name'][:24]:<25} "
              f"{(file['type_name'] or 'unknown')[:34]:<35} {size:<15}")

    # 4. FILE USES
    print("\n🔗 FILE USES (Input/Output)")
    print("-" * 120)
    print(f"{'Job':<6} {'Role':<12} {'Parameter':<20} {'Filename':<35}")
    print("-" * 120)

    cursor.execute("""
        SELECT fu.*, j.number as job_number, f.name as file_name
        FROM ccp4i2_fileuse fu
        JOIN ccp4i2_job j ON fu.job_id = j.id
        JOIN ccp4i2_file f ON fu.file_id = f.id
        WHERE j.project_id = ?
    """, (project['id'],))

    for use in cursor.fetchall():
        role_map = {0: 'Input', 1: 'Output', 2: 'Import'}
        role = role_map.get(use['role'], 'Unknown')
        print(f"{use['job_number']:<6} {role:<12} {(use['job_param_name'] or 'N/A')[:19]:<20} {use['file_name'][:34]:<35}")

    # 5. FILE TYPES
    print("\n🏷️  FILE TYPES")
    print("-" * 120)
    print(f"{'MIME Type':<50} {'Description':<65}")
    print("-" * 120)

    cursor.execute("""
        SELECT DISTINCT ft.*
        FROM ccp4i2_filetype ft
        JOIN ccp4i2_file f ON f.type_id = ft.id
        JOIN ccp4i2_job j ON f.job_id = j.id
        WHERE j.project_id = ?
    """, (project['id'],))

    for ft in cursor.fetchall():
        print(f"{ft['name'][:49]:<50} {(ft['description'] or 'N/A')[:64]:<65}")

    print("\n" + "="*120)
    print("END OF DATABASE INSPECTION")
    print("="*120 + "\n")

    # Show file tree
    project_dir = Path(project['directory'])
    if project_dir.exists():
        print("="*120)
        print("FILE SYSTEM HIERARCHY")
        print("="*120)
        print(f"\n📂 {project_dir}\n")

        def show_tree(path, prefix="", is_last=True):
            items = sorted(path.iterdir(), key=lambda p: (not p.is_dir(), p.name))
            for i, item in enumerate(items):
                is_last_item = (i == len(items) - 1)
                connector = "└── " if is_last_item else "├── "
                new_prefix = prefix + ("    " if is_last_item else "│   ")

                if item.is_file():
                    size = item.stat().st_size
                    if size > 1024 * 1024:
                        size_str = f"{size / (1024*1024):.2f} MB"
                    elif size > 1024:
                        size_str = f"{size / 1024:.2f} KB"
                    else:
                        size_str = f"{size} bytes"
                    print(f"{prefix}{connector}{item.name} ({size_str})")
                else:
                    print(f"{prefix}{connector}{item.name}/")
                    if item.is_dir():
                        show_tree(item, new_prefix, is_last_item)

        show_tree(project_dir)
        print("\n" + "="*120 + "\n")

    conn.close()


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python inspect_test_db.py <database_file>")
        sys.exit(1)

    db_path = sys.argv[1]
    if not Path(db_path).exists():
        print(f"Error: Database file not found: {db_path}")
        sys.exit(1)

    inspect_database(db_path)
