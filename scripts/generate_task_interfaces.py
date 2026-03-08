#!/usr/bin/env python3
"""
Generate React task interface .tsx files from legacy Qt CTaskWidget GUI Python files.

Parses drawContents() method from CTaskWidget subclasses and generates
corresponding React/Material-UI components.

Usage:
    python3 scripts/generate_task_interfaces.py [--task TASKNAME] [--all-missing] [--dry-run]
"""

import ast
import json
import os
import re
import subprocess
import sys
import textwrap
from dataclasses import dataclass, field
from typing import Optional


@dataclass
class Folder:
    title: str
    is_input_data: bool = False
    elements: list = field(default_factory=list)


@dataclass
class Element:
    kind: str  # 'widget', 'label', 'subtitle', 'advice', 'inline'
    item_name: Optional[str] = None
    label: Optional[str] = None
    gui_label: Optional[str] = None
    toggle: Optional[dict] = None  # {param, values}
    children: list = field(default_factory=list)
    browse_db: bool = False


@dataclass
class SubFrame:
    frame: bool = False
    toggle: Optional[dict] = None
    elements: list = field(default_factory=list)


def get_qt_source(git_path: str) -> str:
    """Get file content from main branch."""
    result = subprocess.run(
        ["git", "show", f"main:{git_path}"],
        capture_output=True, text=True
    )
    if result.returncode != 0:
        raise FileNotFoundError(f"Cannot read main:{git_path}")
    return result.stdout


def find_class_and_taskname(source: str) -> tuple[str, str]:
    """Find CTaskWidget subclass and extract TASKNAME."""
    # Find TASKNAME
    m = re.search(r"TASKNAME\s*=\s*['\"](\w+)['\"]", source)
    if not m:
        raise ValueError("No TASKNAME found")
    taskname = m.group(1)

    # Find class name
    m = re.search(r"class\s+(\w+)\s*\(", source)
    if not m:
        raise ValueError("No class found")
    classname = m.group(1)

    return classname, taskname


def parse_draw_contents(source: str) -> list[Folder]:
    """Parse the drawContents method to extract GUI structure."""
    # Extract drawContents method body
    m = re.search(r'def drawContents\(self\):\s*\n(.*?)(?=\n    def |\n  def |\nclass |\Z)',
                  source, re.DOTALL)
    if not m:
        raise ValueError("No drawContents method found")

    body = m.group(1)
    lines = body.split('\n')

    folders = []
    current_folder = None
    current_subframe = None
    toggle_stack = []

    for line in lines:
        stripped = line.strip()
        if not stripped or stripped.startswith('#'):
            continue

        # openFolder
        folder_match = re.search(
            r"self\.openFolder\s*\((.*?)\)", stripped)
        if folder_match:
            args = folder_match.group(1)
            title = "Input Data"
            is_input = False
            t_match = re.search(r"title\s*=\s*['\"](.+?)['\"]", args)
            if t_match:
                title = t_match.group(1)
            if "folderFunction='inputData'" in args or 'folderFunction="inputData"' in args:
                is_input = True
                if not t_match:
                    title = "Input Data"
            current_folder = Folder(title=title, is_input_data=is_input)
            folders.append(current_folder)
            current_subframe = None
            continue

        # openSubFrame
        subframe_match = re.search(
            r"self\.openSubFrame\s*\((.*?)\)", stripped)
        if subframe_match:
            args = subframe_match.group(1)
            toggle = parse_toggle_arg(args)
            current_subframe = SubFrame(frame=True, toggle=toggle)
            if current_folder:
                current_folder.elements.append(current_subframe)
            continue

        # closeSubFrame
        if 'self.closeSubFrame' in stripped:
            current_subframe = None
            continue

        # createLine
        line_match = re.search(r"self\.createLine\s*\(\s*\[(.*?)\]\s*(?:,\s*(.*?))?\s*\)",
                               stripped, re.DOTALL)
        if line_match:
            items_str = line_match.group(1)
            extra_args = line_match.group(2) or ""
            elements = parse_create_line(items_str, extra_args)
            target = current_subframe.elements if current_subframe else (
                current_folder.elements if current_folder else [])
            target.extend(elements)
            continue

        # setMenuText - we note it but don't need it for the generic interface
        # since the backend handles enumerators
        if 'self.setMenuText' in stripped:
            continue

        # setProgramHelpFile - skip
        if 'self.setProgramHelpFile' in stripped:
            continue

        # setToolTip - skip
        if 'self.setToolTip' in stripped:
            continue

        # Signal connections - skip
        if '.dataChanged.connect' in stripped or '.clicked.connect' in stripped:
            continue

        # getWidget calls - skip
        if 'self.getWidget' in stripped:
            continue

        # container access - skip
        if 'self.container.' in stripped:
            continue

        # Method calls on self that aren't GUI building - skip
        if re.match(r'self\.\w+\(\)', stripped):
            continue

    # If no folders were created, make a default one
    if not folders:
        folders.append(Folder(title="Input Data", is_input_data=True))

    return folders


def parse_toggle_arg(args_str: str) -> Optional[dict]:
    """Parse toggle argument from openSubFrame or createLine."""
    toggle_match = re.search(
        r"toggle\s*=\s*\[([^\]]+)\]", args_str)
    if not toggle_match:
        return None
    toggle_content = toggle_match.group(1)
    # Parse ['PARAM', 'open', ['VALUE1', 'VALUE2']]
    parts = re.findall(r"'([^']+)'", toggle_content)
    if len(parts) >= 3:
        param = parts[0]
        condition = parts[1]  # 'open' or 'close'
        values = parts[2:]
        return {"param": param, "condition": condition, "values": values}
    return None


def parse_create_line(items_str: str, extra_args: str = "") -> list[Element]:
    """Parse a createLine call to extract elements."""
    elements = []
    # Tokenize the items string
    tokens = re.findall(r"'([^']*)'|\"([^\"]*)\"", items_str)
    tokens = [t[0] or t[1] for t in tokens]

    toggle = parse_toggle_arg(extra_args) if extra_args else None
    # Also check for toggle in the items themselves (some put it inline)
    if not toggle:
        toggle = parse_toggle_arg(items_str)

    i = 0
    while i < len(tokens):
        token = tokens[i]

        if token == 'widget':
            # Next token might be a modifier like '-guiMode', '-browseDb', '-title'
            i += 1
            browse_db = False
            gui_label = None
            while i < len(tokens) and tokens[i].startswith('-'):
                modifier = tokens[i]
                i += 1
                if modifier == '-browseDb' and i < len(tokens):
                    browse_db = True
                    # Skip the True/False value
                    if i < len(tokens) and tokens[i] in ('True', 'False'):
                        i += 1
                elif modifier == '-guiMode' and i < len(tokens):
                    # Skip the mode value (combo, multiLine, multiLineRadio, etc.)
                    i += 1
                elif modifier == '-title' and i < len(tokens):
                    gui_label = tokens[i]
                    i += 1
            if i < len(tokens):
                item_name = tokens[i]
                elem = Element(kind='widget', item_name=item_name,
                               toggle=toggle, browse_db=browse_db,
                               gui_label=gui_label)
                elements.append(elem)
                i += 1
            continue

        elif token == 'label':
            i += 1
            if i < len(tokens):
                label_text = tokens[i]
                # Check if next token is 'widget'
                if i + 1 < len(tokens) and tokens[i + 1] == 'widget':
                    # This is an inline label+widget pattern
                    i += 2  # skip 'widget'
                    # Handle modifiers
                    while i < len(tokens) and tokens[i].startswith('-'):
                        i += 1
                        if i < len(tokens):
                            i += 1
                    if i < len(tokens):
                        item_name = tokens[i]
                        elem = Element(kind='widget', item_name=item_name,
                                       gui_label=label_text.strip(),
                                       toggle=toggle)
                        elements.append(elem)
                        i += 1
                    continue
                else:
                    # Standalone label - usually part of inline pattern
                    # We'll capture it but it often pairs with widget
                    i += 1
                    continue

        elif token == 'subtitle':
            i += 1
            if i < len(tokens):
                elements.append(Element(kind='subtitle', label=tokens[i]))
                i += 1
            continue

        elif token == 'advice':
            i += 1
            if i < len(tokens):
                elements.append(Element(kind='advice', label=tokens[i]))
                i += 1
            continue

        elif token == 'tip':
            i += 1
            if i < len(tokens):
                # Tips are tooltips, skip them - next token is usually 'widget'
                i += 1
            continue

        elif token == 'stretch':
            i += 1
            continue

        else:
            i += 1

    return elements


def collect_widgets(folders: list[Folder]) -> list[dict]:
    """Flatten folder structure into a list of renderable items."""
    result = []

    for folder in folders:
        folder_items = []
        for elem in folder.elements:
            if isinstance(elem, SubFrame):
                for sub_elem in elem.elements:
                    if isinstance(sub_elem, Element) and sub_elem.kind == 'widget':
                        folder_items.append({
                            'type': 'widget',
                            'name': sub_elem.item_name,
                            'guiLabel': sub_elem.gui_label,
                            'toggle': sub_elem.toggle or elem.toggle,
                        })
                    elif isinstance(sub_elem, Element) and sub_elem.kind == 'subtitle':
                        folder_items.append({
                            'type': 'subtitle',
                            'label': sub_elem.label,
                        })
            elif isinstance(elem, Element):
                if elem.kind == 'widget':
                    folder_items.append({
                        'type': 'widget',
                        'name': elem.item_name,
                        'guiLabel': elem.gui_label,
                        'toggle': elem.toggle,
                    })
                elif elem.kind == 'subtitle':
                    folder_items.append({
                        'type': 'subtitle',
                        'label': elem.label,
                    })

        if folder_items:
            result.append({
                'folder_title': folder.title,
                'items': folder_items,
            })

    return result


def generate_tsx(task_name: str, folders_data: list[dict],
                 toggle_params: set[str]) -> str:
    """Generate the .tsx file content."""
    # Determine imports needed
    needs_bool_toggle = len(toggle_params) > 0
    needs_linear_progress = True  # always include for loading state

    imports = ['import { LinearProgress, Paper } from "@mui/material";']
    imports.append('import { CCP4i2TaskInterfaceProps } from "./task-container";')
    imports.append('import { CCP4i2TaskElement } from "../task-elements/task-element";')
    imports.append('import { CCP4i2ContainerElement } from "../task-elements/ccontainer";')
    imports.append('import { useJob } from "../../../utils";')

    if needs_bool_toggle:
        imports.append('import { useBoolToggle } from "../task-elements/shared-hooks";')

    # Build component body
    hooks = ['  const { useTaskItem, container } = useJob(props.job.id);']

    for param in sorted(toggle_params):
        hooks.append(f'  const {_camel(param)} = useBoolToggle(useTaskItem, "{param}");')

    hooks.append('')
    hooks.append('  if (!container) return <LinearProgress />;')

    # Build JSX
    jsx_parts = []
    for folder_data in folders_data:
        title = folder_data['folder_title']
        items = folder_data['items']

        folder_jsx = []
        folder_jsx.append(f'      {{/* {title} */}}')
        folder_jsx.append(f'      <CCP4i2ContainerElement')
        folder_jsx.append(f'        {{...props}}')
        folder_jsx.append(f'        itemName=""')
        folder_jsx.append(f'        qualifiers={{{{ guiLabel: "{_escape_jsx(title)}" }}}}')
        folder_jsx.append(f'        containerHint="FolderLevel"')
        folder_jsx.append(f'      >')

        for item in items:
            if item['type'] == 'widget':
                widget_jsx = _render_widget(item, toggle_params)
                folder_jsx.extend(widget_jsx)

        folder_jsx.append(f'      </CCP4i2ContainerElement>')
        jsx_parts.append('\n'.join(folder_jsx))

    jsx_body = '\n\n'.join(jsx_parts)

    # Assemble full component
    component = f"""{chr(10).join(imports)}

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {{
{chr(10).join(hooks)}

  return (
    <Paper sx={{{{ display: "flex", flexDirection: "column", gap: 1, p: 1 }}}}>
{jsx_body}
    </Paper>
  );
}};

export default TaskInterface;
"""
    return component


def _render_widget(item: dict, toggle_params: set[str]) -> list[str]:
    """Render a single widget element as JSX lines."""
    name = item['name']
    gui_label = item.get('guiLabel')
    toggle = item.get('toggle')

    lines = []
    indent = '        '

    # Handle toggle wrapping
    if toggle:
        param = toggle['param']
        condition = toggle.get('condition', 'open')
        values = toggle.get('values', [])
        if param in toggle_params:
            camel = _camel(param)
            if condition == 'open' and len(values) > 0:
                # For boolean toggles
                lines.append(f'{indent}{{{camel}.value && (')
                indent = '          '
            elif condition == 'close':
                lines.append(f'{indent}{{!{camel}.value && (')
                indent = '          '

    widget_line = f'{indent}<CCP4i2TaskElement itemName="{name}" {{...props}}'
    if gui_label:
        widget_line += f' qualifiers={{{{ guiLabel: "{_escape_jsx(gui_label)}" }}}}'
    widget_line += ' />'
    lines.append(widget_line)

    # Close toggle wrapper
    if toggle and toggle.get('param') in toggle_params:
        lines.append(f'        )}}'  )

    return lines


def _camel(name: str) -> str:
    """Convert PARAM_NAME to a camelCase variable name for hooks."""
    parts = name.lower().split('_')
    return parts[0] + ''.join(p.capitalize() for p in parts[1:])


def _escape_jsx(s: str) -> str:
    """Escape string for use in JSX."""
    return s.replace('"', '\\"').replace("'", "\\'")


def find_toggle_params(folders_data: list[dict]) -> set[str]:
    """Find all parameters used in toggle conditions."""
    params = set()
    for folder in folders_data:
        for item in folder['items']:
            toggle = item.get('toggle')
            if toggle:
                params.add(toggle['param'])
    return params


def process_task(git_path: str, task_name_override: str = None) -> dict:
    """Process a single Qt GUI file and generate a .tsx file."""
    source = get_qt_source(git_path)
    classname, task_name = find_class_and_taskname(source)
    if task_name_override:
        task_name = task_name_override

    folders = parse_draw_contents(source)
    folders_data = collect_widgets(folders)
    toggle_params = find_toggle_params(folders_data)

    # Only include boolean toggles (not menu enumerators)
    # We can't easily distinguish, so include all for now
    # The useBoolToggle hook handles non-booleans gracefully

    tsx_content = generate_tsx(task_name, folders_data, toggle_params)

    return {
        'task_name': task_name,
        'source': git_path,
        'tsx': tsx_content,
        'element_count': sum(
            len([i for i in f['items'] if i['type'] == 'widget'])
            for f in folders_data
        ),
    }


# Map of missing tasks to their Qt GUI file paths on main branch
MISSING_TASKS = {
    'MakeMonster': 'wrappers/MakeMonster/script/MakeMonster_gui.py',
    'ProvideTLS': 'wrappers/ProvideTLS/script/ProvideTLS_GUI.py',
    'arp_warp_classic': 'wrappers/arp_warp_classic/script/CArpWarpClassic.py',
    'buccaneer_build_refine_mr': 'pipelines/buccaneer_build_refine_mr/script/CTaskbuccaneer_build_refine_mr.py',
    'buster': 'wrappers/buster/script/CTaskbuster.py',
    'chainsaw': 'wrappers/chainsaw/script/CTaskChainsaw.py',
    'cif2mtz': 'wrappers/cif2mtz/script/CTaskCif2mtz.py',
    'dials_image': 'wrappers/dials_image/script/CTaskDials_image.py',
    'dials_rlattice': 'wrappers/dials_rlattice/script/CTaskDials_rlattice.py',
    'dui': 'wrappers/dui/script/CTaskdui.py',
    'edstats': 'wrappers/edstats/script/CTaskEdstats.py',
    'findmyseq': 'wrappers/findmyseq/script/CTaskfindmyseq.py',
    'imosflm': 'wrappers/imosflm/script/CTaskimosflm.py',
    'import_serial': 'wrappers/import_serial/script/CTaskimport_serial.py',
    'import_serial_pipe': 'pipelines/import_serial_pipe/script/CTaskimport_serial_pipe.py',
    'import_xia2': 'tasks/import_xia2/CTaskimport_xia2.py',
    'matthews': 'wrappers/matthews/script/matthews.py',
    'molrep_den': 'wrappers/molrep_den/script/Cmolrep_den.py',
    'mrparse': 'wrappers/mrparse/script/CTaskmrparse.py',
    'mtzutils': 'wrappers/mtzutils/script/CTaskMtzutils.py',
    'nautilus_build_refine': 'pipelines/nautilus_build_refine/script/nautilus_build_refine_gui.py',
    'newProject_fromMerged': 'tasks/newProject_fromMerged/CTaskNewProject_fromMerged.py',
    'pairef': 'wrappers/pairef/script/CTaskpairef.py',
    'pdbset_ui': 'wrappers/pdbset_ui/script/CTaskpdbset_ui.py',
    'phaser_mr': 'tasks/phaser_mr/CTaskPhaser_mr.py',
    'phaser_singleMR': 'wrappers/phaser_singleMR/script/CTaskPhaserSingleMR.py',
    'prosmart': 'wrappers/prosmart/script/CTaskProsmart.py',
    'pyphaser_mr': 'wrappers/pyphaser_mr/script/CTaskPyphaser_mr.py',
    'scaleit': 'wrappers/scaleit/script/CTaskscaleit.py',
    'sculptor': 'wrappers/sculptor/script/CTaskSculptor.py',
    'slicendice': 'wrappers/slicendice/script/Cslicendice.py',
    'tableone': 'pipelines/tableone/script/CTasktableone.py',
    'unique': 'tasks/unique/CTaskUnique.py',
}


def main():
    import argparse
    parser = argparse.ArgumentParser(description='Generate React task interfaces from Qt GUI files')
    parser.add_argument('--task', help='Generate for a specific task name')
    parser.add_argument('--all-missing', action='store_true', help='Generate all missing tasks')
    parser.add_argument('--dry-run', action='store_true', help='Print output without writing files')
    parser.add_argument('--output-dir', default='client/renderer/components/task/task-interfaces',
                        help='Output directory for .tsx files')
    args = parser.parse_args()

    tasks = {}
    if args.task:
        if args.task in MISSING_TASKS:
            tasks[args.task] = MISSING_TASKS[args.task]
        else:
            print(f"Unknown task: {args.task}")
            sys.exit(1)
    elif args.all_missing:
        tasks = MISSING_TASKS
    else:
        parser.print_help()
        sys.exit(1)

    results = {'success': [], 'failed': []}

    for task_name, git_path in sorted(tasks.items()):
        try:
            result = process_task(git_path, task_name_override=task_name)
            tsx_content = result['tsx']

            output_file = os.path.join(args.output_dir, f"{task_name}.tsx")

            if args.dry_run:
                print(f"\n{'='*60}")
                print(f"=== {task_name} ({result['element_count']} elements) ===")
                print(f"{'='*60}")
                print(tsx_content)
            else:
                with open(output_file, 'w') as f:
                    f.write(tsx_content)
                print(f"  OK: {task_name} ({result['element_count']} elements) -> {output_file}")

            results['success'].append({
                'task_name': task_name,
                'source': git_path,
                'output': f"{task_name}.tsx",
                'element_count': result['element_count'],
            })

        except Exception as e:
            print(f"FAIL: {task_name} - {e}")
            results['failed'].append({
                'task_name': task_name,
                'source': git_path,
                'reason': str(e),
            })

    print(f"\nResults: {len(results['success'])} succeeded, {len(results['failed'])} failed")

    if results['failed']:
        print("\nFailed tasks:")
        for f in results['failed']:
            print(f"  {f['task_name']}: {f['reason']}")

    return results


if __name__ == '__main__':
    main()
