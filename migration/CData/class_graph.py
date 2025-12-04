"""
Dependency graph and topological sorting for CData classes.

Ensures classes are generated in correct order (dependencies before dependents).
"""

from typing import Dict, List, Tuple, Set
from pathlib import Path
from collections import defaultdict, deque


class ClassDependencyGraph:
    """Build and analyze class dependency graph for proper ordering."""

    def __init__(self, classes_data: Dict[str, dict], type_resolver):
        """
        Initialize with classes data and type resolver.

        Args:
            classes_data: Dict mapping class names to their metadata
            type_resolver: TypeResolver instance for resolving type references
        """
        self.classes_data = classes_data
        self.type_resolver = type_resolver
        self.graph = {}  # class_name → set of dependencies
        self.class_to_file = {}  # class_name → output filename

        self._build_class_to_file_map()
        self._build_dependency_graph()

    def _build_class_to_file_map(self):
        """Build mapping of class name to output file."""
        for class_name, class_info in self.classes_data.items():
            file_path = class_info.get('file_path', 'unknown.py')
            basename = Path(file_path).name
            self.class_to_file[class_name] = basename

    def _build_dependency_graph(self):
        """Build complete dependency graph for all classes."""
        self.graph = {name: set() for name in self.classes_data.keys()}

        for class_name, class_data in self.classes_data.items():
            deps = self.graph[class_name]

            # Add parent class as dependency
            parent = class_data.get('immediate_parent')
            if parent and parent in self.classes_data:
                # Only add if parent will be generated (not a base class)
                module, category = self.type_resolver.get_type_location(parent)
                if category == 'custom':
                    deps.add(parent)

            # Add attribute types as dependencies
            for attr_name, attr_info in class_data.get('CONTENTS', {}).items():
                attr_type_str = attr_info.get('class', '')
                attr_type = self.type_resolver.resolve_type_name(attr_type_str)

                # Only add if it's a custom class that will be generated
                if attr_type in self.classes_data:
                    module, category = self.type_resolver.get_type_location(attr_type)
                    if category == 'custom':
                        deps.add(attr_type)

    def topological_sort_all(self) -> List[str]:
        """
        Perform topological sort on all classes.

        Returns:
            List of class names in dependency order (dependencies first)

        Raises:
            ValueError: If circular dependency detected
        """
        # Kahn's algorithm
        in_degree = {name: 0 for name in self.graph}

        # Calculate in-degrees
        for name in self.graph:
            for dep in self.graph[name]:
                if dep in in_degree:  # Only count dependencies within our classes
                    in_degree[dep] += 1

        # Start with nodes that have no dependencies
        queue = deque([name for name, degree in in_degree.items() if degree == 0])
        sorted_classes = []

        while queue:
            # Sort queue for deterministic output
            current = queue.popleft()
            sorted_classes.append(current)

            # Reduce in-degree for dependent classes
            for name in self.graph:
                if current in self.graph[name]:
                    in_degree[name] -= 1
                    if in_degree[name] == 0:
                        queue.append(name)

        # Check for cycles
        if len(sorted_classes) != len(self.graph):
            remaining = set(self.graph.keys()) - set(sorted_classes)
            raise ValueError(f"Circular dependency detected among classes: {remaining}")

        # Reverse because we counted in-degrees backwards
        return list(reversed(sorted_classes))

    def get_sorted_classes_by_file(self) -> Dict[str, List[Tuple[str, dict]]]:
        """
        Get classes grouped by output file, sorted by dependencies.

        Returns:
            Dict mapping filename to list of (class_name, class_data) tuples,
            where classes within each file are in dependency order.
        """
        # Get global topological sort
        try:
            sorted_all = self.topological_sort_all()
        except ValueError as e:
            print(f"Warning: {e}")
            print("Falling back to best-effort sorting...")
            sorted_all = self._best_effort_sort()

        # Group by file while preserving order
        file_to_classes = defaultdict(list)
        for class_name in sorted_all:
            filename = self.class_to_file.get(class_name, 'unknown.py')
            class_data = self.classes_data[class_name]
            file_to_classes[filename].append((class_name, class_data))

        # Re-sort classes within each file to ensure parent classes come first
        for filename in file_to_classes:
            file_to_classes[filename] = self._sort_classes_in_file(file_to_classes[filename])

        return dict(file_to_classes)

    def _sort_classes_in_file(self, classes_in_file: List[Tuple[str, dict]]) -> List[Tuple[str, dict]]:
        """
        Sort classes within a single file so parents come before children.

        Args:
            classes_in_file: List of (class_name, class_data) tuples

        Returns:
            Sorted list where dependencies within the file come first
        """
        class_names = {name for name, _ in classes_in_file}
        class_map = {name: data for name, data in classes_in_file}

        # Build local dependency graph (only within this file)
        local_graph = {name: set() for name in class_names}
        for name in class_names:
            parent = class_map[name].get('immediate_parent')
            if parent and parent in class_names:
                local_graph[name].add(parent)

            # Also check attribute types
            for attr_info in class_map[name].get('CONTENTS', {}).values():
                attr_type = self.type_resolver.resolve_type_name(attr_info.get('class', ''))
                if attr_type in class_names:
                    local_graph[name].add(attr_type)

        # Topological sort within file
        # Calculate in-degree: how many dependencies each class has
        in_degree = {name: len(local_graph[name]) for name in local_graph}

        # Start with classes that have no dependencies
        queue = deque([name for name, degree in in_degree.items() if degree == 0])
        sorted_names = []

        while queue:
            # Sort queue for deterministic output
            current = queue.popleft()
            sorted_names.append(current)

            # Reduce in-degree for classes that depend on current
            for name in local_graph:
                if current in local_graph[name]:
                    in_degree[name] -= 1
                    if in_degree[name] == 0:
                        queue.append(name)

        # If there are remaining (circular deps within file), add them sorted
        remaining = set(class_names) - set(sorted_names)
        if remaining:
            print(f"    Warning: Circular dependencies within file: {remaining}")
            sorted_names.extend(sorted(remaining))

        return [(name, class_map[name]) for name in sorted_names]

    def _best_effort_sort(self) -> List[str]:
        """
        Best-effort sorting when circular dependencies exist.

        Uses dependency depth to approximate correct order.
        """
        def get_depth(name, visited=None):
            """Calculate maximum dependency depth for a class."""
            if visited is None:
                visited = set()

            if name in visited:
                return 0  # Circular reference, break cycle

            visited.add(name)
            max_depth = 0

            for dep in self.graph.get(name, set()):
                if dep in self.classes_data:
                    depth = get_depth(dep, visited.copy())
                    max_depth = max(max_depth, depth + 1)

            return max_depth

        # Sort by depth (deeper dependencies come first)
        classes_with_depth = [(name, get_depth(name)) for name in self.graph.keys()]
        classes_with_depth.sort(key=lambda x: (-x[1], x[0]))  # Desc by depth, then by name

        return [name for name, _ in classes_with_depth]

    def get_file_dependencies(self) -> Dict[str, Set[str]]:
        """
        Get file-level dependencies.

        Returns:
            Dict mapping filename to set of filenames it depends on
        """
        file_deps = defaultdict(set)

        for class_name, deps in self.graph.items():
            source_file = self.class_to_file.get(class_name, 'unknown.py')

            for dep in deps:
                dep_file = self.class_to_file.get(dep)
                if dep_file and dep_file != source_file:
                    file_deps[source_file].add(dep_file)

        return dict(file_deps)

    def get_class_dependencies(self, class_name: str) -> Set[str]:
        """
        Get direct dependencies for a class.

        Args:
            class_name: Name of class to query

        Returns:
            Set of class names this class depends on
        """
        return self.graph.get(class_name, set()).copy()

    def find_circular_dependencies(self) -> List[List[str]]:
        """
        Find all circular dependency chains.

        Returns:
            List of cycles, where each cycle is a list of class names
        """
        cycles = []
        visited = set()
        rec_stack = set()

        def dfs(node, path):
            visited.add(node)
            rec_stack.add(node)
            path.append(node)

            for neighbor in self.graph.get(node, set()):
                if neighbor not in self.classes_data:
                    continue

                if neighbor not in visited:
                    dfs(neighbor, path.copy())
                elif neighbor in rec_stack:
                    # Found a cycle
                    cycle_start = path.index(neighbor)
                    cycle = path[cycle_start:] + [neighbor]
                    cycles.append(cycle)

            rec_stack.remove(node)

        for node in self.graph:
            if node not in visited:
                dfs(node, [])

        return cycles

    def print_dependency_report(self):
        """Print dependency analysis report."""
        print("=" * 80)
        print("CLASS DEPENDENCY ANALYSIS")
        print("=" * 80)

        # File-level dependencies
        print("\nFile Dependencies:")
        print("-" * 80)
        file_deps = self.get_file_dependencies()
        for filename in sorted(file_deps.keys()):
            deps = file_deps[filename]
            if deps:
                print(f"\n{filename} depends on:")
                for dep in sorted(deps):
                    print(f"  - {dep}")

        # Check for circular dependencies
        print("\n" + "=" * 80)
        print("Circular Dependency Check:")
        print("-" * 80)
        cycles = self.find_circular_dependencies()
        if cycles:
            print(f"WARNING: Found {len(cycles)} circular dependency chain(s):")
            for i, cycle in enumerate(cycles, 1):
                print(f"\n  Cycle {i}: {' -> '.join(cycle)}")
        else:
            print("✓ No circular dependencies detected")

        # Class count by file
        print("\n" + "=" * 80)
        print("Classes per File:")
        print("-" * 80)
        file_counts = defaultdict(int)
        for class_name in self.classes_data:
            filename = self.class_to_file.get(class_name, 'unknown.py')
            file_counts[filename] += 1

        for filename in sorted(file_counts.keys()):
            print(f"  {filename}: {file_counts[filename]} classes")

        print("\n" + "=" * 80)
