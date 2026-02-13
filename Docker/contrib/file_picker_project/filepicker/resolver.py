"""
File resolver - resolve file templates to actual files.

Handles:
- Glob pattern expansion
- Symlink dereferencing
- Sibling file discovery (for grabbing related files from autoprocessing dirs)
"""

import os
import glob as globmodule
from typing import List, Dict, Any, Optional, Tuple, Set


# Pipeline preference order for autoprocessing fallback (higher = better)
PIPELINE_PREFERENCE = [
    ('xia2-dials', 3),
    ('xia2-3dii', 2),
    ('autoPROC_staraniso', 1),
    ('autoPROC', 0),
]


class ResolvedFile:
    """Represents a resolved file to be transferred."""

    def __init__(self,
                 source_path,      # type: str
                 output_path,      # type: str
                 original_pattern, # type: str
                 is_sibling=False  # type: bool
                 ):
        # type: (...) -> None
        self.source_path = source_path        # Actual file path (dereferenced)
        self.output_path = output_path        # Path in output archive/destination
        self.original_pattern = original_pattern  # Pattern that matched this file
        self.is_sibling = is_sibling          # True if grabbed as sibling file

    def __repr__(self):
        return "ResolvedFile({!r} -> {!r})".format(
            self.source_path, self.output_path)


class FileResolver:
    """
    Resolves file templates to actual files for a crystal directory.
    """

    def __init__(self, crystal_path, crystal_name):
        # type: (str, str) -> None
        """
        Initialize the resolver.

        Args:
            crystal_path: Full path to the crystal directory
            crystal_name: Name of the crystal directory (e.g., "CDK4CyclinD1-x0001")
        """
        self.crystal_path = os.path.abspath(crystal_path)
        self.crystal_name = crystal_name

    def _resolve_symlink(self, path):
        # type: (str) -> str
        """
        Fully resolve a symlink to its final target.

        Handles:
        - Direct symlinks
        - Relative symlinks
        - Chained symlinks
        """
        if not os.path.islink(path):
            return path

        # os.path.realpath resolves all symlinks
        return os.path.realpath(path)

    def _get_sibling_files(self, resolved_path, sibling_patterns):
        # type: (str, List[str]) -> List[Tuple[str, str]]
        """
        Get sibling files from the same directory as a resolved symlink target.

        Args:
            resolved_path: The resolved (dereferenced) path of the main file
            sibling_patterns: Glob patterns for sibling files to grab

        Returns:
            List of (source_path, relative_name) tuples
        """
        siblings = []
        parent_dir = os.path.dirname(resolved_path)

        for pattern in sibling_patterns:
            full_pattern = os.path.join(parent_dir, pattern)
            matches = globmodule.glob(full_pattern)

            for match in matches:
                if os.path.isfile(match):
                    rel_name = os.path.basename(match)
                    siblings.append((match, rel_name))

        return siblings

    def resolve_spec(self, spec):
        # type: (Dict[str, Any]) -> List[ResolvedFile]
        """
        Resolve a single file spec to actual files.

        Args:
            spec: A file spec dictionary from config

        Returns:
            List of ResolvedFile objects
        """
        pattern = spec['pattern']
        output_name = spec.get('output_name')
        dereference = spec.get('dereference', True)
        grab_siblings = spec.get('grab_siblings')
        required = spec.get('required', False)

        # Expand the glob pattern relative to crystal directory
        full_pattern = os.path.join(self.crystal_path, pattern)
        matches = globmodule.glob(full_pattern)

        if not matches and required:
            raise FileNotFoundError(
                "Required file not found: {} in {}".format(
                    pattern, self.crystal_name))

        results = []
        seen_sources = set()  # type: Set[str]

        for match_path in matches:
            if not os.path.exists(match_path):
                continue

            # Skip directories
            if os.path.isdir(match_path) and not os.path.islink(match_path):
                continue

            # Determine the actual source file
            if dereference and os.path.islink(match_path):
                source_path = self._resolve_symlink(match_path)
            else:
                source_path = match_path

            # Skip if we've already seen this source
            if source_path in seen_sources:
                continue
            seen_sources.add(source_path)

            # Verify the source exists
            if not os.path.isfile(source_path):
                if required:
                    raise FileNotFoundError(
                        "Symlink target not found: {} -> {}".format(
                            match_path, source_path))
                continue

            # Determine output path
            if output_name:
                out_path = os.path.join(self.crystal_name, output_name)
            else:
                # Preserve relative path structure
                rel_path = os.path.relpath(match_path, self.crystal_path)
                out_path = os.path.join(self.crystal_name, rel_path)

            results.append(ResolvedFile(
                source_path=source_path,
                output_path=out_path,
                original_pattern=pattern,
            ))

            # Handle sibling files
            if grab_siblings and os.path.islink(match_path):
                siblings = self._get_sibling_files(source_path, grab_siblings)
                for sib_source, sib_name in siblings:
                    if sib_source in seen_sources:
                        continue
                    seen_sources.add(sib_source)

                    # Put siblings in an 'autoprocessing' subdirectory
                    sib_out_path = os.path.join(
                        self.crystal_name, 'autoprocessing', sib_name)

                    results.append(ResolvedFile(
                        source_path=sib_source,
                        output_path=sib_out_path,
                        original_pattern=pattern,
                        is_sibling=True,
                    ))

        return results

    def _score_pipeline(self, dirname):
        # type: (str) -> int
        """Score an autoprocessing directory name by pipeline preference."""
        for pipeline_name, score in PIPELINE_PREFERENCE:
            if pipeline_name in dirname:
                return score
        return -1

    def _find_best_unmerged(self, seen_sources):
        # type: (Set[str]) -> List[ResolvedFile]
        """
        Fallback: scan autoprocessing/*/ for *_scaled_unmerged.mtz when
        grab_siblings didn't find one (e.g., main symlink points to autoPROC).

        Ranks by pipeline: xia2-dials > xia2-3dii > autoPROC.
        Within same pipeline, prefers most recent (by mtime).
        Also grabs xia2.mmcif.bz2 from the same directory if present.
        """
        autoprocessing_path = os.path.join(self.crystal_path, 'autoprocessing')
        if not os.path.isdir(autoprocessing_path):
            return []

        candidates = []  # type: List[Tuple[int, float, str, str, str]]
        for subdir in os.listdir(autoprocessing_path):
            subdir_path = os.path.join(autoprocessing_path, subdir)
            if not os.path.isdir(subdir_path):
                continue

            pipeline_score = self._score_pipeline(subdir)

            for entry in os.listdir(subdir_path):
                if entry.endswith('_scaled_unmerged.mtz'):
                    full_path = os.path.join(subdir_path, entry)
                    if full_path in seen_sources:
                        continue
                    if os.path.isfile(full_path):
                        mtime = os.path.getmtime(full_path)
                        candidates.append(
                            (pipeline_score, mtime, full_path, entry, subdir_path))

        if not candidates:
            return []

        candidates.sort(key=lambda x: (x[0], x[1]), reverse=True)
        best = candidates[0]
        _, _, source_path, filename, best_dir = best

        results = []
        results.append(ResolvedFile(
            source_path=source_path,
            output_path=os.path.join(self.crystal_name, 'autoprocessing', filename),
            original_pattern='autoprocessing fallback',
            is_sibling=True,
        ))

        # Also grab xia2.mmcif.bz2 from the same directory
        mmcif_path = os.path.join(best_dir, 'xia2.mmcif.bz2')
        if mmcif_path not in seen_sources and os.path.isfile(mmcif_path):
            results.append(ResolvedFile(
                source_path=mmcif_path,
                output_path=os.path.join(
                    self.crystal_name, 'autoprocessing', 'xia2.mmcif.bz2'),
                original_pattern='autoprocessing fallback',
                is_sibling=True,
            ))

        return results

    def resolve_template(self, template):
        # type: (List[Dict[str, Any]]) -> List[ResolvedFile]
        """
        Resolve a full file template to actual files.

        After resolving all specs, checks whether an unmerged MTZ was found
        via grab_siblings. If not (e.g., main symlink pointed to autoPROC),
        falls back to scanning autoprocessing directories directly, ranking
        by pipeline preference (xia2-dials > xia2-3dii > autoPROC).

        Args:
            template: List of file spec dictionaries

        Returns:
            List of ResolvedFile objects
        """
        all_results = []
        seen_outputs = set()  # type: Set[str]

        for spec in template:
            try:
                files = self.resolve_spec(spec)
                for f in files:
                    # Avoid duplicate outputs
                    if f.output_path not in seen_outputs:
                        all_results.append(f)
                        seen_outputs.add(f.output_path)
            except FileNotFoundError as e:
                # Re-raise if required, otherwise skip
                if spec.get('required', False):
                    raise
                # Could log warning here

        # Fallback: if no *_scaled_unmerged.mtz was found via grab_siblings,
        # scan autoprocessing directories directly
        has_unmerged = any(
            f.is_sibling and f.output_path.endswith('_scaled_unmerged.mtz')
            for f in all_results
        )
        if not has_unmerged:
            seen_sources = {f.source_path for f in all_results}
            fallback_files = self._find_best_unmerged(seen_sources)
            for f in fallback_files:
                if f.output_path not in seen_outputs:
                    all_results.append(f)
                    seen_outputs.add(f.output_path)

        return all_results


def resolve_crystals(selector, indices, template):
    # type: (Any, List[int], List[Dict[str, Any]]) -> Dict[int, List[ResolvedFile]]
    """
    Resolve files for multiple crystals.

    Args:
        selector: CrystalSelector instance
        indices: List of crystal indices to process
        template: File template to apply

    Returns:
        Dict mapping crystal index to list of ResolvedFile objects
    """
    results = {}

    for index in indices:
        crystal_path = selector.get_crystal_path(index)
        crystal_name = selector.get_crystal_dir_name(index)

        if crystal_path is None:
            continue

        resolver = FileResolver(crystal_path, crystal_name)
        try:
            files = resolver.resolve_template(template)
            results[index] = files
        except FileNotFoundError as e:
            # Could collect errors for reporting
            raise RuntimeError(
                "Failed to resolve files for crystal x{:04d}: {}".format(
                    index, str(e)))

    return results


def calculate_total_size(resolved_files):
    # type: (Dict[int, List[ResolvedFile]]) -> int
    """Calculate total size of all resolved files in bytes."""
    total = 0
    seen = set()  # type: Set[str]

    for index, files in resolved_files.items():
        for f in files:
            if f.source_path not in seen:
                seen.add(f.source_path)
                try:
                    total += os.path.getsize(f.source_path)
                except OSError:
                    pass

    return total


def format_size(size_bytes):
    # type: (int) -> str
    """Format byte size as human-readable string."""
    for unit in ['B', 'KB', 'MB', 'GB', 'TB']:
        if size_bytes < 1024:
            return "{:.1f} {}".format(size_bytes, unit)
        size_bytes /= 1024.0
    return "{:.1f} PB".format(size_bytes)
