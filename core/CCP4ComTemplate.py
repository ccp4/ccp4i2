"""
Modern CCP4ComTemplate - Qt-free template expansion for command lines and scripts.

This is a clean reimplementation of the legacy CCP4ComTemplate.py, designed to:
- Work without Qt dependencies
- Support template variable substitution ($VAR syntax)
- Handle dotted paths for nested attributes ($HKLIN.fullPath)
- Integrate with modern CData/HierarchicalObject system

Template Syntax:
- $VARIABLE - Expands to str(container.find("VARIABLE"))
- $VARIABLE.attribute - Expands to str(container.find("VARIABLE").attribute)
- Lines starting with # are comments (ignored)
"""

from typing import Tuple, Optional
import re


class CComTemplate:
    """Template processor for CCP4 command line and script generation."""

    def __init__(self, parent=None, template: str = None):
        """
        Initialize template processor.

        Args:
            parent: Parent object (for compatibility, not used)
            template: Template string with $VARIABLE substitutions
        """
        self.parent = parent
        self.template = template if template else ""

    def makeComScript(self, container) -> Tuple[str, 'CErrorReport']:
        """
        Process template and substitute variables from container.

        This is the main entry point called by CPluginScript.makeCommandAndScript().

        Args:
            container: Container object with find() method for variable lookup

        Returns:
            Tuple of (expanded_text, error_report)

        Example:
            >>> template = CComTemplate(template='1 HKLIN $HKLIN\\n1 END')
            >>> text, err = template.makeComScript(container)
            >>> # text might be: '1 HKLIN /path/to/file.mtz\\n1 END'
        """
        from core.base_object.error_reporting import CErrorReport

        error = CErrorReport()

        if not self.template:
            return "", error

        # Split into lines and process each
        lines = self.template.split('\n')
        output_lines = []

        for line in lines:
            # Skip empty lines and comments
            stripped = line.strip()
            if not stripped or stripped.startswith('#'):
                continue

            # Expand variables in this line
            try:
                expanded = self._expand_line(line, container)
                if expanded:  # Only add non-empty lines
                    output_lines.append(expanded)
            except Exception as e:
                error.append(
                    klass=self.__class__.__name__,
                    code=101,
                    details=f"Error expanding template line '{line}': {e}"
                )
                # Continue processing other lines

        # Join with spaces (for command line) or newlines (for scripts)
        # Command line templates typically on one line, scripts on multiple
        if len(output_lines) == 1:
            # Single line - return as is (command line case)
            result = output_lines[0]
        else:
            # Multiple lines - join with newlines (script case)
            result = '\n'.join(output_lines)

        return result, error

    def _expand_line(self, line: str, container) -> str:
        """
        Expand all $VARIABLE references in a line.

        Args:
            line: Template line with $VARIABLE references
            container: Container for variable lookup

        Returns:
            Expanded line with variables substituted

        Raises:
            ValueError: If a variable cannot be found or converted

        Notes:
            Legacy CCP4 templates use numeric prefixes like "1 HKLIN $HKLIN"
            where the leading number is a priority/line number indicator.
            We strip this prefix to get the actual command content.
        """
        # Strip leading numeric prefix (e.g., "1 HKLIN" -> "HKLIN")
        # This matches lines like: "1 HKLIN $HKLIN" or "2 END"
        line = line.strip()
        if line and line[0].isdigit():
            # Find where the number ends and content begins
            space_idx = line.find(' ')
            if space_idx > 0:
                line = line[space_idx+1:].strip()

        # Find all $VARIABLE or $VARIABLE.path patterns
        # Pattern matches: $word or $word.word.word
        pattern = r'\$([A-Za-z_][A-Za-z0-9_]*(?:\.[A-Za-z_][A-Za-z0-9_]*)*)'

        def replace_variable(match):
            var_path = match.group(1)  # e.g., "HKLIN" or "HKLIN.fullPath"
            return self._get_value(var_path, container)

        expanded = re.sub(pattern, replace_variable, line)
        return expanded

    def _get_value(self, var_path: str, container) -> str:
        """
        Get value for a variable path from container.

        Args:
            var_path: Variable path like "HKLIN" or "HKLIN.fullPath"
            container: Container for variable lookup

        Returns:
            String representation of the value

        Raises:
            ValueError: If variable not found or cannot be converted to string

        Example:
            >>> self._get_value("HKLIN", container)
            '/path/to/file.mtz'
            >>> self._get_value("CELL.a", container)
            '42.7'
        """
        parts = var_path.split('.')
        obj_name = parts[0]
        attr_path = parts[1:] if len(parts) > 1 else []

        # Look up object in container
        if not hasattr(container, 'find'):
            raise ValueError(f"Container does not have find() method (has type {type(container).__name__})")

        try:
            obj = container.find(obj_name)
        except Exception as e:
            raise ValueError(f"Error calling find('{obj_name}'): {e}")

        if obj is None:
            raise ValueError(f"Variable '{obj_name}' not found in container")

        # Navigate attribute path if present
        current = obj
        for attr in attr_path:
            if hasattr(current, attr):
                current = getattr(current, attr)
            elif hasattr(current, 'get') and callable(current.get):
                # Try CData.get() if available
                data = current.get()
                if isinstance(data, dict) and attr in data:
                    current = data[attr]
                else:
                    raise ValueError(f"Attribute '{attr}' not found on '{obj_name}'")
            else:
                raise ValueError(f"Attribute '{attr}' not found on '{obj_name}'")

        # Convert to string
        try:
            # Special handling for CDataFile: use str(current) to get full path
            # CDataFile.__str__() returns getFullPath() which computes the absolute path
            # from project/relPath/baseName using dbHandler
            from core.base_object.cdata_file import CDataFile
            if isinstance(current, CDataFile):
                # Temporarily set the file's plugin reference so it can find dbHandler
                # The template's parent is the plugin that owns the container
                if hasattr(self, 'parent') and self.parent is not None:
                    # Store plugin reference on the file object temporarily
                    current._temp_plugin_ref = self.parent

                try:
                    full_path = str(current)
                    if not full_path or len(full_path.strip()) == 0:
                        raise ValueError(f"CDataFile '{var_path}' has no path set")
                    return full_path
                finally:
                    # Clean up temporary reference
                    if hasattr(current, '_temp_plugin_ref'):
                        delattr(current, '_temp_plugin_ref')
            # Handle other CData objects with .value attribute
            elif hasattr(current, 'value'):
                value = current.value
                if value is None:
                    raise ValueError(f"Variable '{var_path}' has no value set")
                return str(value)
            else:
                # Direct string conversion
                return str(current)
        except Exception as e:
            raise ValueError(f"Cannot convert '{var_path}' to string: {e}")


# Legacy compatibility - some old code might expect these classes
class CComTemplateElement:
    """Legacy compatibility stub."""
    pass


class CComTemplateLine:
    """Legacy compatibility stub."""
    pass


class CComTemplateIf:
    """Legacy compatibility stub."""
    pass


class CComTemplateLoop:
    """Legacy compatibility stub."""
    pass


class CComTemplateCase:
    """Legacy compatibility stub."""
    pass
