Validation Patterns
===================

Content-aware validation using ``validity()`` method overrides.

Overview
--------

Override ``validity()`` to implement custom validation logic:

.. code-block:: python

    def validity(self):
        # Adjust qualifiers based on mode
        if self.container.inputData.MODE == 'A':
            self.container.inputData.FIELD_B.setQualifier('allowUndefined', True)

        # Call parent validation
        return super().validity()

Common Patterns
---------------

1. **Conditional Requirements** - Make fields optional based on mode selection
2. **Filtering Expected Errors** - Remove errors for programmatically-populated fields
3. **Embedded Wrapper Adjustments** - Adjust child wrapper qualifiers

See ``mddocs/pipeline/VALIDITY_PATTERNS.md`` for comprehensive documentation.
