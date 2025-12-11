Fundamental Types
=================

.. automodule:: core.base_object.fundamental_types
   :members:
   :undoc-members:
   :show-inheritance:

Fundamental types are the leaf nodes in the CCP4i2 data hierarchy, holding actual values.

Available Types
---------------

CString
~~~~~~~

String values with optional constraints.

.. code-block:: python

    title = CString(value="My Title")
    title.value = "New Title"

CInt
~~~~

Integer values with optional min/max constraints.

.. code-block:: python

    cycles = CInt(value=10)
    cycles.setQualifiers({'min': 1, 'max': 100})

CFloat
~~~~~~

Floating-point values with optional range constraints.

.. code-block:: python

    resolution = CFloat(value=2.5)
    resolution.setQualifiers({'min': 0.1, 'max': 100.0})

CBoolean
~~~~~~~~

Boolean values (True/False).

.. code-block:: python

    use_tls = CBoolean(value=True)

Common Interface
----------------

All fundamental types share:

- ``value`` property for get/set
- ``isSet()`` to check if value is defined
- ``unSet()`` to clear value
- ``setToDefault()`` to reset to default
