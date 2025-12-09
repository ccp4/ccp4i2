Containers
==========

CCP4i2 plugins organize parameters into logical containers.

Standard Containers
-------------------

inputData
~~~~~~~~~

Input files and parameters:

.. code-block:: python

    plugin.container.inputData.XYZIN  # Input coordinate file
    plugin.container.inputData.HKLIN  # Input reflection file

outputData
~~~~~~~~~~

Output files:

.. code-block:: python

    plugin.container.outputData.XYZOUT  # Output coordinate file
    plugin.container.outputData.HKLOUT  # Output reflection file

controlParameters
~~~~~~~~~~~~~~~~~

Processing options:

.. code-block:: python

    plugin.container.controlParameters.NCYCLES  # Number of cycles
    plugin.container.controlParameters.WEIGHT_TYPE  # Weighting scheme

Accessing Parameters
--------------------

.. code-block:: python

    # Direct access
    value = plugin.container.inputData.NCYCLES.value

    # Using find()
    obj = plugin.container.find("inputData.NCYCLES")

    # Check if set
    if plugin.container.inputData.NCYCLES.isSet():
        ...
