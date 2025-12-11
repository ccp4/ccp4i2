CContainer
==========

.. automodule:: core.CCP4Container
   :members:
   :undoc-members:
   :show-inheritance:

CContainer is the base class for container objects that hold collections of CData children.

Key Features
------------

- Hierarchical organization of parameters
- Child object management
- Navigation via ``find()`` method
- XML serialization

Common Usage
------------

.. code-block:: python

    # Access child objects
    hklin = container.inputData.HKLIN

    # Navigate by path
    obj = container.find("inputData.HKLIN")

    # Iterate over children
    for child in container.children():
        print(child.objectName())
