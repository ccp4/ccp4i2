CLI Overview
============

The CCP4i2 CLI provides a modern command-line interface for crystallographic computing.

Command Structure
-----------------

.. code-block:: bash

    ccp4i2 <resource> <action> [<identifier>] [options]

Resources
~~~~~~~~~

- ``projects`` - Manage projects
- ``jobs`` - Manage jobs within projects
- ``files`` - File operations
- ``plugins`` - Plugin/task discovery

Quick Examples
--------------

.. code-block:: bash

    # List projects
    ccp4i2 projects list

    # Create a project
    ccp4i2 projects create my_project

    # Create and run a job
    ccp4i2 run my_project refmac5

    # Check job status
    ccp4i2 jobs show my_project 1

For detailed documentation, see ``mddocs/cli/CLI.md``.
