Signal System
=============

.. automodule:: core.base_object.signal_system
   :members:
   :undoc-members:
   :show-inheritance:

The signal system provides Qt-compatible signal/slot communication without Qt dependencies.

Key Classes
-----------

Signal
~~~~~~

Emits notifications to connected callbacks.

.. code-block:: python

    from core.base_object.signal_system import Signal

    # Create signal
    finished = Signal("finished", dict)

    # Connect callback
    finished.connect(my_callback)

    # Emit signal
    finished.emit({'status': 0})

SignalManager
~~~~~~~~~~~~~

Manages signals for an object.

.. code-block:: python

    class MyClass:
        def __init__(self):
            self._signal_manager = SignalManager(self)
            self.finished = self._signal_manager.create_signal("finished", dict)

Usage in Pipelines
------------------

.. code-block:: python

    class MyPipeline(CPluginScript):
        def process(self):
            # Create sub-plugin
            self.sub = self.makePluginObject('refmac5')

            # Connect to finished signal
            self.connectSignal(self.sub, 'finished', self.on_sub_finished)

            # Run sub-plugin
            self.sub.process()

        def on_sub_finished(self, status_dict):
            if status_dict['finishStatus'] == self.SUCCEEDED:
                self.reportStatus(self.SUCCEEDED)
            else:
                self.reportStatus(self.FAILED)
