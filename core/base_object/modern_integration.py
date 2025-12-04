"""
Integration example showing how to replace Qt-based CCP4i2 classes
with modern Python equivalents.

This demonstrates practical replacements for:
- CPluginScript -> ModernPluginScript
- Qt signal/slots -> Signal system
- Qt parent/child -> HierarchicalObject
- QEventLoop -> EventLoop
- Database API -> SignalEmitter integration
"""

import asyncio
import logging
from pathlib import Path
from typing import Any, Dict, Optional

# Import our new base systems
try:
    from .hierarchy_system import HierarchicalObject, DataContainer
    from .event_system import Application
except ImportError:
    from hierarchy_system import HierarchicalObject, DataContainer
    from event_system import Application

logger = logging.getLogger(__name__)


class ModernDbApi(HierarchicalObject):
    """
    Modern replacement for CCP4i2DjangoDbApi with proper signal system.

    Instead of FakeSignal, use real signals for database events.
    """

    def __init__(self, parent: Optional[HierarchicalObject] = None):
        super().__init__(parent, "ModernDbApi")

        # Replace Qt signals with our signal system
        self.project_reset = self.create_signal("project_reset", dict)
        self.job_status_changed = self.create_signal("job_status_changed", dict)
        self.file_imported = self.create_signal("file_imported", dict)
        self.data_changed = self.create_signal("data_changed", dict)

    def emit_project_reset(self, project_data: Dict[str, Any]):
        """Emit project reset signal with proper data."""
        logger.info("Project reset: %s", project_data.get("name", "Unknown"))
        self.project_reset.emit(project_data)

    def emit_job_status_change(self, job_id: str, old_status: str, new_status: str):
        """Emit job status change signal."""
        data = {"job_id": job_id, "old_status": old_status, "new_status": new_status}
        logger.info("Job status changed: %s -> %s", old_status, new_status)
        self.job_status_changed.emit(data)


class ModernPluginScript(HierarchicalObject):
    """
    Modern replacement for CCP4PluginScript with integrated signal/hierarchy system.

    Provides:
    - Built-in signal system for task monitoring
    - Proper parent/child relationships
    - Async task execution
    - Error handling and recovery
    """

    def __init__(
        self,
        task_name: str,
        work_directory: str = None,
        parent: Optional[HierarchicalObject] = None,
        db_handler: Optional[ModernDbApi] = None,
    ):
        super().__init__(parent, f"PluginScript_{task_name}")

        self.task_name = task_name
        self.work_directory = Path(work_directory) if work_directory else Path.cwd()
        self.db_handler = db_handler

        # Create container for parameters (replaces Qt container)
        self.container = DataContainer(parent=self, name="ParameterContainer")

        # Plugin-specific signals
        self.task_started = self.create_signal("task_started", dict)
        self.task_progress = self.create_signal("task_progress", dict)
        self.task_finished = self.create_signal("task_finished", dict)
        self.task_failed = self.create_signal("task_failed", dict)
        self.close_app = self.create_signal("close_app", dict)

        # Connect to database signals if handler provided
        if self.db_handler:
            self.db_handler.data_changed.connect(self._on_db_data_changed)

    def _on_db_data_changed(self, data: Dict[str, Any]):
        """Handle database changes that affect this plugin."""
        logger.info("Plugin %s received DB change: %s", self.task_name, data)
        # Handle database updates relevant to this plugin

    async def process(self) -> Dict[str, Any]:
        """
        Main processing method - replaces Qt-based process().

        Override this in subclasses to implement specific task logic.
        """
        self.task_started.emit({"task_name": self.task_name})

        try:
            # Simulate some work
            await asyncio.sleep(0.1)

            result = await self._execute_task()

            self.task_finished.emit({"task_name": self.task_name, "result": result})

            # Signal completion (replaces Qt's closeApp signal)
            self.close_app.emit(
                {"jobId": self.name, "status": "completed", "result": result}
            )

            return result

        except Exception as e:
            error_data = {
                "task_name": self.task_name,
                "error": str(e),
                "error_type": type(e).__name__,
            }

            self.task_failed.emit(error_data)
            self.close_app.emit(
                {"jobId": self.name, "status": "failed", "error": error_data}
            )

            raise

    async def _execute_task(self) -> Dict[str, Any]:
        """
        Override this method in subclasses to implement actual task logic.
        """
        # Default implementation - just return success
        return {"status": "completed", "message": "Task executed successfully"}

    def update_progress(self, progress: float, message: str = ""):
        """Update task progress (0.0 to 1.0)."""
        self.task_progress.emit(
            {"task_name": self.task_name, "progress": progress, "message": message}
        )

    def set_parameter(self, path: str, value: Any):
        """Set a parameter in the container."""
        self.container.set_data(path, value)

    def get_parameter(self, path: str, default: Any = None) -> Any:
        """Get a parameter from the container."""
        return self.container.get_data(path, default)


class ModernTaskManager(HierarchicalObject):
    """
    Modern replacement for task management functionality.

    Manages plugin instances and task execution without Qt dependencies.
    """

    def __init__(self, parent: Optional[HierarchicalObject] = None):
        super().__init__(parent, "ModernTaskManager")

        self._plugin_classes: Dict[str, type] = {}
        self._active_plugins: Dict[str, ModernPluginScript] = {}

        # Signals
        self.plugin_registered = self.create_signal("plugin_registered", dict)
        self.plugin_started = self.create_signal("plugin_started", ModernPluginScript)
        self.plugin_completed = self.create_signal("plugin_completed", dict)

    def register_plugin_class(self, task_name: str, plugin_class: type):
        """Register a plugin class for a task name."""
        self._plugin_classes[task_name] = plugin_class
        self.plugin_registered.emit(
            {"task_name": task_name, "class_name": plugin_class.__name__}
        )
        logger.info(
            "Registered plugin class %s for task %s", plugin_class.__name__, task_name
        )

    def get_plugin_class(self, task_name: str) -> Optional[type]:
        """Get the plugin class for a task name."""
        return self._plugin_classes.get(task_name)

    def create_plugin(
        self,
        task_name: str,
        work_directory: str = None,
        db_handler: Optional[ModernDbApi] = None,
    ) -> Optional[ModernPluginScript]:
        """Create a plugin instance."""
        plugin_class = self.get_plugin_class(task_name)
        if not plugin_class:
            logger.error("No plugin class registered for task %s", task_name)
            return None

        try:
            plugin = plugin_class(
                task_name=task_name,
                work_directory=work_directory,
                parent=self,
                db_handler=db_handler,
            )

            self._active_plugins[plugin.name] = plugin
            self.plugin_started.emit(plugin)

            # Connect to completion signals
            plugin.close_app.connect(
                lambda data: self._on_plugin_completed(plugin.name, data)
            )

            return plugin

        except Exception as e:
            logger.error("Failed to create plugin for task %s: %s", task_name, e)
            return None

    def _on_plugin_completed(self, plugin_name: str, completion_data: Dict[str, Any]):
        """Handle plugin completion."""
        plugin = self._active_plugins.pop(plugin_name, None)
        if plugin:
            self.plugin_completed.emit(
                {"plugin_name": plugin_name, "completion_data": completion_data}
            )
            logger.info("Plugin %s completed", plugin_name)


class ModernApplication(Application):
    """
    Modern CCP4i2 application that integrates all the new systems.

    Replaces Qt application structure with:
    - Modern event loop
    - Signal-based communication
    - Hierarchical object management
    - Async task execution
    """

    def __init__(self):
        super().__init__("CCP4i2ModernApp")

        # Core components
        self.db_api = ModernDbApi(parent=self)
        self.task_manager = ModernTaskManager(parent=self)

        # Application signals
        self.job_submitted = self.create_signal("job_submitted", dict)
        self.job_completed = self.create_signal("job_completed", dict)

        # Connect inter-component signals
        self._setup_signal_connections()

    def _setup_signal_connections(self):
        """Setup signal connections between components."""
        # Connect task manager to database
        self.task_manager.plugin_completed.connect(self._on_job_completed)

        # Connect database events to application events
        self.db_api.job_status_changed.connect(self._on_job_status_changed)

    def _on_job_completed(self, completion_data: Dict[str, Any]):
        """Handle job completion."""
        self.job_completed.emit(completion_data)

        # Update database
        if "plugin_name" in completion_data:
            self.db_api.emit_job_status_change(
                job_id=completion_data["plugin_name"],
                old_status="running",
                new_status="completed",
            )

    def _on_job_status_changed(self, status_data: Dict[str, Any]):
        """Handle job status changes."""
        logger.info("Job status changed: %s", status_data)

    async def submit_job(
        self,
        task_name: str,
        work_directory: str = None,
        parameters: Dict[str, Any] = None,
    ) -> Optional[str]:
        """
        Submit a job for execution.

        Returns:
            Job ID if successful, None if failed
        """
        # Create plugin instance
        plugin = self.task_manager.create_plugin(
            task_name=task_name, work_directory=work_directory, db_handler=self.db_api
        )

        if not plugin:
            return None

        # Set parameters
        if parameters:
            for key, value in parameters.items():
                plugin.set_parameter(key, value)

        # Emit job submission signal
        job_data = {
            "job_id": plugin.name,
            "task_name": task_name,
            "work_directory": work_directory,
            "parameters": parameters or {},
        }
        self.job_submitted.emit(job_data)

        # Schedule the job for execution
        task_id = self.schedule_task(plugin.process, name=f"Job_{plugin.name}")

        logger.info("Submitted job %s (task_id: %s)", plugin.name, task_id)
        return plugin.name


# Example concrete plugin implementation
class ExampleRefmacPlugin(ModernPluginScript):
    """
    Example plugin that demonstrates the new architecture.
    """

    async def _execute_task(self) -> Dict[str, Any]:
        """Implement actual refinement task."""
        logger.info("Starting Refmac refinement...")

        # Simulate refinement steps
        steps = ["Reading input", "Setting up", "Refinement", "Writing output"]

        for i, step in enumerate(steps):
            self.update_progress((i + 1) / len(steps), step)
            await asyncio.sleep(0.5)  # Simulate work
            logger.info("Refmac: %s", step)

        # Simulate results
        result = {
            "status": "completed",
            "r_factor": 0.185,
            "r_free": 0.221,
            "output_files": ["refined.pdb", "refined.mtz"],
        }

        return result


# Usage example
async def main():
    """Demonstrate the modern CCP4i2 architecture."""
    logging.basicConfig(level=logging.INFO)

    # Create application
    app = ModernApplication()

    # Register plugin
    app.task_manager.register_plugin_class("refmac", ExampleRefmacPlugin)

    # Setup monitoring
    def on_job_completed(completion_data):
        print(f"Job completed: {completion_data}")

    app.job_completed.connect(on_job_completed)

    # Start application
    await app.start()

    # Submit a job
    job_id = await app.submit_job(
        task_name="refmac",
        work_directory="/tmp/test_job",
        parameters={"input_pdb": "input.pdb", "input_mtz": "input.mtz", "cycles": 10},
    )

    print(f"Submitted job: {job_id}")

    # Wait for completion
    await asyncio.sleep(3.0)

    # Shutdown
    await app.shutdown()
    print("Application shutdown complete")


if __name__ == "__main__":
    asyncio.run(main())
