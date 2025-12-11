from ccp4i2.core.CCP4TaskManager import CTaskManager


def get_task_tree():
    """
    Build the task tree structure for the frontend.

    Returns:
        Dict with:
        - tree: List of [module_name, module_title, [task_names...]] tuples
        - lookup: Dict mapping taskName -> {version: {metadata...}}
        - iconLookup: Dict mapping module name to icon path
    """
    task_manager = CTaskManager()
    result = {
        "tree": task_manager.task_tree(),
        "lookup": task_manager.task_lookup,
        "iconLookup": task_manager.task_icon_lookup,
    }
    return result
