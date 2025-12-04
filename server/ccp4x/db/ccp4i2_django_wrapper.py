import logging
import traceback
from core import CCP4ProjectsManager
from .ccp4i2_django_projects_manager import CCP4i2DjangoProjectsManager

logger = logging.getLogger(f"ccp4x:{__name__}")


# Decoorator to install and use FakeProjectManager
def using_django_pm(func):
    """
    Decorator to temporarily replace the CCP4ProjectsManager instance with a
    CCP4i2DjangoProjectsManager instance while the decorated function is executed.
    This decorator is useful for ensuring that the decorated function uses the
    Django-based project manager instead of the default one.
    Args:
        func (callable): The function to be decorated.
    Returns:
        callable: The wrapped function with the temporary project manager replacement.
    The decorator performs the following steps:
    1. Logs a debug message before the function is called.
    2. Saves the current CCP4ProjectsManager instance.
    3. Replaces the CCP4ProjectsManager instance with a CCP4i2DjangoProjectsManager instance.
    4. Executes the decorated function.
    5. If an exception occurs, logs the exception and prints the traceback.
    6. Restores the original CCP4ProjectsManager instance.
    7. Logs a debug message after the function is called.
    """

    def wrapper(*args, **kwargs):
        logger.debug("Something is happening before the function is called.")
        oldPM = CCP4ProjectsManager.CProjectsManager.insts
        result = None  # Initialize result before try block
        try:
            CCP4ProjectsManager.CProjectsManager.insts = CCP4i2DjangoProjectsManager()
            result = func(*args, **kwargs)
        except Exception as err:
            logging.exception(
                "Encountered issue while in FakePM decorator", exc_info=err
            )
            traceback.print_exc()
        finally:
            if oldPM is not None:
                CCP4ProjectsManager.CProjectsManager.insts = oldPM
            logger.debug("Something is happening after the function is called.")
        return result if result is not None else ""

    return wrapper
