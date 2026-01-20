"""
Shared utilities for the compounds app.
"""

import logging
import os

logger = logging.getLogger(__name__)


def delete_file_field(file_field, save=False):
    """
    Delete a file from a FileField, handling both local and cloud storage.

    This function safely deletes files regardless of the storage backend:
    - Local filesystem: Uses os.remove() after checking file exists
    - Cloud storage (Azure Blob, S3, etc.): Uses storage backend's delete()

    Args:
        file_field: A Django FileField instance (e.g., instance.genbank_file)
        save: Whether to save the model after clearing the field (default False)

    Returns:
        True if file was deleted, False if no file to delete
    """
    if not file_field:
        return False

    try:
        # Try local filesystem deletion first
        file_path = file_field.path
        if os.path.isfile(file_path):
            os.remove(file_path)
            logger.debug(f"Deleted local file: {file_path}")
            return True
    except NotImplementedError:
        # Cloud storage backends don't support .path
        # Use the storage backend's delete method instead
        try:
            file_field.delete(save=save)
            logger.debug(f"Deleted cloud file: {file_field.name}")
            return True
        except Exception as e:
            logger.warning(f"Failed to delete cloud file {file_field.name}: {e}")
            return False
    except Exception as e:
        logger.warning(f"Failed to delete file: {e}")
        return False

    return False
