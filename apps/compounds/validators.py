"""
File upload validators for compounds app.

Provides reusable validation for file uploads including:
- File extension restrictions
- File size limits
- Content type validation
"""

import os
from django.core.exceptions import ValidationError


# Maximum file sizes (in bytes)
MAX_QC_FILE_SIZE = 50 * 1024 * 1024  # 50 MB for QC files (NMR, LCMS PDFs)
MAX_DOCUMENT_SIZE = 25 * 1024 * 1024  # 25 MB for protocol documents
MAX_DATA_FILE_SIZE = 100 * 1024 * 1024  # 100 MB for assay data files
MAX_GENBANK_SIZE = 10 * 1024 * 1024  # 10 MB for GenBank files
MAX_SEQUENCE_FILE_SIZE = 5 * 1024 * 1024  # 5 MB for sequencing files

# Common document and image extensions (shared across QC files and protocol documents)
COMMON_DOCUMENT_EXTENSIONS = {
    # PDF
    '.pdf',
    # Microsoft Office
    '.doc', '.docx', '.xls', '.xlsx', '.ppt', '.pptx',
    # OpenDocument (LibreOffice, OpenOffice)
    '.odt', '.ods', '.odp', '.odg',
    # Text and data
    '.txt', '.csv', '.tsv', '.md', '.rtf',
    # Images
    '.png', '.jpg', '.jpeg', '.gif', '.tif', '.tiff', '.bmp', '.webp', '.svg',
    # Archives (for bundled data)
    '.zip',
}

# Allowed file extensions by category
QC_FILE_EXTENSIONS = COMMON_DOCUMENT_EXTENSIONS | {
    # NMR formats
    '.mnova', '.jdf', '.fid',
    # MS formats
    '.raw', '.d', '.wiff', '.mzml', '.mzxml',
    # HPLC/chromatography
    '.cdf', '.aia',
}

DOCUMENT_EXTENSIONS = COMMON_DOCUMENT_EXTENSIONS

DATA_FILE_EXTENSIONS = {
    '.xlsx', '.xls', '.csv', '.txt', '.tsv',
    '.xml', '.json',
}

GENBANK_EXTENSIONS = {
    '.gb', '.gbk', '.genbank', '.gbff',
    '.fasta', '.fa', '.fna', '.seq',
}

SEQUENCE_FILE_EXTENSIONS = {
    '.ab1', '.abi',  # ABI sequencing traces
    '.seq', '.fasta', '.fa',
    '.txt', '.pdf',
}


def validate_file_extension(file, allowed_extensions, file_type_name="file"):
    """
    Validate that uploaded file has an allowed extension.

    Args:
        file: Uploaded file object
        allowed_extensions: Set of allowed extensions (lowercase, with dot)
        file_type_name: Human-readable name for error messages

    Raises:
        ValidationError if extension is not allowed
    """
    if not file:
        return

    ext = os.path.splitext(file.name)[1].lower()
    if ext not in allowed_extensions:
        allowed_list = ', '.join(sorted(allowed_extensions))
        raise ValidationError(
            f"Invalid {file_type_name} extension '{ext}'. "
            f"Allowed extensions: {allowed_list}"
        )


def validate_file_size(file, max_size, file_type_name="file"):
    """
    Validate that uploaded file doesn't exceed maximum size.

    Args:
        file: Uploaded file object
        max_size: Maximum size in bytes
        file_type_name: Human-readable name for error messages

    Raises:
        ValidationError if file is too large
    """
    if not file:
        return

    if file.size > max_size:
        max_mb = max_size / (1024 * 1024)
        file_mb = file.size / (1024 * 1024)
        raise ValidationError(
            f"{file_type_name.capitalize()} is too large ({file_mb:.1f} MB). "
            f"Maximum size is {max_mb:.0f} MB."
        )


def validate_qc_file(file):
    """Validate a batch QC file upload."""
    validate_file_extension(file, QC_FILE_EXTENSIONS, "QC file")
    validate_file_size(file, MAX_QC_FILE_SIZE, "QC file")


def validate_protocol_document(file):
    """Validate a protocol document upload."""
    validate_file_extension(file, DOCUMENT_EXTENSIONS, "document")
    validate_file_size(file, MAX_DOCUMENT_SIZE, "document")


def validate_assay_data_file(file):
    """Validate an assay data file upload."""
    validate_file_extension(file, DATA_FILE_EXTENSIONS, "data file")
    validate_file_size(file, MAX_DATA_FILE_SIZE, "data file")


def validate_genbank_file(file):
    """Validate a GenBank/sequence file upload."""
    validate_file_extension(file, GENBANK_EXTENSIONS, "GenBank file")
    validate_file_size(file, MAX_GENBANK_SIZE, "GenBank file")


def validate_sequencing_file(file):
    """Validate a sequencing result file upload."""
    validate_file_extension(file, SEQUENCE_FILE_EXTENSIONS, "sequencing file")
    validate_file_size(file, MAX_SEQUENCE_FILE_SIZE, "sequencing file")
