"""
Job report generation utilities.

Provides access to various job reports: parameters, execution results, diagnostics.
"""

import logging
from pathlib import Path
from xml.etree import ElementTree as ET
from ccp4i2.db import models
from ccp4i2.lib.response import Result
from ..reporting.i2_report import generate_job_report

logger = logging.getLogger(f"ccp4i2:{__name__}")


def get_job_params_xml(job: models.Job) -> Result[str]:
    """
    Get job parameters XML.

    Returns the job's parameter configuration in XML format, either from
    params.xml (for completed jobs) or input_params.xml (for pending jobs).

    Args:
        job: Job model instance

    Returns:
        Result containing XML string

    Example:
        >>> result = get_job_params_xml(job)
        >>> if result.success:
        ...     with open('output.xml', 'w') as f:
        ...         f.write(result.data)
    """
    try:
        # Determine which file to check first based on status
        if job.status in [models.Job.Status.UNKNOWN, models.Job.Status.PENDING]:
            primary_path = job.directory / "input_params.xml"
            fallback_path = job.directory / "params.xml"
        else:
            primary_path = job.directory / "params.xml"
            fallback_path = job.directory / "input_params.xml"

        # Try primary path
        try:
            with open(primary_path, "r", encoding="UTF-8") as f:
                return Result.ok(f.read())
        except FileNotFoundError:
            pass

        # Try fallback path
        try:
            with open(fallback_path, "r", encoding="UTF-8") as f:
                return Result.ok(f.read())
        except FileNotFoundError:
            return Result.fail(
                f"No params file found for job {job.uuid}",
                details={
                    "tried_paths": [str(primary_path), str(fallback_path)],
                    "job_status": models.Job.Status(job.status).label
                }
            )

    except Exception as err:
        logger.exception("Failed to get params XML for job %s", job.uuid, exc_info=err)
        return Result.fail(
            f"Failed to read params XML: {str(err)}",
            details={"job_id": str(job.uuid), "error_type": type(err).__name__}
        )


def get_job_report_xml(job: models.Job, regenerate: bool = False) -> Result[bytes]:
    """
    Generate and retrieve XML report for a job.

    Creates a comprehensive XML report containing job results, statistics,
    and analysis outcomes. Reports are cached unless regenerate=True.

    Args:
        job: Job model instance
        regenerate: If True, regenerate report even if cached

    Returns:
        Result containing XML bytes

    Example:
        >>> result = get_job_report_xml(job)
        >>> if result.success:
        ...     with open('report.xml', 'wb') as f:
        ...         f.write(result.data)
    """
    try:
        # For pending jobs, return simple status
        if job.status in [models.Job.Status.PENDING, models.Job.Status.UNKNOWN]:
            report_xml = ET.Element("report")
            ET.SubElement(report_xml, "status").text = "PENDING"
            ET.indent(report_xml, space="\t", level=0)
            return Result.ok(ET.tostring(report_xml))

        # For completed/failed jobs
        if job.status in [
            models.Job.Status.FAILED,
            models.Job.Status.UNSATISFACTORY,
            models.Job.Status.INTERRUPTED,
            models.Job.Status.FINISHED,
        ]:
            report_xml_path = job.directory / "report_xml.xml"

            # Use cached report if exists and not regenerating
            if report_xml_path.exists() and not regenerate:
                with open(report_xml_path, "rb") as f:
                    return Result.ok(f.read())

            # Generate new report
            report_xml = generate_job_report(job)
            ET.indent(report_xml, space="\t", level=0)
            xml_bytes = ET.tostring(report_xml)

            # Only cache successful reports (not error/placeholder reports)
            # Check if this is a failed report by looking for error indicators
            xml_str = xml_bytes.decode('utf-8', errors='ignore')
            is_error_report = 'No report because' in xml_str or 'Report generation failed' in xml_str

            if not is_error_report:
                with open(report_xml_path, "wb") as f:
                    f.write(xml_bytes)
            else:
                # Delete any stale cached report so we regenerate next time
                if report_xml_path.exists():
                    report_xml_path.unlink()

            return Result.ok(xml_bytes)

        # For running jobs
        report_xml = generate_job_report(job)
        ET.indent(report_xml, space="\t", level=0)
        return Result.ok(ET.tostring(report_xml))

    except Exception as err:
        logger.exception("Failed to generate report for job %s", job.uuid, exc_info=err)
        return Result.fail(
            f"Failed to generate report: {str(err)}",
            details={"job_id": str(job.uuid), "error_type": type(err).__name__}
        )


def get_job_diagnostic_xml(job: models.Job) -> Result[str]:
    """
    Retrieve diagnostic information as XML.

    Returns detailed diagnostic data generated during job execution,
    including error messages, warnings, and debugging information.

    Args:
        job: Job model instance

    Returns:
        Result containing diagnostic XML string

    Example:
        >>> result = get_job_diagnostic_xml(job)
        >>> if result.success:
        ...     print(result.data)
    """
    try:
        diagnostic_path = job.directory / "diagnostic.xml"

        if not diagnostic_path.exists():
            return Result.fail(
                f"Diagnostic file not found for job {job.uuid}",
                details={
                    "path": str(diagnostic_path),
                    "job_status": models.Job.Status(job.status).label
                }
            )

        with open(diagnostic_path, "r", encoding="UTF-8") as f:
            return Result.ok(f.read())

    except Exception as err:
        logger.exception("Failed to read diagnostic XML for job %s", job.uuid, exc_info=err)
        return Result.fail(
            f"Failed to read diagnostic XML: {str(err)}",
            details={"job_id": str(job.uuid), "error_type": type(err).__name__}
        )
