"""
Ctruncate-based French-Wilson converter for CObsDataFile.

This module provides legacy French-Wilson intensity-to-amplitude conversions
using the ctruncate wrapper. It provides the SAME API as ServalcatConverter,
allowing easy A/B comparison testing between ctruncate and servalcat fw.

Ctruncate handles:
- IPAIR (I+/I-) → FPAIR (F+/F-) + FMEAN
- IPAIR → IMEAN + FMEAN
- IMEAN → FMEAN

Usage:
    # Drop-in replacement for ServalcatConverter
    from ccp4i2.core.conversions.ctruncate_converter import CtruncateConverter

    # Same API as ServalcatConverter
    output_path = CtruncateConverter.to_fmean(obs_file, work_directory)
    output_path = CtruncateConverter.ipair_to_fpair(obs_file, work_directory)
    output_path = CtruncateConverter.ipair_to_imean(obs_file, work_directory)
"""

from pathlib import Path
from typing import Optional, Any
from ccp4i2.core.CCP4ErrorHandling import CException, SEVERITY_ERROR


class CtruncateConverter:
    """
    Converter using ctruncate wrapper for French-Wilson conversions.

    Provides the same API as ServalcatConverter to allow easy comparison
    between ctruncate (legacy) and servalcat (modern) implementations.
    """

    # Error codes for ctruncate conversions
    ERROR_CODES = {
        20: {
            'description': 'ctruncate plugin execution failed',
            'severity': SEVERITY_ERROR},
        21: {
            'description': 'Expected output file not created by ctruncate',
            'severity': SEVERITY_ERROR},
        22: {
            'description': 'ctruncate plugin not available in registry',
            'severity': SEVERITY_ERROR},
        23: {
            'description': 'Unsupported conversion path for ctruncate',
            'severity': SEVERITY_ERROR},
    }

    @staticmethod
    def _get_ctruncate_work_dir(work_directory: Optional[Any] = None) -> Path:
        """
        Get or create the 'ctruncate' subdirectory for conversion logs.

        Args:
            work_directory: Parent work directory (CCP4_JOBS/job_X)

        Returns:
            Path to ctruncate subdirectory
        """
        if work_directory:
            ctruncate_dir = Path(work_directory) / 'ctruncate'
        else:
            # Fallback to temp directory
            import tempfile
            ctruncate_dir = Path(tempfile.mkdtemp(prefix='ctruncate_'))

        ctruncate_dir.mkdir(parents=True, exist_ok=True)
        return ctruncate_dir

    @staticmethod
    def _run_ctruncate(
        input_mtz: str,
        output_mtz: str,
        content_flag_in: int,
        content_flag_out: int,
        work_dir: Path
    ) -> None:
        """
        Execute ctruncate plugin for French-Wilson conversion.

        Args:
            input_mtz: Path to input MTZ file
            output_mtz: Path to output MTZ file
            content_flag_in: Input contentFlag (1=IPAIR, 2=FPAIR, 3=IMEAN, 4=FMEAN)
            content_flag_out: Output contentFlag
            work_dir: Working directory for ctruncate execution

        Raises:
            CException: If ctruncate plugin fails or output not created
        """
        from ccp4i2.core.task_manager.plugin_registry import get_plugin_class
        from ccp4i2.core.CCP4XtalData import CObsDataFile

        # Get ctruncate plugin class
        try:
            CtruncatePlugin = get_plugin_class('ctruncate')
        except Exception as e:
            raise CException(
                CtruncateConverter, 22,
                details=f"Failed to import ctruncate plugin: {str(e)}"
            )

        # Create ctruncate plugin instance
        ctruncate_plugin = CtruncatePlugin(parent=None, name='ctruncate_converter')
        ctruncate_plugin.workDirectory = work_dir

        # Set up input file
        input_obs = CObsDataFile(parent=ctruncate_plugin.container.inputData, name='OBSIN')
        input_obs.setFullPath(str(input_mtz))
        input_obs.contentFlag.value = content_flag_in
        ctruncate_plugin.container.inputData.OBSIN = input_obs

        # Set up output parameters
        ctruncate_plugin.container.controlParameters.OUTPUTMINIMTZ = True
        ctruncate_plugin.container.controlParameters.OUTPUTMINIMTZCONTENTFLAG = content_flag_out

        # Output intensities if converting to IMEAN
        if content_flag_out == 3:  # IMEAN
            ctruncate_plugin.container.controlParameters.OUTPUT_INTENSITIES = True

        # Execute the conversion
        try:
            # Build command
            ctruncate_plugin.makeCommandAndScript()

            # Run the process
            status = ctruncate_plugin.process()
            if status != ctruncate_plugin.SUCCEEDED:
                # Check error report
                if ctruncate_plugin.errorReport.maxSeverity() >= SEVERITY_ERROR:
                    error_details = ctruncate_plugin.errorReport.report()
                    raise CException(
                        CtruncateConverter, 20,
                        details=f"ctruncate execution failed:\n{error_details}"
                    )
                raise CException(
                    CtruncateConverter, 20,
                    details="ctruncate.process() failed"
                )

            # Process output files
            ctruncate_plugin.processOutputFiles()

            # Get the output mini-MTZ file
            if not ctruncate_plugin.container.outputData.HKLOUT.isSet():
                raise CException(
                    CtruncateConverter, 21,
                    details="ctruncate did not create HKLOUT"
                )

            # Copy to desired output location
            import shutil
            ctruncate_output = Path(ctruncate_plugin.container.outputData.HKLOUT.getFullPath())
            if not ctruncate_output.exists():
                raise CException(
                    CtruncateConverter, 21,
                    details=f"ctruncate output not found: {ctruncate_output}"
                )

            shutil.copy2(ctruncate_output, output_mtz)

        except CException:
            raise
        except Exception as e:
            raise CException(
                CtruncateConverter, 20,
                details=f"Unexpected error running ctruncate: {str(e)}"
            ) from e

    @staticmethod
    def ipair_to_fpair(obs_file, work_directory: Optional[Any] = None) -> str:
        """
        Convert IPAIR (I+/I-) to FPAIR (F+/F-) using ctruncate.

        Uses French-Wilson conversion. This method provides the same API
        as ServalcatConverter.ipair_to_fpair() for easy comparison testing.

        Args:
            obs_file: CObsDataFile instance with IPAIR data
            work_directory: Optional directory for ctruncate working files

        Returns:
            Full path to converted FPAIR file

        Raises:
            CException: If ctruncate fails or output not created
        """
        # Get output path
        output_path = obs_file._get_conversion_output_path(
            'FPAIR', work_directory=work_directory)

        # Create ctruncate work directory
        ctruncate_dir = CtruncateConverter._get_ctruncate_work_dir(work_directory)

        # Run ctruncate (IPAIR=1 → FPAIR=2)
        CtruncateConverter._run_ctruncate(
            input_mtz=obs_file.getFullPath(),
            output_mtz=output_path,
            content_flag_in=obs_file.CONTENT_FLAG_IPAIR,
            content_flag_out=obs_file.CONTENT_FLAG_FPAIR,
            work_dir=ctruncate_dir
        )

        print(f"✅ Ctruncate IPAIR → FPAIR conversion: {output_path}")
        return output_path

    @staticmethod
    def ipair_to_imean(obs_file, work_directory: Optional[Any] = None) -> str:
        """
        Convert IPAIR (I+/I-) to IMEAN (I, SIGI) using ctruncate.

        Uses French-Wilson conversion. This method provides the same API
        as ServalcatConverter.ipair_to_imean() for easy comparison testing.

        Args:
            obs_file: CObsDataFile instance with IPAIR data
            work_directory: Optional directory for ctruncate working files

        Returns:
            Full path to converted IMEAN file

        Raises:
            CException: If ctruncate fails or output not created
        """
        # Get output path
        output_path = obs_file._get_conversion_output_path(
            'IMEAN', work_directory=work_directory)

        # Create ctruncate work directory
        ctruncate_dir = CtruncateConverter._get_ctruncate_work_dir(work_directory)

        # Run ctruncate (IPAIR=1 → IMEAN=3)
        CtruncateConverter._run_ctruncate(
            input_mtz=obs_file.getFullPath(),
            output_mtz=output_path,
            content_flag_in=obs_file.CONTENT_FLAG_IPAIR,
            content_flag_out=obs_file.CONTENT_FLAG_IMEAN,
            work_dir=ctruncate_dir
        )

        print(f"✅ Ctruncate IPAIR → IMEAN conversion: {output_path}")
        return output_path

    @staticmethod
    def to_fmean(obs_file, work_directory: Optional[Any] = None) -> str:
        """
        Convert to FMEAN (F, SIGF) using ctruncate.

        Handles:
        - IPAIR → FMEAN: French-Wilson via ctruncate
        - IMEAN → FMEAN: French-Wilson via ctruncate

        This method provides the same API as ServalcatConverter.to_fmean()
        for easy comparison testing.

        Args:
            obs_file: CObsDataFile instance
            work_directory: Optional directory for ctruncate working files

        Returns:
            Full path to converted FMEAN file

        Raises:
            CException: If conversion fails
        """
        current_flag = int(obs_file.contentFlag)

        # Validate input format
        if current_flag not in (obs_file.CONTENT_FLAG_IPAIR, obs_file.CONTENT_FLAG_IMEAN):
            raise CException(
                CtruncateConverter,
                23,
                details=f"contentFlag {current_flag} → FMEAN not supported via ctruncate. "
                        f"Only IPAIR (1) and IMEAN (3) supported."
            )

        # Get output path
        output_path = obs_file._get_conversion_output_path(
            'FMEAN', work_directory=work_directory)

        # Create ctruncate work directory
        ctruncate_dir = CtruncateConverter._get_ctruncate_work_dir(work_directory)

        # Run ctruncate (IPAIR=1 or IMEAN=3 → FMEAN=4)
        CtruncateConverter._run_ctruncate(
            input_mtz=obs_file.getFullPath(),
            output_mtz=output_path,
            content_flag_in=current_flag,
            content_flag_out=obs_file.CONTENT_FLAG_FMEAN,
            work_dir=ctruncate_dir
        )

        print(f"✅ Ctruncate → FMEAN conversion: {output_path}")
        return output_path
