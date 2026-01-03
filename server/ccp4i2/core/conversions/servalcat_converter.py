"""
Servalcat-based French-Wilson converter for CObsDataFile.

This module provides modern French-Wilson intensity-to-amplitude conversions
using servalcat fw (servalcat French-Wilson) as a replacement for ctruncate.

Servalcat fw provides improved statistics based on sounder algorithms for:
- IPAIR (I+/I-) → FPAIR (F+/F-) + FMEAN
- IPAIR → IMEAN + FMEAN
- IMEAN → FMEAN

Key differences from ctruncate:
- Uses subprocess calls instead of wrapper classes
- Executes in 'servalcat' subdirectory for log capture
- More accurate French-Wilson statistics
- Simpler invocation pattern
"""

import subprocess
from pathlib import Path
from typing import Optional, Any, Tuple
from ccp4i2.core.CCP4ErrorHandling import CException, SEVERITY_ERROR


class ServalcatConverter:
    """
    Converter using servalcat fw for French-Wilson intensity-to-amplitude conversions.

    This replaces ctruncate-based conversions with modern servalcat fw subcommand.
    """

    # Error codes for servalcat conversions
    ERROR_CODES = {
        10: {
            'description': 'servalcat executable not found in PATH',
            'severity': SEVERITY_ERROR},
        11: {
            'description': 'servalcat fw conversion failed',
            'severity': SEVERITY_ERROR},
        12: {
            'description': 'Expected output file not created by servalcat',
            'severity': SEVERITY_ERROR},
        13: {
            'description': 'Unsupported conversion path for servalcat',
            'severity': SEVERITY_ERROR},
    }

    @staticmethod
    def _get_servalcat_work_dir(work_directory: Optional[Any] = None) -> Path:
        """
        Get or create the 'servalcat' subdirectory for conversion logs.

        Similar to how ctruncate uses a 'ctruncate' subdirectory, we use 'servalcat'
        to isolate program.xml and log files.

        Args:
            work_directory: Parent work directory (CCP4_JOBS/job_X)

        Returns:
            Path to servalcat subdirectory
        """
        if work_directory:
            servalcat_dir = Path(work_directory) / 'servalcat'
        else:
            # Fallback to temp directory
            import tempfile
            servalcat_dir = Path(tempfile.mkdtemp(prefix='servalcat_'))

        servalcat_dir.mkdir(parents=True, exist_ok=True)
        return servalcat_dir

    @staticmethod
    def _run_servalcat_fw(
        input_mtz: str,
        labin: str,
        output_prefix: str,
        work_dir: Path
    ) -> Tuple[int, str, str, Path]:
        """
        Execute servalcat fw command.

        Servalcat fw produces a monolithic MTZ with multiple outputs:
        - F, SIGF (FMEAN)
        - I, SIGI (IMEAN)
        - F(+), SIGF(+), F(-), SIGF(-) (FPAIR)

        Args:
            input_mtz: Path to input MTZ file
            labin: Column specification (e.g., "I,SIGI" or "Iplus,SIGIplus,Iminus,SIGIminus")
            output_prefix: Output file prefix
            work_dir: Working directory for servalcat execution

        Returns:
            Tuple of (return_code, stdout, stderr, monolithic_mtz_path)
        """
        cmd = [
            'servalcat', 'fw',
            '--hklin', str(input_mtz),
            '--labin', labin,
            '-o', output_prefix
        ]

        # Run in servalcat subdirectory to capture logs
        result = subprocess.run(
            cmd,
            cwd=str(work_dir),
            capture_output=True,
            text=True
        )

        # Servalcat creates a monolithic MTZ file
        monolithic_mtz = work_dir / f"{Path(output_prefix).name}.mtz"

        return result.returncode, result.stdout, result.stderr, monolithic_mtz

    @staticmethod
    def _split_monolithic_mtz(
        monolithic_mtz: Path,
        output_type: str,
        output_path: str,
        obs_file
    ) -> str:
        """
        Split monolithic servalcat MTZ to extract specific content type.

        Servalcat fw produces a monolithic MTZ with columns:
        - F, SIGF (FMEAN)
        - I, SIGI (IMEAN)
        - F(+), SIGF(+), F(-), SIGF(-) (FPAIR)

        Args:
            monolithic_mtz: Path to servalcat's monolithic output
            output_type: Target content type ('FPAIR', 'IMEAN', 'FMEAN')
            output_path: Desired output file path
            obs_file: CObsDataFile instance (unused, kept for API compatibility)

        Returns:
            Path to extracted mini-MTZ
        """
        from ccp4i2.core.conversions.mtz_utils import extract_columns

        # Column mappings for servalcat fw output
        # Maps output_type → list of (source_col, target_col, mtz_type) tuples
        column_mappings = {
            'FMEAN': [
                ('F', 'F', 'F'),
                ('SIGF', 'SIGF', 'Q')
            ],
            'IMEAN': [
                ('I', 'I', 'J'),
                ('SIGI', 'SIGI', 'Q')
            ],
            'FPAIR': [
                ('F(+)', 'Fplus', 'G'),
                ('SIGF(+)', 'SIGFplus', 'L'),
                ('F(-)', 'Fminus', 'G'),
                ('SIGF(-)', 'SIGFminus', 'L')
            ],
        }

        if output_type not in column_mappings:
            raise CException(
                ServalcatConverter, 13,
                details=f"Unknown output type: {output_type}"
            )

        column_spec = column_mappings[output_type]

        try:
            extract_columns(
                str(monolithic_mtz),
                str(output_path),
                column_spec
            )
        except RuntimeError as e:
            raise CException(
                ServalcatConverter, 12,
                details=f"Failed to extract {output_type} columns: {str(e)}"
            ) from e

        return str(output_path)

    @staticmethod
    def ipair_to_fpair(obs_file, work_directory: Optional[Any] = None) -> str:
        """
        Convert IPAIR (I+/I-) to FPAIR (F+/F-) using servalcat fw.

        Uses French-Wilson conversion. Servalcat fw produces a monolithic MTZ
        with FMEAN, IMEAN, and FPAIR - we extract just FPAIR.

        Args:
            obs_file: CObsDataFile instance with IPAIR data
            work_directory: Optional directory for servalcat working files

        Returns:
            Full path to converted FPAIR file

        Raises:
            CException: If servalcat fails or output not created
        """
        # Get column names from CONTENT_SIGNATURE_LIST
        column_names = obs_file.CONTENT_SIGNATURE_LIST[obs_file.CONTENT_FLAG_IPAIR - 1]

        # Servalcat expects: Iplus,SIGIplus,Iminus,SIGIminus
        labin = f"{column_names[0]},{column_names[1]},{column_names[2]},{column_names[3]}"

        # Get output path
        output_path = obs_file._get_conversion_output_path(
            'FPAIR', work_directory=work_directory)
        output_path_obj = Path(output_path)

        # Create servalcat work directory
        servalcat_dir = ServalcatConverter._get_servalcat_work_dir(
            work_directory)

        # Run servalcat fw with output prefix (use base name without extension)
        output_prefix = output_path_obj.stem

        returncode, stdout, stderr, monolithic_mtz = ServalcatConverter._run_servalcat_fw(
            input_mtz=obs_file.getFullPath(),
            labin=labin,
            output_prefix=str(servalcat_dir / output_prefix),
            work_dir=servalcat_dir
        )

        if returncode != 0:
            error_msg = f"servalcat fw failed with return code {returncode}"
            if stderr:
                error_msg += f"\nStderr: {stderr}"
            if stdout:
                error_msg += f"\nStdout: {stdout}"
            raise CException(ServalcatConverter, 11, details=error_msg)

        # Verify monolithic MTZ exists
        if not monolithic_mtz.exists():
            raise CException(
                ServalcatConverter, 12,
                details=f"Servalcat did not create output: {monolithic_mtz}"
            )

        # Extract FPAIR from monolithic MTZ
        ServalcatConverter._split_monolithic_mtz(
            monolithic_mtz, 'FPAIR', output_path, obs_file)

        print(f"✅ Servalcat IPAIR → FPAIR conversion: {output_path}")
        return output_path

    @staticmethod
    def ipair_to_imean(obs_file, work_directory: Optional[Any] = None) -> str:
        """
        Convert IPAIR (I+/I-) to IMEAN (I, SIGI) using servalcat fw.

        Uses French-Wilson conversion. Servalcat fw produces a monolithic MTZ
        with FMEAN, IMEAN, and FPAIR - we extract just IMEAN.

        Args:
            obs_file: CObsDataFile instance with IPAIR data
            work_directory: Optional directory for servalcat working files

        Returns:
            Full path to converted IMEAN file

        Raises:
            CException: If servalcat fails or output not created
        """
        # Get column names
        column_names = obs_file.CONTENT_SIGNATURE_LIST[obs_file.CONTENT_FLAG_IPAIR - 1]
        labin = f"{column_names[0]},{column_names[1]},{column_names[2]},{column_names[3]}"

        # Get output path
        output_path = obs_file._get_conversion_output_path(
            'IMEAN', work_directory=work_directory)
        output_path_obj = Path(output_path)

        # Create servalcat work directory
        servalcat_dir = ServalcatConverter._get_servalcat_work_dir(
            work_directory)

        # Run servalcat fw
        output_prefix = output_path_obj.stem
        returncode, stdout, stderr, monolithic_mtz = ServalcatConverter._run_servalcat_fw(
            input_mtz=obs_file.getFullPath(),
            labin=labin,
            output_prefix=str(servalcat_dir / output_prefix),
            work_dir=servalcat_dir
        )

        if returncode != 0:
            error_msg = f"servalcat fw failed with return code {returncode}"
            if stderr:
                error_msg += f"\nStderr: {stderr}"
            raise CException(ServalcatConverter, 11, details=error_msg)

        # Verify monolithic MTZ exists
        if not monolithic_mtz.exists():
            raise CException(
                ServalcatConverter, 12,
                details=f"Servalcat did not create output: {monolithic_mtz}"
            )

        # Extract IMEAN from monolithic MTZ
        ServalcatConverter._split_monolithic_mtz(
            monolithic_mtz, 'IMEAN', output_path, obs_file)

        print(f"✅ Servalcat IPAIR → IMEAN conversion: {output_path}")
        return output_path

    @staticmethod
    def to_fmean(obs_file, work_directory: Optional[Any] = None) -> str:
        """
        Convert to FMEAN (F, SIGF) using servalcat fw.

        Handles:
        - IPAIR → FMEAN: French-Wilson via servalcat
        - IMEAN → FMEAN: French-Wilson via servalcat

        Servalcat fw produces a monolithic MTZ with FMEAN, IMEAN, and FPAIR.

        Args:
            obs_file: CObsDataFile instance
            work_directory: Optional directory for servalcat working files

        Returns:
            Full path to converted FMEAN file

        Raises:
            CException: If conversion fails
        """
        current_flag = int(obs_file.contentFlag)

        # Get column names
        column_names = obs_file.CONTENT_SIGNATURE_LIST[current_flag - 1]

        # Build labin string
        if current_flag == obs_file.CONTENT_FLAG_IPAIR:
            # IPAIR: ['Iplus', 'SIGIplus', 'Iminus', 'SIGIminus']
            labin = f"{column_names[0]},{column_names[1]},{column_names[2]},{column_names[3]}"
        elif current_flag == obs_file.CONTENT_FLAG_IMEAN:
            # IMEAN: ['I', 'SIGI']
            labin = f"{column_names[0]},{column_names[1]}"
        else:
            raise CException(
                ServalcatConverter,
                13,
                details=f"contentFlag {current_flag} → FMEAN not supported via servalcat")

        # Get output path
        output_path = obs_file._get_conversion_output_path(
            'FMEAN', work_directory=work_directory)
        output_path_obj = Path(output_path)

        # Create servalcat work directory
        servalcat_dir = ServalcatConverter._get_servalcat_work_dir(
            work_directory)

        # Run servalcat fw
        output_prefix = output_path_obj.stem
        returncode, stdout, stderr, monolithic_mtz = ServalcatConverter._run_servalcat_fw(
            input_mtz=obs_file.getFullPath(),
            labin=labin,
            output_prefix=str(servalcat_dir / output_prefix),
            work_dir=servalcat_dir
        )

        if returncode != 0:
            error_msg = f"servalcat fw failed with return code {returncode}"
            if stderr:
                error_msg += f"\nStderr: {stderr}"
            raise CException(ServalcatConverter, 11, details=error_msg)

        # Verify monolithic MTZ exists
        if not monolithic_mtz.exists():
            raise CException(
                ServalcatConverter, 12,
                details=f"Servalcat did not create output: {monolithic_mtz}"
            )

        # Extract FMEAN from monolithic MTZ
        ServalcatConverter._split_monolithic_mtz(
            monolithic_mtz, 'FMEAN', output_path, obs_file)

        print(f"✅ Servalcat → FMEAN conversion: {output_path}")
        return output_path
