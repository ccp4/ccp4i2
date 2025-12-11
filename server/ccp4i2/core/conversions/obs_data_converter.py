"""
Converter for CObsDataFile MTZ observation data format transformations.

This module handles conversions between different MTZ observation data formats:
- IPAIR (1): Anomalous intensities (I+, SIGI+, I-, SIGI-)
- FPAIR (2): Anomalous structure factors (F+, SIGF+, F-, SIGF-)
- IMEAN (3): Mean intensities (I, SIGI)
- FMEAN (4): Mean structure factors (F, SIGF)

Conversion Matrix:
                    TO
          IPAIR FPAIR IMEAN FMEAN
FROM IPAIR   ✓    ✓     ✓     ✓
     FPAIR   ✗    ✓     ✗     ✓
     IMEAN   ✗    ✗     ✓     ✓
     FMEAN   ✗    ✗     ✗     ✓

Implementation Details:
- IPAIR → FPAIR/IMEAN/FMEAN: Uses servalcat fw (French-Wilson)
- IMEAN → FMEAN: Uses servalcat fw (French-Wilson)
- FPAIR → FMEAN: Uses gemmi + numpy (inverse-variance weighted mean)
"""

from typing import Optional, Any
from ccp4i2.core.CCP4ErrorHandling import CException, CErrorReport, SEVERITY_ERROR


class ObsDataConverter:
    """
    Static converter class for CObsDataFile contentFlag transformations.

    All methods are static and take the CObsDataFile instance as first parameter.
    This design allows the converter to be used independently while keeping
    the intuitive obs_file.as_FMEAN() API through thin wrapper methods.
    """

    # Error codes for observation data conversions
    ERROR_CODES = {
        1: {'description': 'Input file does not exist', 'severity': SEVERITY_ERROR},
        2: {'description': 'Cannot determine contentFlag from input file', 'severity': SEVERITY_ERROR},
        3: {'description': 'Unsupported conversion path', 'severity': SEVERITY_ERROR},
        4: {'description': 'ctruncate plugin not available', 'severity': SEVERITY_ERROR},
        5: {'description': 'ctruncate conversion failed', 'severity': SEVERITY_ERROR},
        6: {'description': 'Output file not created after conversion', 'severity': SEVERITY_ERROR},
        7: {'description': 'Required anomalous data columns not found', 'severity': SEVERITY_ERROR},
        8: {'description': 'Invalid sigma values in anomalous data', 'severity': SEVERITY_ERROR},
    }

    @staticmethod
    def _validate_input_file(obs_file):
        """
        Validate input file before conversion.

        Args:
            obs_file: CObsDataFile instance

        Raises:
            CException: If file doesn't exist or contentFlag cannot be determined
        """
        from pathlib import Path

        input_path = obs_file.getFullPath()

        # Check file exists
        if not Path(input_path).exists():
            raise CException(ObsDataConverter, 1,
                             details=f"File: {input_path}")

        # Check contentFlag is set
        content_flag = int(
            obs_file.contentFlag) if obs_file.contentFlag.isSet() else 0
        if content_flag == 0:
            raise CException(ObsDataConverter, 2,
                             details=f"File: {input_path}")

    @staticmethod
    def to_fpair(obs_file, work_directory: Optional[Any] = None) -> str:
        """
        Convert IPAIR (I+/I-) to FPAIR (F+/F-) using servalcat fw.

        Uses French-Wilson conversion to transform anomalous intensities
        to anomalous structure factor amplitudes.

        Args:
            obs_file: CObsDataFile instance to convert
            work_directory: Optional directory for servalcat working files

        Returns:
            Full path to converted FPAIR file

        Raises:
            CException: If validation fails, unsupported conversion, or servalcat fails
        """
        from pathlib import Path
        import shutil
        from ccp4i2.core.conversions.servalcat_converter import ServalcatConverter

        # Auto-detect content flag and validate
        obs_file.setContentFlag()
        ObsDataConverter._validate_input_file(obs_file)

        current_flag = int(obs_file.contentFlag)

        # If already FPAIR, just copy
        if current_flag == obs_file.CONTENT_FLAG_FPAIR:
            output_path = obs_file._get_conversion_output_path(
                'FPAIR', work_directory=work_directory)
            input_path = Path(obs_file.getFullPath())
            shutil.copy2(input_path, output_path)
            return str(output_path)

        # Only IPAIR can convert to FPAIR
        if current_flag != obs_file.CONTENT_FLAG_IPAIR:
            raise CException(
                ObsDataConverter, 3,
                details=f"Source: contentFlag {current_flag}, Target: FPAIR (2). Only IPAIR (1) → FPAIR supported."
            )

        # Use servalcat for conversion
        return ServalcatConverter.ipair_to_fpair(obs_file, work_directory=work_directory)

    @staticmethod
    def to_imean(obs_file, work_directory: Optional[Any] = None) -> str:
        """
        Convert IPAIR (I+/I-) to IMEAN (I,SIGI) using servalcat fw.

        Averages anomalous intensities I+ and I- to produce mean intensities.

        Args:
            obs_file: CObsDataFile instance to convert
            work_directory: Optional directory for servalcat working files

        Returns:
            Full path to converted IMEAN file

        Raises:
            CException: If validation fails, unsupported conversion, or servalcat fails
        """
        from pathlib import Path
        import shutil
        from ccp4i2.core.conversions.servalcat_converter import ServalcatConverter

        # Auto-detect content flag from file
        obs_file.setContentFlag()
        current_flag = int(obs_file.contentFlag)

        # Validation handles this now, but kept for clarity
        if current_flag == 0:
            raise CException(ObsDataConverter, 2,
                             details=f"File: {obs_file.getFullPath()}")

        # If already IMEAN, just copy the file
        if current_flag == obs_file.CONTENT_FLAG_IMEAN:
            output_path = obs_file._get_conversion_output_path(
                'IMEAN', work_directory=work_directory)
            input_path = Path(obs_file.getFullPath())
            shutil.copy2(input_path, output_path)
            return str(output_path)

        # Only IPAIR can convert to IMEAN
        if current_flag != obs_file.CONTENT_FLAG_IPAIR:
            raise CException(
                ObsDataConverter, 3,
                details=f"Source: contentFlag {current_flag}, Target: IMEAN (3). Only IPAIR (1) → IMEAN supported."
            )

        # Use servalcat for conversion
        return ServalcatConverter.ipair_to_imean(obs_file, work_directory=work_directory)

    @staticmethod
    def to_fmean(obs_file, work_directory: Optional[Any] = None) -> str:
        """
        Convert to FMEAN format (Mean Structure Factors: F, SIGF).

        Handles multiple input formats:
        - IPAIR → FMEAN: French-Wilson via servalcat fw
        - IMEAN → FMEAN: French-Wilson via servalcat fw
        - FPAIR → FMEAN: Inverse-variance weighted mean via gemmi

        Args:
            obs_file: CObsDataFile instance to convert
            work_directory: Optional directory for working files

        Returns:
            Full path to converted FMEAN file

        Raises:
            CException: If validation fails, unsupported conversion, or fails
        """
        from pathlib import Path
        import shutil
        from ccp4i2.core.conversions.servalcat_converter import ServalcatConverter

        # Auto-detect content flag from file
        obs_file.setContentFlag()
        current_flag = int(obs_file.contentFlag)

        # Validation handles this now, but kept for clarity
        if current_flag == 0:
            raise CException(ObsDataConverter, 2,
                             details=f"File: {obs_file.getFullPath()}")

        # If already FMEAN, just copy the file
        if current_flag == obs_file.CONTENT_FLAG_FMEAN:
            output_path = obs_file._get_conversion_output_path(
                'FMEAN', work_directory=work_directory)
            input_path = Path(obs_file.getFullPath())
            shutil.copy2(input_path, output_path)
            return str(output_path)

        # FPAIR → FMEAN: Use gemmi-based weighted mean
        if current_flag == obs_file.CONTENT_FLAG_FPAIR:
            return ObsDataConverter._fpair_to_fmean_gemmi(
                obs_file, work_directory=work_directory)

        # IPAIR/IMEAN → FMEAN: Use servalcat fw
        return ServalcatConverter.to_fmean(obs_file,
                                           work_directory=work_directory)

    @staticmethod
    def _fpair_to_fmean_gemmi(obs_file, work_directory: Optional[Any] = None) -> str:
        """
        Convert FPAIR (F+/F-) to FMEAN (F, SIGF) using gemmi and numpy.

        This replaces the old sftools-based conversion with a pure Python implementation.
        Uses inverse-variance weighted averaging:
        - F = (F+ * SIGF-^2 + F- * SIGF+^2) / (SIGF+^2 + SIGF-^2)
        - SIGF = sqrt((SIGF+^2 * SIGF-^2) / (SIGF+^2 + SIGF-^2))

        For reflections with only F+ or F-, uses the available value directly.

        Args:
            obs_file: CObsDataFile instance with FPAIR data
            work_directory: Directory for output file

        Returns:
            Path to output MTZ file with F, SIGF columns
        """
        import gemmi
        import numpy as np

        # Read input MTZ
        input_path = obs_file.getFullPath()
        mtz = gemmi.read_mtz_file(str(input_path))

        # Get FPAIR column names from the file
        column_names = obs_file.CONTENT_SIGNATURE_LIST[obs_file.CONTENT_FLAG_FPAIR - 1]
        # column_names = ['Fplus', 'SIGFplus', 'Fminus', 'SIGFminus']

        # Find the actual columns (they might have different names in the file)
        fplus_col = mtz.column_with_label(column_names[0])
        sigfplus_col = mtz.column_with_label(column_names[1])
        fminus_col = mtz.column_with_label(column_names[2])
        sigfminus_col = mtz.column_with_label(column_names[3])

        # Extract data as numpy arrays
        fplus = np.array(fplus_col)
        sigfplus = np.array(sigfplus_col)
        fminus = np.array(fminus_col)
        sigfminus = np.array(sigfminus_col)

        # Initialize output arrays
        f_mean = np.full_like(fplus, np.nan)
        sigf_mean = np.full_like(fplus, np.nan)

        # Mask for valid data (not NaN)
        fplus_valid = ~np.isnan(fplus)
        fminus_valid = ~np.isnan(fminus)

        # Case 1: Both F+ and F- present - use weighted mean
        both_present = fplus_valid & fminus_valid
        if np.any(both_present):
            sigfsq_plus = sigfplus[both_present] ** 2
            sigfsq_minus = sigfminus[both_present] ** 2
            sigfsq_sum = sigfsq_plus + sigfsq_minus

            # Weighted mean: F = (F+ * SIG-^2 + F- * SIG+^2) / (SIG+^2 + SIG-^2)
            f_mean[both_present] = (fplus[both_present] * sigfsq_minus +
                                    fminus[both_present] * sigfsq_plus) / sigfsq_sum

            # Combined sigma: SIGF = sqrt((SIG+^2 * SIG-^2) / (SIG+^2 + SIG-^2))
            sigf_mean[both_present] = np.sqrt(
                (sigfsq_plus * sigfsq_minus) / sigfsq_sum)

        # Case 2: Only F+ present
        only_fplus = fplus_valid & ~fminus_valid
        if np.any(only_fplus):
            f_mean[only_fplus] = fplus[only_fplus]
            sigf_mean[only_fplus] = sigfplus[only_fplus]

        # Case 3: Only F- present
        only_fminus = ~fplus_valid & fminus_valid
        if np.any(only_fminus):
            f_mean[only_fminus] = fminus[only_fminus]
            sigf_mean[only_fminus] = sigfminus[only_fminus]

        # Create output MTZ with same structure as input
        mtz_out = gemmi.Mtz()
        mtz_out.spacegroup = mtz.spacegroup
        mtz_out.cell = mtz.cell

        # Add dataset and HKL columns
        mtz_out.add_dataset('crystal')
        mtz_out.add_column('H', 'H')
        mtz_out.add_column('K', 'H')
        mtz_out.add_column('L', 'H')

        # Add F and SIGF columns
        mtz_out.add_column('F', 'F')
        mtz_out.add_column('SIGF', 'Q')

        # Get HKL indices from input and create data array
        miller_array = mtz.make_miller_array()

        # Create data array: [H, K, L, F, SIGF]
        data = np.column_stack([
            miller_array,
            f_mean,
            sigf_mean
        ])

        mtz_out.set_data(data)

        # Determine output path
        output_path = obs_file._get_conversion_output_path(
            'FMEAN', work_directory=work_directory)

        # Write output file
        mtz_out.write_to_file(str(output_path))

        print(f"✅ FPAIR → FMEAN conversion completed: {output_path}")
        return str(output_path)
