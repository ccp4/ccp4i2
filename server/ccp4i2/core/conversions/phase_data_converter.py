"""
Converter for phase data format transformations.

This module handles conversions between different phase data formats in MTZ files:

Phase Data Files (CPhsDataFile):
- HL (1): Hendrickson-Lattman coefficients (HLA, HLB, HLC, HLD)
- PHIFOM (2): Phase + Figure of Merit (PHI, FOM)

Map Coefficients Files (CMapCoeffsDataFile):
- FPHI (1): Structure factors + Phase (F, PHI)

Conversion Matrix (CPhsDataFile):
                TO
          HL   PHIFOM
FROM HL    ✓      ✓
     PHIFOM ✓      ✓

Production Implementation:
- Uses CCP4 chltofom plugin for accurate, validated HL ↔ PHIFOM conversions
- Integrated with CCP4i2 plugin framework for consistent error handling
- Ensures compatibility with CCP4 ecosystem
- Provides reliable, well-tested transformations
- Automatic direction detection based on input contentFlag

The module also includes experimental gemmi-based implementations for future
development, but these are not used in production due to FOM calculation
differences compared to the CCP4 reference.

Dependencies:
- CCP4 chltofom (required for production conversions)
- gemmi (used for MTZ I/O and experimental implementations)
- numpy (used for experimental numerical calculations)

References:
- Read, R.J. (1986). Acta Cryst. A42, 140-149. (Phase probability distributions)
- Hendrickson & Lattman (1970). Acta Cryst. B26, 136-143. (HL coefficients)
"""

from typing import Optional, Any
from ccp4i2.core.CCP4ErrorHandling import CException, CErrorReport, SEVERITY_ERROR


class PhaseDataConverter:
    """
    Static converter class for phase data format transformations.

    Handles conversions between HL coefficients, PHI/FOM, and F/PHI formats.
    All methods are static and take the file instance as first parameter.
    """

    # Error codes for phase data conversions
    ERROR_CODES = {
        1: {'description': 'Input file does not exist', 'severity': SEVERITY_ERROR},
        2: {'description': 'Cannot determine contentFlag from input file', 'severity': SEVERITY_ERROR},
        3: {'description': 'Unsupported conversion: only HL and PHIFOM formats supported', 'severity': SEVERITY_ERROR},
        4: {'description': 'chltofom plugin not available - ensure CCP4I2_ROOT is set', 'severity': SEVERITY_ERROR},
        5: {'description': 'Cannot determine conversion direction from input columns', 'severity': SEVERITY_ERROR},
        6: {'description': 'chltofom completed but output file not created', 'severity': SEVERITY_ERROR},
        7: {'description': 'Expected column not found in chltofom output', 'severity': SEVERITY_ERROR},
        8: {'description': 'PHI column not found in input MTZ file', 'severity': SEVERITY_ERROR},
        9: {'description': 'FOM column not found in input MTZ file', 'severity': SEVERITY_ERROR},
        10: {'description': 'Invalid number of HL coefficient columns', 'severity': SEVERITY_ERROR},
    }

    @staticmethod
    def _validate_input_file(phase_file):
        """
        Validate input file before conversion.

        Args:
            phase_file: CPhsDataFile or CMapCoeffsDataFile instance

        Raises:
            CException: If file doesn't exist or contentFlag cannot be determined
        """
        from pathlib import Path

        input_path = phase_file.getFullPath()

        # Check file exists
        if not Path(input_path).exists():
            raise CException(PhaseDataConverter, 1, details=f"File: {input_path}")

        # Check contentFlag is set
        content_flag = int(phase_file.contentFlag) if phase_file.contentFlag.isSet() else 0
        if content_flag == 0:
            raise CException(PhaseDataConverter, 2, details=f"File: {input_path}")

    @staticmethod
    def to_hl(phase_file, work_directory: Optional[Any] = None) -> str:
        """
        Convert phase data to HL format (Hendrickson-Lattman coefficients).

        HL format: HLA, HLB, HLC, HLD

        Possible conversions:
        - PHIFOM → HL: Convert phase+FOM to HL coefficients using CCP4 chltofom

        Uses the CCP4 chltofom tool for accurate, validated conversion.

        Args:
            phase_file: CPhsDataFile or CMapCoeffsDataFile instance
            work_directory: Directory for output files

        Returns:
            Full path to converted HL file

        Raises:
            CException: If validation fails, conversion not supported, or chltofom fails
        """
        # Validate input file
        PhaseDataConverter._validate_input_file(phase_file)

        output_path = phase_file._get_conversion_output_path('HL', work_directory=work_directory)
        input_path = phase_file.getFullPath()

        # Detect current content flag and route to appropriate conversion
        content_flag = int(phase_file.contentFlag)

        if content_flag == 1:  # HL format
            # Already in HL format, just return current path
            return input_path
        elif content_flag == 2:  # PHIFOM format
            return PhaseDataConverter._phifom_to_hl_chltofom(
                input_path, output_path)
        else:
            raise CException(
                PhaseDataConverter, 3,
                details=f"contentFlag={content_flag}, target=HL. Only PHIFOM (2) → HL supported."
            )

    @staticmethod
    def to_phifom(phase_file, work_directory: Optional[Any] = None) -> str:
        """
        Convert phase data to PHIFOM format (Phase + Figure of Merit).

        PHIFOM format: PHI, FOM

        Possible conversions:
        - HL → PHIFOM: Convert HL coefficients to best phase estimate + FOM using CCP4 chltofom

        Uses the CCP4 chltofom tool for accurate, validated conversion.

        Args:
            phase_file: CPhsDataFile or CMapCoeffsDataFile instance
            work_directory: Directory for output files

        Returns:
            Full path to converted PHIFOM file

        Raises:
            CException: If validation fails, conversion not supported, or chltofom fails
        """
        # Validate input file
        PhaseDataConverter._validate_input_file(phase_file)

        output_path = phase_file._get_conversion_output_path('PHIFOM', work_directory=work_directory)
        input_path = phase_file.getFullPath()

        # Detect current content flag and route to appropriate conversion
        content_flag = int(phase_file.contentFlag)

        if content_flag == 1:  # HL format
            return PhaseDataConverter._hl_to_phifom_chltofom(
                input_path, output_path)
        elif content_flag == 2:  # PHIFOM format
            # Already in PHIFOM format, just return current path
            return input_path
        else:
            raise CException(
                PhaseDataConverter, 3,
                details=f"contentFlag={content_flag}, target=PHIFOM. Only HL (1) → PHIFOM supported."
            )

    @staticmethod
    def to_fphi(map_coeffs_file, work_directory: Optional[Any] = None) -> str:
        """
        Return path to FPHI format file.

        Since CMapCoeffsDataFile only has one content type (FPHI),
        this method simply returns the current file path without conversion.

        FPHI format: F, PHI

        Args:
            map_coeffs_file: CMapCoeffsDataFile instance
            work_directory: Ignored for this conversion

        Returns:
            Full path to current file (no conversion needed)
        """
        # No conversion needed - CMapCoeffsDataFile is always FPHI
        return map_coeffs_file.getFullPath()

    # ========================================================================
    # Production implementation using CCP4 chltofom plugin
    # ========================================================================

    @staticmethod
    def _run_chltofom_plugin(input_path: str, output_path: str) -> str:
        """
        Run chltofom plugin for HL ↔ PHIFOM conversion.

        The plugin automatically detects conversion direction based on
        the input file's contentFlag. After running chltofom, this method
        extracts only the converted columns (plus H, K, L) to ensure the
        output file has the correct contentFlag signature.

        Args:
            input_path: Path to input MTZ (either HL or PHIFOM format)
            output_path: Path for output MTZ

        Returns:
            Path to output MTZ file with only converted columns

        Raises:
            RuntimeError: If chltofom plugin fails
        """
        from pathlib import Path
        import gemmi

        # Import the chltofom plugin wrapper
        from ccp4i2.wrappers.chltofom.script import chltofom

        # Detect input format to determine which columns to extract
        input_mtz = gemmi.read_mtz_file(input_path)
        input_labels = [col.label for col in input_mtz.columns]

        # Determine conversion direction
        has_hl = all(col in input_labels for col in ['HLA', 'HLB', 'HLC', 'HLD'])
        has_phifom = all(col in input_labels for col in ['PHI', 'FOM'])

        if has_hl and not has_phifom:
            # HL → PHIFOM
            columns_to_keep = ['PHI', 'FOM']
        elif has_phifom and not has_hl:
            # PHIFOM → HL
            columns_to_keep = ['HLA', 'HLB', 'HLC', 'HLD']
        else:
            raise CException(
                PhaseDataConverter, 5,
                details=f"Input columns: {', '.join(input_labels)}"
            )

        # Create the plugin wrapper with correct workDirectory
        # Use output_path's directory to avoid polluting project root with log files
        work_dir = Path(output_path).parent
        wrapper = chltofom.chltofom(name='phase_conversion', workDirectory=str(work_dir))

        # Set up input
        wrapper.container.inputData.HKLIN.setFullPath(input_path)

        # Set up output - write directly to output path
        # (we'll post-process to extract only desired columns)
        wrapper.container.outputData.HKLOUT.setFullPath(output_path)

        # Note: Don't use OUTPUTMINIMTZ as it requires CCP4Utils.makeTmpFile
        wrapper.container.controlParameters.OUTPUTMINIMTZ.set(False)

        # Run the plugin
        try:
            from ccp4i2.core.CCP4Modules import PROCESSMANAGER
            PROCESSMANAGER().setWaitForFinished(10000)  # 10 second timeout per step
        except Exception as e:
            # PROCESSMANAGER not available in minimal environment - continue anyway
            print(f"Warning: Could not set PROCESSMANAGER timeout: {e}")

        pid = wrapper.process()

        # Check if output was created
        if not Path(output_path).exists():
            # Check if plugin reported errors
            error_details = f"Output file: {output_path}"
            if hasattr(wrapper, 'errorReport') and wrapper.errorReport.count() > 0:
                error_details += f"\nPlugin errors: {wrapper.errorReport.report()}"
            raise CException(PhaseDataConverter, 6, details=error_details)

        # Now extract only the desired columns using gemmi
        # Read the file that has all columns
        mtz_full = gemmi.read_mtz_file(output_path)
        mtz_mini = gemmi.Mtz()

        # Copy crystal and dataset info
        mtz_mini.cell = mtz_full.cell
        mtz_mini.spacegroup = mtz_full.spacegroup
        mtz_mini.add_dataset('HKL_base')

        # Collect columns to copy (H, K, L + converted columns)
        columns_to_copy = []

        # Copy H, K, L columns
        for col_label in ['H', 'K', 'L']:
            col = mtz_full.column_with_label(col_label)
            mtz_mini.add_column(col_label, col.type)
            columns_to_copy.append(col)

        # Copy converted columns
        for col_label in columns_to_keep:
            # Find the column (may have a complex label from chltofom)
            col = None
            for c in mtz_full.columns:
                if c.label == col_label or c.label.startswith(col_label):
                    col = c
                    break

            if col is None:
                available_cols = [c.label for c in mtz_full.columns]
                raise CException(
                    PhaseDataConverter, 7,
                    details=f"Column '{col_label}' not in output. Available: {', '.join(available_cols)}"
                )

            # Add column with simple label (not the complex chltofom label)
            mtz_mini.add_column(col_label, col.type)
            columns_to_copy.append(col)

        # Copy data row by row
        import numpy as np
        n_reflections = mtz_full.nreflections
        data = np.zeros((n_reflections, len(columns_to_copy)), dtype=np.float32)

        for i, col in enumerate(columns_to_copy):
            data[:, i] = np.array(list(col), dtype=np.float32)

        mtz_mini.set_data(data)

        # Overwrite the output file with mini version
        mtz_mini.write_to_file(output_path)

        return output_path

    @staticmethod
    def _hl_to_phifom_chltofom(input_path: str, output_path: str) -> str:
        """
        Convert HL coefficients to PHIFOM format using CCP4 chltofom plugin.

        This is the production implementation that uses the validated
        CCP4 plugin framework for accurate conversion.

        Args:
            input_path: Path to input MTZ with HL coefficients
            output_path: Path for output MTZ with PHI/FOM

        Returns:
            Path to output MTZ file

        Raises:
            RuntimeError: If chltofom plugin fails
        """
        return PhaseDataConverter._run_chltofom_plugin(input_path, output_path)

    @staticmethod
    def _phifom_to_hl_chltofom(input_path: str, output_path: str) -> str:
        """
        Convert PHIFOM format to HL coefficients using CCP4 chltofom plugin.

        This is the production implementation that uses the validated
        CCP4 plugin framework for accurate conversion.

        Args:
            input_path: Path to input MTZ with PHI/FOM
            output_path: Path for output MTZ with HL coefficients

        Returns:
            Path to output MTZ file

        Raises:
            RuntimeError: If chltofom plugin fails
        """
        return PhaseDataConverter._run_chltofom_plugin(input_path, output_path)

    # ========================================================================
    # Experimental gemmi-based implementation (for future development)
    # ========================================================================
    #
    # NOTE: These methods implement HL ↔ PHIFOM conversions using pure Python
    # with gemmi and numpy. However, benchmark testing against CCP4 chltofom
    # shows FOM correlation of only 0.08 (expected >0.99), indicating our
    # numerical method for calculating FOM from the phase probability
    # distribution differs significantly from the CCP4 reference.
    #
    # The gemmi-based approach is preserved here for future investigation and
    # refinement. Possible issues to investigate:
    # - Phase probability distribution calculation
    # - Numerical integration method for FOM
    # - Handling of HLC/HLD terms in the calculation
    #
    # For production use, the chltofom-based methods above should be used.
    # ========================================================================

    @staticmethod
    def _phifom_to_hl_gemmi_experimental(input_path: str, output_path: str, mtz) -> str:
        """
        Convert PHIFOM format to HL coefficients using gemmi.

        Uses centrosymmetric approximation:
        - HLA = FOM * cos(PHI)
        - HLB = FOM * sin(PHI)
        - HLC = 0
        - HLD = 0

        Args:
            input_path: Path to input MTZ with PHI/FOM
            output_path: Path for output MTZ with HL coefficients
            mtz: gemmi.Mtz object (already loaded)

        Returns:
            Path to output MTZ file
        """
        import gemmi
        import numpy as np

        # Find PHI and FOM columns
        phi_col = None
        fom_col = None

        for col in mtz.columns:
            label = col.label.upper()
            if col.type == 'P' or 'PHI' in label:  # Phase column
                phi_col = col
            elif col.type == 'W' or 'FOM' in label:  # Weight/FOM column
                fom_col = col

        if phi_col is None:
            raise CException(PhaseDataConverter, 8, details=f"Input file: {input_path}")
        if fom_col is None:
            raise CException(PhaseDataConverter, 9, details=f"Input file: {input_path}")

        # Extract PHI and FOM as numpy arrays
        phi = np.array(phi_col)  # Degrees
        fom = np.array(fom_col)

        # Convert PHI to radians for calculation
        phi_rad = np.radians(phi)

        # Calculate HL coefficients using centrosymmetric approximation
        hla = fom * np.cos(phi_rad)
        hlb = fom * np.sin(phi_rad)
        hlc = np.zeros_like(fom)
        hld = np.zeros_like(fom)

        # Create output MTZ with same structure
        out_mtz = gemmi.Mtz(with_base=False)
        out_mtz.spacegroup = mtz.spacegroup
        out_mtz.set_cell_for_all(mtz.cell)

        # Add a single dataset for all columns
        out_mtz.add_dataset('hl_data')

        # Add H, K, L columns
        out_mtz.add_column('H', 'H')
        out_mtz.add_column('K', 'H')
        out_mtz.add_column('L', 'H')

        # Add HL coefficient columns
        out_mtz.add_column('HLA', 'A')  # HL coefficient type
        out_mtz.add_column('HLB', 'A')
        out_mtz.add_column('HLC', 'A')
        out_mtz.add_column('HLD', 'A')

        # Set data
        out_mtz.set_data(np.column_stack([
            np.array(mtz.column_with_label('H')),
            np.array(mtz.column_with_label('K')),
            np.array(mtz.column_with_label('L')),
            hla,
            hlb,
            hlc,
            hld
        ]))

        # Write output
        out_mtz.write_to_file(output_path)
        return output_path

    @staticmethod
    def _hl_to_phifom_gemmi_experimental(input_path: str, output_path: str, mtz) -> str:
        """
        EXPERIMENTAL: Convert HL coefficients to PHIFOM format using gemmi.

        Reads HLA, HLB, HLC, HLD columns and calculates best phase (PHI) and
        figure of merit (FOM) using numerical optimization of the phase
        probability distribution.

        Args:
            input_path: Path to input MTZ with HL coefficients
            output_path: Path for output MTZ with PHI/FOM
            mtz: gemmi.Mtz object (already loaded)

        Returns:
            Path to output MTZ file

        References:
            - Read, R.J. (1986). Acta Cryst. A42, 140-149.
            - Hendrickson & Lattman (1970). Acta Cryst. B26, 136-143.
        """
        import gemmi
        import numpy as np

        # Find HL coefficient columns (column type 'A')
        hl_cols = {}
        for col in mtz.columns:
            if col.type == 'A':  # HL coefficient type
                label = col.label
                if 'HLA' in label:
                    hl_cols['HLA'] = col
                elif 'HLB' in label:
                    hl_cols['HLB'] = col
                elif 'HLC' in label:
                    hl_cols['HLC'] = col
                elif 'HLD' in label:
                    hl_cols['HLD'] = col

        if len(hl_cols) != 4:
            raise CException(
                PhaseDataConverter, 10,
                details=f"Found {len(hl_cols)} columns: {', '.join(hl_cols.keys())}"
            )

        # Extract HL coefficients as numpy arrays
        hla = np.array(hl_cols['HLA'])
        hlb = np.array(hl_cols['HLB'])
        hlc = np.array(hl_cols['HLC'])
        hld = np.array(hl_cols['HLD'])

        # Calculate PHI and FOM from HL coefficients
        phi, fom = PhaseDataConverter._hl_to_phifom_calculation(hla, hlb, hlc, hld)

        # Create output MTZ with same structure
        out_mtz = gemmi.Mtz(with_base=False)
        out_mtz.spacegroup = mtz.spacegroup
        out_mtz.set_cell_for_all(mtz.cell)

        # Add a single dataset for all columns
        out_mtz.add_dataset('phase_data')

        # Add H, K, L columns
        out_mtz.add_column('H', 'H')
        out_mtz.add_column('K', 'H')
        out_mtz.add_column('L', 'H')

        # Add PHI and FOM columns
        out_mtz.add_column('PHI', 'P')  # Phase column type
        out_mtz.add_column('FOM', 'W')  # Weight/FOM column type

        # Set data
        out_mtz.set_data(np.column_stack([
            np.array(mtz.column_with_label('H')),
            np.array(mtz.column_with_label('K')),
            np.array(mtz.column_with_label('L')),
            phi,
            fom
        ]))

        # Write output
        out_mtz.write_to_file(output_path)
        return output_path

    @staticmethod
    def _hl_to_phifom_calculation(hla, hlb, hlc, hld):
        """
        Calculate best phase and FOM from Hendrickson-Lattman coefficients.

        The phase probability distribution P(phi) is proportional to:
        exp(A*cos(phi) + B*sin(phi) + C*cos(2*phi) + D*sin(2*phi))

        where A=HLA, B=HLB, C=HLC, D=HLD

        Best phase (phi_best) is found by maximizing P(phi) numerically.
        FOM = <cos(phi - phi_best)> is calculated by numerical integration.

        Args:
            hla, hlb, hlc, hld: Hendrickson-Lattman coefficient arrays (numpy)

        Returns:
            tuple: (phi_best, fom) - Best phase estimate and figure of merit
                   Both in degrees for phi_best, FOM in [0, 1]

        References:
            - Read, R.J. (1986). Acta Cryst. A42, 140-149.
            - Hendrickson & Lattman (1970). Acta Cryst. B26, 136-143.
        """
        import numpy as np

        n_reflections = len(hla)
        phi_best = np.zeros(n_reflections)
        fom = np.zeros(n_reflections)

        # Sample phases from 0 to 360 degrees
        # Use fine grid for accurate maximum finding
        phi_samples = np.linspace(0, 2 * np.pi, 360)

        for i in range(n_reflections):
            # Calculate phase probability distribution
            # P(phi) ∝ exp(A*cos(phi) + B*sin(phi) + C*cos(2*phi) + D*sin(2*phi))
            exponent = (hla[i] * np.cos(phi_samples) +
                        hlb[i] * np.sin(phi_samples) +
                        hlc[i] * np.cos(2 * phi_samples) +
                        hld[i] * np.sin(2 * phi_samples))

            # Prevent overflow by subtracting maximum
            exponent_max = np.max(exponent)
            prob = np.exp(exponent - exponent_max)

            # Find best phase (maximum probability)
            max_idx = np.argmax(prob)
            phi_best[i] = phi_samples[max_idx]

            # Calculate FOM = <cos(phi - phi_best)> by numerical integration
            # FOM = ∫ cos(phi - phi_best) * P(phi) dphi / ∫ P(phi) dphi
            cos_term = np.cos(phi_samples - phi_best[i])
            numerator = np.trapezoid(cos_term * prob, phi_samples)
            denominator = np.trapezoid(prob, phi_samples)

            fom[i] = numerator / denominator if denominator > 0 else 0.0

            # Ensure FOM is in valid range [0, 1]
            fom[i] = np.clip(fom[i], 0.0, 1.0)

        # Convert phi_best from radians to degrees
        phi_best_deg = np.degrees(phi_best)

        return phi_best_deg, fom
