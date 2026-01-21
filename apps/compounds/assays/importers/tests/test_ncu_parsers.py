"""
Tests for NCU ADME parsers.

Uses fixture files in ./fixtures/ directory.
"""

import pytest
from pathlib import Path

from ..registry import detect_parser, parse_adme_file
from ..ncu.lm import LiverMicrosomeParser
from ..ncu.bs import BloodSerumStabilityParser
from ..ncu.gsh import GSHStabilityParser
from ..ncu.caco2 import Caco2PermeabilityParser


FIXTURES_DIR = Path(__file__).parent / 'fixtures'


class TestParserDetection:
    """Test that parsers are correctly detected for each file type."""

    def test_detect_lm_parser(self):
        """LM file should be detected by LiverMicrosomeParser."""
        filepath = FIXTURES_DIR / 'ADME-NCU-LM-20231218.xlsx'
        parser = detect_parser(filepath)
        assert parser is not None, f"No parser detected for {filepath.name}"
        assert isinstance(parser, LiverMicrosomeParser)

    def test_detect_bs_parser(self):
        """BS file should be detected by BloodSerumStabilityParser."""
        filepath = FIXTURES_DIR / 'ADME-NCU-BS-20240906.xlsx'
        parser = detect_parser(filepath)
        assert parser is not None, f"No parser detected for {filepath.name}"
        assert isinstance(parser, BloodSerumStabilityParser)

    def test_detect_gsh_parser(self):
        """GSH file should be detected by GSHStabilityParser."""
        filepath = FIXTURES_DIR / 'ADME-NCU-GSH Stability-20241106.xlsx'
        parser = detect_parser(filepath)
        assert parser is not None, f"No parser detected for {filepath.name}"
        assert isinstance(parser, GSHStabilityParser)

    def test_detect_caco2_parser(self):
        """Caco-2 file should be detected by Caco2PermeabilityParser."""
        filepath = FIXTURES_DIR / 'ADME-NCU-Caco-2 Permeability-20231219.xlsx'
        parser = detect_parser(filepath)
        assert parser is not None, f"No parser detected for {filepath.name}"
        assert isinstance(parser, Caco2PermeabilityParser)


class TestLiverMicrosomeParser:
    """Tests for NCU Liver Microsome (LM) parser."""

    @pytest.fixture
    def filepath(self):
        return FIXTURES_DIR / 'ADME-NCU-LM-20231218.xlsx'

    @pytest.fixture
    def parser(self):
        return LiverMicrosomeParser()

    def test_parse_returns_results(self, parser, filepath):
        """Parser should return ParseResult with compounds."""
        result = parser.parse(filepath)

        # Should have no critical errors
        critical_errors = [e for e in result.errors if e.severity == 'error']
        assert not critical_errors, f"Critical errors: {[e.message for e in critical_errors]}"

        # Should have results
        assert result.results, "No results parsed"
        print(f"\nLM: Parsed {len(result.results)} results")
        for r in result.results:
            print(f"  - {r.compound_id} ({r.species}): {r.results.get('clint_ul_min_mg')}")

    def test_parse_extracts_compound_ids(self, parser, filepath):
        """Parser should extract compound IDs correctly."""
        result = parser.parse(filepath)

        compound_ids = {r.compound_id for r in result.results}
        # Should have the test compound and control
        assert 'PH-NTU-00-00-000-0' in compound_ids or 'Verapamil' in compound_ids, \
            f"Expected compound IDs not found. Got: {compound_ids}"

    def test_parse_extracts_species(self, parser, filepath):
        """Parser should extract species (Human, Mouse) correctly."""
        result = parser.parse(filepath)

        species_set = {r.species for r in result.results}
        assert 'Human' in species_set or 'Mouse' in species_set, \
            f"Expected species not found. Got: {species_set}"

    def test_parse_extracts_clint(self, parser, filepath):
        """Parser should extract CLint values."""
        result = parser.parse(filepath)

        # At least one result should have CLint
        clint_values = [r.results.get('clint_ul_min_mg') for r in result.results]
        assert any(v is not None for v in clint_values), \
            f"No CLint values found. Got: {clint_values}"

    def test_metadata_contains_sheet_info(self, parser, filepath):
        """Metadata should contain sheet information."""
        result = parser.parse(filepath)

        assert 'sheets_found' in result.metadata
        assert 'summary_sheet_used' in result.metadata
        assert 'Summary' in result.metadata['sheets_found']


class TestBloodSerumStabilityParser:
    """Tests for NCU Blood/Serum Stability (BS) parser."""

    @pytest.fixture
    def filepath(self):
        return FIXTURES_DIR / 'ADME-NCU-BS-20240906.xlsx'

    @pytest.fixture
    def parser(self):
        return BloodSerumStabilityParser()

    def test_parse_returns_results(self, parser, filepath):
        """Parser should return ParseResult with compounds."""
        result = parser.parse(filepath)

        # Should have no critical errors
        critical_errors = [e for e in result.errors if e.severity == 'error']
        assert not critical_errors, f"Critical errors: {[e.message for e in critical_errors]}"

        # Should have results
        assert result.results, "No results parsed"
        print(f"\nBS: Parsed {len(result.results)} results")
        for r in result.results:
            print(f"  - {r.compound_id} ({r.species}): t1/2={r.results.get('t1_2_min')}")

    def test_parse_extracts_compound_ids(self, parser, filepath):
        """Parser should extract NCL compound IDs."""
        result = parser.parse(filepath)

        compound_ids = {r.compound_id for r in result.results}
        # Should have NCL compounds (anonymized to NCL-0000000X)
        ncl_ids = [cid for cid in compound_ids if cid.startswith('NCL-')]
        assert ncl_ids, f"No NCL compound IDs found. Got: {compound_ids}"

    def test_parse_extracts_t1_2(self, parser, filepath):
        """Parser should extract t1/2 values."""
        result = parser.parse(filepath)

        t12_values = [r.results.get('t1_2_min') for r in result.results]
        assert any(v is not None for v in t12_values), \
            f"No t1/2 values found"

    def test_parse_extracts_time_course(self, parser, filepath):
        """Parser should extract time course data."""
        result = parser.parse(filepath)

        # At least one result should have time course
        has_time_course = any('time_course' in r.results for r in result.results)
        assert has_time_course, "No time course data found"


class TestGSHStabilityParser:
    """Tests for NCU GSH Stability parser."""

    @pytest.fixture
    def filepath(self):
        return FIXTURES_DIR / 'ADME-NCU-GSH Stability-20241106.xlsx'

    @pytest.fixture
    def parser(self):
        return GSHStabilityParser()

    def test_parse_returns_results(self, parser, filepath):
        """Parser should return ParseResult with compounds."""
        result = parser.parse(filepath)

        # Should have no critical errors
        critical_errors = [e for e in result.errors if e.severity == 'error']
        assert not critical_errors, f"Critical errors: {[e.message for e in critical_errors]}"

        # Should have results
        assert result.results, "No results parsed"
        print(f"\nGSH: Parsed {len(result.results)} results")
        for r in result.results:
            print(f"  - {r.compound_id}: t1/2={r.results.get('t1_2_min')}, gsh_reactive={r.results.get('gsh_reactive')}")

    def test_parse_extracts_compound_ids(self, parser, filepath):
        """Parser should extract NCL compound IDs."""
        result = parser.parse(filepath)

        compound_ids = {r.compound_id for r in result.results}
        ncl_ids = [cid for cid in compound_ids if cid.startswith('NCL-')]
        assert ncl_ids, f"No NCL compound IDs found. Got: {compound_ids}"

    def test_parse_extracts_gsh_conditions(self, parser, filepath):
        """Parser should extract both +GSH and -GSH conditions."""
        result = parser.parse(filepath)

        # Check for time_course with both conditions
        for r in result.results:
            if 'time_course' in r.results:
                tc = r.results['time_course']
                assert 'with_gsh' in tc or 'without_gsh' in tc, \
                    f"Expected GSH conditions in time_course"


class TestCaco2PermeabilityParser:
    """Tests for NCU Caco-2 Permeability parser."""

    @pytest.fixture
    def filepath(self):
        return FIXTURES_DIR / 'ADME-NCU-Caco-2 Permeability-20231219.xlsx'

    @pytest.fixture
    def parser(self):
        return Caco2PermeabilityParser()

    def test_parse_returns_results(self, parser, filepath):
        """Parser should return ParseResult with compounds."""
        result = parser.parse(filepath)

        # Should have no critical errors
        critical_errors = [e for e in result.errors if e.severity == 'error']
        assert not critical_errors, f"Critical errors: {[e.message for e in critical_errors]}"

        # Should have results
        assert result.results, "No results parsed"
        print(f"\nCaco-2: Parsed {len(result.results)} results")
        for r in result.results:
            print(f"  - {r.compound_id}: Papp A-B={r.results.get('papp_ab')}, Efflux={r.results.get('efflux_ratio')}")

    def test_parse_extracts_compound_ids(self, parser, filepath):
        """Parser should extract compound IDs."""
        result = parser.parse(filepath)

        compound_ids = {r.compound_id for r in result.results}
        # Should have test compound and controls
        assert compound_ids, f"No compound IDs found"

    def test_parse_extracts_permeability_values(self, parser, filepath):
        """Parser should extract Papp and efflux ratio values."""
        result = parser.parse(filepath)

        papp_values = [r.results.get('papp_ab') for r in result.results]
        efflux_values = [r.results.get('efflux_ratio') for r in result.results]

        assert any(v is not None for v in papp_values), "No Papp A-B values found"
        assert any(v is not None for v in efflux_values), "No efflux ratio values found"

    def test_parse_extracts_monolayer_integrity(self, parser, filepath):
        """Parser should extract TEER and LY leakage data."""
        result = parser.parse(filepath)

        has_integrity = any('monolayer_integrity' in r.results for r in result.results)
        assert has_integrity, "No monolayer integrity data found"


class TestParseAdmeFile:
    """Test the high-level parse_adme_file function."""

    def test_parse_bs_file(self):
        """parse_adme_file should work for BS files."""
        filepath = FIXTURES_DIR / 'ADME-NCU-BS-20240906.xlsx'
        result = parse_adme_file(filepath)

        assert result.results, f"No results. Errors: {[e.message for e in result.errors]}"

    def test_parse_lm_file(self):
        """parse_adme_file should work for LM files."""
        filepath = FIXTURES_DIR / 'ADME-NCU-LM-20231218.xlsx'
        result = parse_adme_file(filepath)

        assert result.results, f"No results. Errors: {[e.message for e in result.errors]}"

    def test_parse_gsh_file(self):
        """parse_adme_file should work for GSH files."""
        filepath = FIXTURES_DIR / 'ADME-NCU-GSH Stability-20241106.xlsx'
        result = parse_adme_file(filepath)

        assert result.results, f"No results. Errors: {[e.message for e in result.errors]}"

    def test_parse_caco2_file(self):
        """parse_adme_file should work for Caco-2 files."""
        filepath = FIXTURES_DIR / 'ADME-NCU-Caco-2 Permeability-20231219.xlsx'
        result = parse_adme_file(filepath)

        assert result.results, f"No results. Errors: {[e.message for e in result.errors]}"
