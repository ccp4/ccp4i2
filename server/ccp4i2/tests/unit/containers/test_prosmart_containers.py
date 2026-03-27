#!/usr/bin/env python3
"""
Test that prosmart_refmac loads all containers from .def.xml
"""
from ccp4i2.pipelines.prosmart_refmac.script.prosmart_refmac import prosmart_refmac


def test_prosmart_refmac_has_standard_containers():
    """Test that prosmart_refmac plugin has standard container sections."""
    plugin = prosmart_refmac()

    assert hasattr(plugin.container, 'inputData'), "Missing inputData container"
    assert hasattr(plugin.container, 'outputData'), "Missing outputData container"
    assert hasattr(plugin.container, 'controlParameters'), "Missing controlParameters container"


def test_prosmart_refmac_has_prosmart_protein():
    """Test that prosmart_refmac plugin has the prosmartProtein sub-container."""
    plugin = prosmart_refmac()

    assert hasattr(plugin.container, 'prosmartProtein'), \
        "prosmartProtein container does not exist"

    assert hasattr(plugin.container.prosmartProtein, 'TOGGLE'), \
        "TOGGLE not found in prosmartProtein"


def test_prosmart_refmac_has_specific_containers():
    """Test that prosmart_refmac plugin has task-specific sub-containers."""
    plugin = prosmart_refmac()

    assert hasattr(plugin.container, 'prosmartProtein'), "Missing prosmartProtein"
    assert hasattr(plugin.container, 'prosmartNucleicAcid'), "Missing prosmartNucleicAcid"
    assert hasattr(plugin.container, 'platonyzer'), "Missing platonyzer"
