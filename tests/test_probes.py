"""Tests for source-neutral probe-property containers."""

import math

import pytest

from igc_sea.analysis.probes import ProbeError, ProbeProperties


def _properties(**overrides) -> ProbeProperties:
    values = {
        "name": "SYNTHETIC PROBE",
        "cross_section_m2": 6.3e-19,
        "antoine_c1": 96.084,
        "antoine_c2": -7900.2,
        "antoine_c3": -11.0,
        "antoine_c4": 7.0e-6,
        "antoine_c5": 2.0,
    }
    values.update(overrides)
    return ProbeProperties(**values)


def test_antoine_mapping_is_explicit():
    properties = _properties()
    assert properties.antoine == {
        "C1": 96.084,
        "C2": -7900.2,
        "C3": -11.0,
        "C4": 7.0e-6,
        "C5": 2.0,
    }


def test_positive_saturation_pressure():
    pressure = _properties().p_sat(303.15)
    assert math.isfinite(pressure)
    assert pressure > 0


def test_valid_properties_pass_validation():
    _properties().validate()


@pytest.mark.parametrize("cross_section", [0.0, -1.0, float("nan")])
def test_invalid_cross_section_is_rejected(cross_section):
    with pytest.raises(ProbeError, match="cross_section_m2"):
        _properties(cross_section_m2=cross_section).validate()


def test_nonfinite_coefficient_is_rejected():
    with pytest.raises(ProbeError, match="coefficient C1"):
        _properties(antoine_c1=float("nan")).validate()
