import pytest
from pyomo.environ import value, units as pyunits
from pyomo.util.check_units import assert_units_consistent
from wrd.components.ro import main, default_ro_config
from watertap.unit_models.reverse_osmosis_1D import (
    ReverseOsmosis1D,
    PressureChangeType,
    MassTransferCoefficient,
    ConcentrationPolarizationType,
)

test_default_ro_config = dict(
    has_pressure_change=True,
    pressure_change_type=PressureChangeType.fixed_per_stage,
    mass_transfer_coefficient=MassTransferCoefficient.calculated,
    concentration_polarization_type=ConcentrationPolarizationType.calculated,
    transformation_scheme="BACKWARD",
    transformation_method="dae.finite_difference",
    module_type="spiral_wound",
    finite_elements=10,
    has_full_reporting=True,
)


@pytest.mark.component
def test_ro_main():
    m = main()
    assert isinstance(m.fs.ro.unit, ReverseOsmosis1D)


@pytest.mark.component
def test_ro_default_config():
    assert len(default_ro_config) == len(test_default_ro_config)


import pprint

pprint.pprint(default_ro_config)
pprint.pprint(test_default_ro_config)
