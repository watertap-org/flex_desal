import pytest

from pyomo.environ import (
    ConcreteModel,
    TransformationFactory,
    assert_optimal_termination,
    units as pyunits,
)
from pyomo.network import Arc

from idaes.core import FlowsheetBlock
import idaes.core.util.scaling as iscale
from idaes.core.util.testing import initialization_tester
from idaes.core.util.initialization import propagate_state
from idaes.models.unit_models import Feed, Product
from idaes.core.util.model_statistics import degrees_of_freedom

from watertap.property_models.seawater_prop_pack import SeawaterParameterBlock
from watertap.core.control_volume_isothermal import ControlVolume0DBlock
from watertap.core.solvers import get_solver

from models.head_loss import HeadLoss

solver = get_solver()


def build_head_loss():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = SeawaterParameterBlock()
    m.fs.feed = Feed(property_package=m.fs.properties)
    m.fs.unit = HeadLoss(property_package=m.fs.properties)
    m.fs.product = Product(property_package=m.fs.properties)

    m.fs.feed_to_unit = Arc(source=m.fs.feed.outlet, destination=m.fs.unit.inlet)
    m.fs.unit_to_product = Arc(source=m.fs.unit.outlet, destination=m.fs.product.inlet)

    m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"].fix(0.965)  # kg/s
    m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "TDS"].fix(0.035)  # kg/s
    m.fs.feed.properties[0].pressure.fix(1 * pyunits.bar)  # Pa
    m.fs.feed.properties[0].temperature.fix(298.15)  # K

    m.fs.unit.control_volume.deltaP.fix(10 * pyunits.bar)  # Pa

    m.fs.properties.set_default_scaling("flow_mass_phase_comp", 1, index=("Liq", "H2O"))
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e2, index=("Liq", "TDS")
    )

    TransformationFactory("network.expand_arcs").apply_to(m)

    iscale.calculate_scaling_factors(m)

    return m


@pytest.mark.unit
def test_head_loss():
    m = build_head_loss()

    assert isinstance(m.fs.unit.control_volume, ControlVolume0DBlock)
    assert hasattr(m.fs.unit, "inlet")
    assert hasattr(m.fs.unit, "outlet")
    assert hasattr(m.fs.unit, "deltaP")  # this is just a reference
    assert hasattr(m.fs.unit.control_volume, "deltaP")

    m.fs.feed.initialize()
    propagate_state(m.fs.feed_to_unit)
    initialization_tester(m, unit=m.fs.unit)
    propagate_state(m.fs.unit_to_product)
    m.fs.product.initialize()

    assert degrees_of_freedom(m) == 0

    results = solver.solve(m)
    assert_optimal_termination(results)
