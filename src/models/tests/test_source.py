import pytest

from pyomo.environ import (
    ConcreteModel,
    TransformationFactory,
    assert_optimal_termination,
    units as pyunits,
    value,
)
from pyomo.network import Arc

from idaes.core import FlowsheetBlock, UnitModelCostingBlock
import idaes.core.util.scaling as iscale
from idaes.core.util.testing import initialization_tester
from idaes.core.util.initialization import propagate_state
from idaes.models.unit_models import Feed, Product
from idaes.core.util.model_statistics import degrees_of_freedom

from watertap.property_models.seawater_prop_pack import SeawaterParameterBlock
from watertap.core.solvers import get_solver
from watertap.costing import WaterTAPCosting

from models.source import Source

solver = get_solver()


def build_source():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = SeawaterParameterBlock()
    m.fs.costing = WaterTAPCosting()
    m.fs.costing.base_currency = pyunits.USD_2021

    m.fs.unit = Source(property_package=m.fs.properties)
    m.fs.unit.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
    m.fs.product = Product(property_package=m.fs.properties)

    m.fs.unit_to_product = Arc(source=m.fs.unit.outlet, destination=m.fs.product.inlet)

    m.fs.unit.properties[0].flow_mass_phase_comp["Liq", "H2O"].fix(0.965)  # kg/s
    m.fs.unit.properties[0].flow_mass_phase_comp["Liq", "TDS"].fix(0.035)  # kg/s
    m.fs.unit.properties[0].pressure.fix(1 * pyunits.bar)  # Pa
    m.fs.unit.properties[0].temperature.fix(298.15)  # K

    m.fs.properties.set_default_scaling("flow_mass_phase_comp", 1, index=("Liq", "H2O"))
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e2, index=("Liq", "TDS")
    )

    TransformationFactory("network.expand_arcs").apply_to(m)

    iscale.calculate_scaling_factors(m)

    m.fs.costing.cost_process()
    m.fs.costing.add_LCOW(m.fs.product.properties[0].flow_vol_phase["Liq"])

    return m


@pytest.mark.component
def test_head_loss():
    m = build_source()

    assert not hasattr(m.fs.unit, "inlet")  # only has an outlet
    assert hasattr(m.fs.unit, "outlet")

    initialization_tester(m, unit=m.fs.unit)
    propagate_state(m.fs.unit_to_product)
    m.fs.product.initialize()

    assert degrees_of_freedom(m) == 0

    results = solver.solve(m)
    assert_optimal_termination(results)

    assert pytest.approx(value(m.fs.costing.LCOW), rel=1e-5) == value(
        m.fs.costing.source.unit_cost
    )
