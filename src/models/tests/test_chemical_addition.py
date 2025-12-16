import pytest

from pyomo.environ import (
    ConcreteModel,
    TransformationFactory,
    assert_optimal_termination,
    units as pyunits,
)
from pyomo.network import Arc

from idaes.core import FlowsheetBlock, UnitModelCostingBlock
import idaes.core.util.scaling as iscale
from idaes.core.util.testing import initialization_tester
from idaes.core.util.initialization import propagate_state
from idaes.models.unit_models import Feed, Product
from idaes.core.util.model_statistics import degrees_of_freedom

from watertap.property_models.seawater_prop_pack import SeawaterParameterBlock
from watertap.core.control_volume_isothermal import ControlVolume0DBlock
from watertap.core.solvers import get_solver
from watertap.costing import WaterTAPCosting

from models.chemical_addition import ChemicalAddition

solver = get_solver()


def build_chem_addition_model(chemical="ammonia"):
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties = SeawaterParameterBlock()

    m.fs.unit = ChemicalAddition(property_package=m.fs.properties, chemical=chemical)
    m.fs.unit.dose.fix(0.1)  # g/L

    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e-1, index=("Liq", "H2O")
    )
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e2, index=("Liq", "TDS")
    )

    m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(0.965)  # kg/s
    m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "TDS"].fix(0.035)  # kg/s
    m.fs.unit.inlet.temperature[0].fix(293)
    m.fs.unit.inlet.pressure[0].fix(101325)

    iscale.calculate_scaling_factors(m)

    return m


@pytest.mark.component
def test_chemical_addition():

    for chem in ["ammonia", "ferric_chloride", "sulfuric_acid"]:
        m = build_chem_addition_model(chemical=chem)

        # Check unit configuration
        assert len(m.fs.unit.config) == 6

        # Check that variables are constructed
        assert hasattr(m.fs.unit, "dose")
        assert hasattr(m.fs.unit, "chemical_soln_flow_vol")
        assert hasattr(m.fs.unit, "pumping_power")

        initialization_tester(m)

        # Solve model
        results = solver.solve(m)
        # Check for optimal solution
        assert_optimal_termination(results)


@pytest.mark.component
def test_chemical_addition_costing():
    for chem in ["ammonia", "ferric_chloride", "sulfuric_acid"]:
        m = build_chem_addition_model(chemical=chem)

        m.fs.costing = WaterTAPCosting()
        m.fs.unit.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)

        m.fs.costing.cost_process()
        initialization_tester(m)

        # Solve model
        results = solver.solve(m)
        # Check for optimal solution
        assert_optimal_termination(results)
