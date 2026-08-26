from watertap.core.solvers import get_solver
from parameter_sweep import LinearSample
from wrd.components.detailed_pump import *
from pyomo.environ import (
    assert_optimal_termination,
    ConcreteModel,
    TransformationFactory,
    units as pyunits,
)
from pyomo.network import Arc
from srp.utils import touch_flow_and_conc
from idaes.models.unit_models import (
    Feed,
    Product,
)
from idaes.core import FlowsheetBlock
from watertap.property_models.NaCl_T_dep_prop_pack import NaClParameterBlock
from idaes.core.util.scaling import calculate_scaling_factors
from idaes.core.util.initialization import propagate_state
from idaes.core.util.model_statistics import degrees_of_freedom
from wrd.utilities import *
from idaes.core.util.model_diagnostics import DiagnosticsToolbox

__all__ = [
    "build_flowsheet",
    "build_sweep_params",
    "initialize_model",
    "optimize",
    "build_outputs",
]


# Build flowsheet function
def build_flowsheet(op_limits=None, scenario=None):
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = NaClParameterBlock()

    m.fs.feed = Feed(property_package=m.fs.properties)
    touch_flow_and_conc(m.fs.feed)
    m.fs.pump = FlowsheetBlock(dynamic=False)
    build_pump(
        m.fs.pump,
        stage_num=1,
        file="wrd_inputs_8_19_21.yaml",
        prop_package=m.fs.properties,
        uf=True,
    )

    m.fs.product = Product(property_package=m.fs.properties)
    touch_flow_and_conc(m.fs.product)

    # Arcs to connect the unit models
    m.fs.feed_to_pump = Arc(
        source=m.fs.feed.outlet,
        destination=m.fs.pump.feed.inlet,
    )
    m.fs.pump_to_product = Arc(
        source=m.fs.pump.product.outlet,
        destination=m.fs.product.inlet,
    )

    TransformationFactory("network.expand_arcs").apply_to(m)

    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e-1, index=("Liq", "H2O")
    )
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e2, index=("Liq", "NaCl")
    )

    add_pump_scaling(m.fs.pump)
    calculate_scaling_factors(m)

    # Set feed and operational conditions
    m.fs.feed.properties[0].pressure.setlb(0)
    m.fs.pump.unit.control_volume.properties_in[0].pressure.setlb(0)

    m.fs.feed.properties.calculate_state(
        var_args={
            ("flow_vol_phase", ("Liq")): 3000 * pyunits.gallons / pyunits.minute,
            ("conc_mass_phase_comp", ("Liq", "NaCl")): 0.5 * pyunits.g / pyunits.L,
            ("pressure", None): 0.1 * pyunits.psi,
            ("temperature", None): 298.15 * pyunits.K,
        },
        hold_state=True,
    )

    set_pump_op_conditions(m.fs.pump)
    print(degrees_of_freedom(m))  # Should be zero

    return m


# Build sweep parameters function
def build_sweep_params(
    m,
    num_samples=10,
    var_lims=None,
    scenario=None,
):
    sweep_params = {}
    if scenario == "default":
        sweep_params["Feed Flow"] = LinearSample(
            m.fs.feed.properties[0].flow_vol_phase["Liq"],
            var_lims["Feed Flow"]["Qin_min"],
            var_lims["Feed Flow"]["Qin_max"],
            num_samples,
        )
    else:
        raise KeyError(
            f"Scenario {scenario} not recognized. Please choose from: 'default'"
        )
    return sweep_params


def initialize_model(m):
    if m.fs.pump.uf:
        # Allow negative suction pressure for UF configuration
        m.fs.feed.properties[0].pressure.setlb(0)
        m.fs.pump.feed.properties[0].pressure.setlb(0)

        # Change the bounds for the pump inlet pressure
        m.fs.pump.unit.control_volume.properties_in[0].pressure.setlb(0)
    # Initialize system
    m.fs.feed.initialize()
    propagate_state(m.fs.feed_to_pump)
    initialize_pump(m.fs.pump)
    propagate_state(m.fs.pump_to_product)
    m.fs.product.initialize()

    assert degrees_of_freedom(m) == 0
    # --solve---
    solver = get_solver()
    results = solver.solve(m)
    assert_optimal_termination(results)
    return results


# Optimizaiton function
def optimize(m, solver=None, check_termination=True):
    print(f"Degrees of freedom start of optimize: {degrees_of_freedom(m)}")
    # Fix the mass fraction and unfix the flow mass so that the flow volume can be fixed by param sweep
    m.fs.feed.properties[0].mass_frac_phase_comp["Liq", "NaCl"].fix()
    m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"].unfix()
    m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "NaCl"].unfix()

    print(
        f"Degrees of freedom right before solve: {degrees_of_freedom(m)}"
    )  # Should be 0
    assert degrees_of_freedom(m) == 0
    # --solve---
    solver = get_solver()
    results = solver.solve(m, tee=True)
    # assert_optimal_termination(results)
    return results


def build_outputs(m):
    outputs = {}
    outputs["Stage1 Power"] = m.fs.pump.unit.work_mechanical[0]
    return outputs


# Debugging
if __name__ == "__main__":
    m = build_flowsheet(scenario=None)
    initialize_model(m)
    # Dummy version of fixing value
    print(f"Degrees of freedom before fixing: {degrees_of_freedom(m)}")  # Should be 0
    # Param Sweep will fix these two variables
    m.fs.feed.properties[0].flow_vol_phase["Liq"].fix(0.158)

    print(f"Degrees of freedom after fixing: {degrees_of_freedom(m)}")
    results = optimize(m)

    dt = DiagnosticsToolbox(m)
    assert_optimal_termination(results)
    print(m.fs.pump.unit.work_mechanical[0]())
