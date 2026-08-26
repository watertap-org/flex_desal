from watertap.core.solvers import get_solver
from parameter_sweep import LinearSample
from wrd.components.ro_stage import initialize_ro_stage
from wrd.components.ro_train import *
from models import HeadLoss
from pyomo.environ import (
    assert_optimal_termination,
    ConcreteModel,
    TransformationFactory,
    units as pyunits,
    value,
    Objective,
    minimize,
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


# Debugging Tool
def print_arcs_to_port(model, port, port_name=""):
    """Print all arcs connected to a given port (as source or destination)"""
    from pyomo.network import Arc

    print(f"\n{'='*60}")
    print(f"Arcs connected to {port_name or port}:")
    print(f"{'='*60}")

    found = False
    for component in model.component_data_objects(Arc, descend_into=True):
        name = component.getname(fully_qualified=True)
        if component.destination is port:
            print(f"  ← {name}: {component.source} → {port_name or port}")
            found = True
        elif component.source is port:
            print(f"  → {name}: {port_name or port} → {component.destination}")
            found = True

    if not found:
        print(f"  (No arcs found)")
    print()


# Debugging Tool
def report_ports_with_multiple_arcs(model):
    """Report ports with more than one incoming or outgoing Arc."""
    from collections import defaultdict
    from pyomo.network import Arc

    incoming = defaultdict(list)
    outgoing = defaultdict(list)

    for arc in model.component_data_objects(Arc, descend_into=True):
        arc_name = arc.getname(fully_qualified=True)
        outgoing[id(arc.source)].append((arc.source, arc_name))
        incoming[id(arc.destination)].append((arc.destination, arc_name))

    found_issue = False
    print(f"\n{'='*60}")
    print("Port arc cardinality check")
    print(f"{'='*60}")

    for entries in incoming.values():
        if len(entries) > 1:
            port = entries[0][0]
            arcs = [name for _, name in entries]
            print(f"IN  > 1 : {port} <- {arcs}")
            found_issue = True

    for entries in outgoing.values():
        if len(entries) > 1:
            port = entries[0][0]
            arcs = [name for _, name in entries]
            print(f"OUT > 1 : {port} -> {arcs}")
            found_issue = True

    if not found_issue:
        print("No ports with multiple incoming/outgoing arcs found.")
    print()
    return found_issue


# Build flowsheet function
def build_flowsheet(op_limits=None, scenario=None):
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = NaClParameterBlock()

    m.fs.feed = Feed(property_package=m.fs.properties)
    touch_flow_and_conc(m.fs.feed)

    m.fs.ro_train = FlowsheetBlock(dynamic=False)
    build_ro_train(
        m.fs.ro_train,
        num_stages=3,
        prop_package=m.fs.properties,
        file="wrd_inputs_8_19_21.yaml",
    )

    m.fs.product = Product(property_package=m.fs.properties)
    touch_flow_and_conc(m.fs.product)
    m.fs.brine = Product(property_package=m.fs.properties)
    touch_flow_and_conc(m.fs.brine)

    # Arcs to connect the unit models
    m.fs.feed_to_train = Arc(
        source=m.fs.feed.outlet,
        destination=m.fs.ro_train.feed.inlet,
    )
    m.fs.train_to_product = Arc(
        source=m.fs.ro_train.product.outlet,
        destination=m.fs.product.inlet,
    )
    m.fs.train_to_brine = Arc(
        source=m.fs.ro_train.disposal.outlet,
        destination=m.fs.brine.inlet,
    )

    # Add the pressure drop between PRO and TSRO
    m.fs.ro_train.tsro_header = HeadLoss(property_package=m.fs.properties)
    touch_flow_and_conc(m.fs.ro_train.tsro_header)
    m.fs.ro_train.tsro_header.control_volume.deltaP[0].fix(
        get_config_value(m.fs.ro_train.config_data, "tsro_header_loss", "headers")
    )

    # Undo connection between the second and third stage to add in the head loss unit.
    # Deactivate (do not delete) to keep port arc bookkeeping consistent.
    if hasattr(m.fs.ro_train, "stage_2_to_stage_3"):
        m.fs.ro_train.stage_2_to_stage_3.deactivate()
    if hasattr(m.fs.ro_train, "stage_2_to_stage_3_expanded"):
        m.fs.ro_train.stage_2_to_stage_3_expanded.deactivate()
    m.fs.ro_train.stage_2_to_tsro_header = Arc(
        source=m.fs.ro_train.stage[2].disposal.outlet,
        destination=m.fs.ro_train.tsro_header.inlet,
    )
    m.fs.ro_train.tsro_header_to_stage_3 = Arc(
        source=m.fs.ro_train.tsro_header.outlet,
        destination=m.fs.ro_train.stage[3].feed.inlet,
    )

    # report_ports_with_multiple_arcs(m)

    TransformationFactory("network.expand_arcs").apply_to(m)

    # Below are a few different options to address the issue of negative 3rd stage power
    # Add a lower bound on the stage 3 pump work
    m.fs.ro_train.stage[3].pump.unit.work_mechanical[0].setlb(0)
    # Increase the pump efficiency so stage 2 isn't so heavily favored
    m.fs.ro_train.stage[3].pump.unit.efficiency_pump.fix(0.7)
    # Provide a lower bound on the third stage recovery
    # m.fs.ro_train.stage[3].ro.unit.recovery_vol_phase[0,'Liq'].setlb(.5)

    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e-1, index=("Liq", "H2O")  # changed from 1
    )
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e2, index=("Liq", "NaCl")
    )

    set_ro_train_scaling(m.fs.ro_train)
    calculate_scaling_factors(m)

    # Set feed and operational conditions
    m.fs.feed.properties.calculate_state(
        var_args={
            ("flow_vol_phase", ("Liq")): 2500 * pyunits.gallons / pyunits.minute,
            ("conc_mass_phase_comp", ("Liq", "NaCl")): 0.5 * pyunits.g / pyunits.L,
            ("pressure", None): 35.4 * pyunits.psi,
            ("temperature", None): 298.15 * pyunits.K,
        },
        hold_state=True,
    )
    set_ro_train_op_conditions(m.fs.ro_train)
    # Set limits on each stage recovery and flowrate
    for i in m.fs.ro_train.stages:
        # m.fs.ro_train.stage[i].feed.properties[0].flow_vol_phase["Liq"].setlb(
        #     op_limits[f"Stage {i}"]["Qin_min"]
        # )
        # m.fs.ro_train.stage[i].feed.properties[0].flow_vol_phase["Liq"].setub(
        #     op_limits[f"Stage {i}"]["Qin_max"]
        # )
        m.fs.ro_train.stage[i].disposal.properties[0].flow_vol_phase["Liq"].setlb(
            op_limits[f"Stage {i}"]["Qout_min"]
        )
        # m.fs.ro_train.stage[i].ro.unit.recovery_vol_phase[0, "Liq"].setlb(op_limits[f"Stage {i}"]["RR_min"])
        # m.fs.ro_train.stage[i].ro.unit.recovery_vol_phase[0, "Liq"].setub(op_limits[f"Stage {i}"]["RR_max"])
    print(f"Degrees of freedom: {degrees_of_freedom(m)}")  # Should be zero.
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
        sweep_params["Recovery"] = LinearSample(
            m.fs.ro_train.recovery_vol,
            var_lims["RR_lb"],
            var_lims["RR_ub"],
            num_samples,
        )

        sweep_params["Feed Flow"] = LinearSample(
            m.fs.feed.properties[0].flow_vol_phase["Liq"],
            var_lims["Qin_lb"],
            var_lims["Qin_ub"],
            num_samples,
        )
    else:
        raise KeyError(
            f"Scenario {scenario} not recognized. Please choose from: 'default'"
        )
    return sweep_params


def initialize_model(m):
    # Expanding arcs here instead

    # Initialize system
    # Changing the default pressure to have a more medium recovery
    m.fs.ro_train.stage[1].pump.unit.control_volume.properties_out[0].pressure.fix(
        134 * pyunits.psi
    )
    m.fs.ro_train.stage[2].pump.unit.control_volume.properties_out[0].pressure.fix(
        150 * pyunits.psi
    )
    m.fs.ro_train.stage[3].pump.unit.control_volume.properties_out[0].pressure.fix(
        135 * pyunits.psi
    )  # This value could be increased if another water perm value is used for S3
    assert degrees_of_freedom(m) == 0
    m.fs.feed.initialize()
    propagate_state(m.fs.feed_to_train)
    # Initialize train, including the tsro header
    blk = m.fs.ro_train
    a = blk.find_component("feed_to_stage_1")
    propagate_state(a)
    initialize_ro_stage(blk.stage[1])
    a = blk.find_component("stage_1_to_product")
    propagate_state(a)
    a = blk.find_component("stage_1_to_stage_2")
    propagate_state(a)
    initialize_ro_stage(blk.stage[2])
    a = blk.find_component("stage_2_to_tsro_header")
    propagate_state(a)
    a = blk.find_component("stage_2_to_product")
    propagate_state(a)
    blk.tsro_header.initialize()
    a = blk.find_component("tsro_header_to_stage_3")
    propagate_state(a)
    initialize_ro_stage(blk.stage[3])
    a = blk.find_component("stage_3_to_product")
    propagate_state(a)
    a = blk.find_component("stage_3_to_brine")
    propagate_state(a)
    # a = blk.find_component("stage_3_to_product")
    # propagate_state(a)
    blk.mixer.initialize()
    propagate_state(blk.mixer_to_product)
    blk.product.initialize()

    propagate_state(m.fs.train_to_product)
    m.fs.product.initialize()
    propagate_state(m.fs.train_to_brine)
    m.fs.brine.initialize()

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
    # Pump outlet pressures unfixed
    m.fs.ro_train.stage[1].pump.unit.control_volume.properties_out[0].pressure.unfix()
    m.fs.ro_train.stage[2].pump.unit.control_volume.properties_out[0].pressure.unfix()
    m.fs.ro_train.stage[3].pump.unit.control_volume.properties_out[0].pressure.unfix()
    # Reduce DOF by assume S1 and S2 have same recovery, which is about typical from data, but not necessarily optimal...
    # m.fs.ro_train.eq_recovery = Constraint(
    #     expr= m.fs.ro_train.stage[1].ro.unit.recovery_vol_phase[0, "Liq"] == m.fs.ro_train.stage[2].ro.unit.recovery_vol_phase[0, "Liq"])

    m.fs.obj = Objective(
        expr=m.fs.ro_train.total_pump_power, sense=minimize
    )  # Dummy objective to trigger solve

    print(
        f"Degrees of freedom right before solve: {degrees_of_freedom(m)}"
    )  # Should be 0
    # assert degrees_of_freedom(m) == 0
    # --solve---
    solver = get_solver()
    try:
        results = solver.solve(m)
        assert_optimal_termination(results)
    except Exception as e:
        print("-----FAILED TO SOLVE-----")
        print(f"Exception: {e}")
        print(f"RR = {m.fs.ro_train.stage[1].ro.unit.recovery_vol_phase[0, 'Liq']()}%")
        print(f"Qin = {m.fs.feed.properties[0].flow_vol_phase['Liq']()}")
        # dt = DiagnosticsToolbox(m)
        # dt.report_numerical_issues()
        # dt.display_constraints_with_large_residuals()
    return results


def build_outputs(m):
    def _pump_unit(stage_blk):
        # Some workflows expose the pump model directly, others nest it under stage.pump.unit
        return (
            stage_blk.pump.unit if hasattr(stage_blk.pump, "unit") else stage_blk.pump
        )

    outputs = {}
    # RR and Flows
    outputs["Stage1 RR"] = m.fs.ro_train.stage[1].ro.unit.recovery_vol_phase[0, "Liq"]
    outputs["Stage2 RR"] = m.fs.ro_train.stage[2].ro.unit.recovery_vol_phase[0, "Liq"]
    outputs["Stage3 RR"] = m.fs.ro_train.stage[3].ro.unit.recovery_vol_phase[0, "Liq"]
    outputs["Stage1 Qin (m3/s)"] = (
        m.fs.ro_train.stage[1].feed.properties[0].flow_vol_phase["Liq"]
    )
    outputs["Stage2 Qin (m3/s)"] = (
        m.fs.ro_train.stage[2].feed.properties[0].flow_vol_phase["Liq"]
    )
    outputs["Stage3 Qin (m3/s)"] = (
        m.fs.ro_train.stage[3].feed.properties[0].flow_vol_phase["Liq"]
    )
    outputs["Stage1 Brine (m3/s)"] = (
        m.fs.ro_train.stage[1].disposal.properties[0].flow_vol_phase["Liq"]
    )
    outputs["Stage2 Brine (m3/s)"] = (
        m.fs.ro_train.stage[2].disposal.properties[0].flow_vol_phase["Liq"]
    )
    outputs["Stage3 Brine (m3/s)"] = (
        m.fs.ro_train.stage[3].disposal.properties[0].flow_vol_phase["Liq"]
    )

    pump1 = _pump_unit(m.fs.ro_train.stage[1])
    pump2 = _pump_unit(m.fs.ro_train.stage[2])
    pump3 = _pump_unit(m.fs.ro_train.stage[3])

    # Pressures
    outputs["Stage1 Pin (Pa)"] = pump1.control_volume.properties_in[0].pressure
    outputs["Stage2 Pin (Pa)"] = pump2.control_volume.properties_in[0].pressure
    outputs["Stage3 Pin (Pa)"] = pump3.control_volume.properties_in[0].pressure
    outputs["Stage1 Pout (Pa)"] = pump1.control_volume.properties_out[0].pressure
    outputs["Stage2 Pout (Pa)"] = pump2.control_volume.properties_out[0].pressure
    outputs["Stage3 Pout (Pa)"] = pump3.control_volume.properties_out[0].pressure

    # Powers
    outputs["Total Power (W)"] = m.fs.ro_train.total_pump_power
    outputs["Stage1 Power (W)"] = pump1.work_mechanical[0]
    outputs["Stage2 Power (W)"] = pump2.work_mechanical[0]
    outputs["Stage3 Power (W)"] = pump3.work_mechanical[0]

    # Other
    outputs["Stage1 First Elem. Prod. (m3/s)"] = m.fs.ro_train.stage[
        1
    ].ro.unit.first_elem_prod
    outputs["Stage2 First Elem. Prod. (m3/s)"] = m.fs.ro_train.stage[
        2
    ].ro.unit.first_elem_prod
    outputs["Stage3 First Elem. Prod. (m3/s)"] = m.fs.ro_train.stage[
        3
    ].ro.unit.first_elem_prod
    outputs["Stage1 Perm. Conc. (kg/m3)"] = (
        m.fs.ro_train.stage[1]
        .ro.unit.mixed_permeate[0]
        .conc_mass_phase_comp["Liq", "NaCl"]
    )
    outputs["Stage2 Perm. Conc. (kg/m3)"] = (
        m.fs.ro_train.stage[2]
        .ro.unit.mixed_permeate[0]
        .conc_mass_phase_comp["Liq", "NaCl"]
    )
    outputs["Stage3 Perm. Conc. (kg/m3)"] = (
        m.fs.ro_train.stage[3]
        .ro.unit.mixed_permeate[0]
        .conc_mass_phase_comp["Liq", "NaCl"]
    )
    outputs["Stage1 Pump Speed (-)"] = pump1.design_speed_fraction
    outputs["Stage2 Pump Speed (-)"] = pump2.design_speed_fraction
    # outputs["Stage3 Pump Speed (-)"] = pump3.design_speed_fraction # Pump 3 has no speed fraction

    return outputs


# Debugging
if __name__ == "__main__":
    op_limits = {
        "Stage 1": {
            "Qout_min": 3
            * 72
            / 3600,  # This limits will bound the flowrate for a given recovery.  Equal to 3 m3/hr * 72 Pressure Vessels per train
            # "Qin_max": 635 / 3600, # Based on pump limitation
        },
        "Stage 2": {
            "Qout_min": 3 * 30 / 3600,
            # "Qin_max": 251 / 3600, # This came from stage 1 min recovery of 55%
        },
        "Stage 3": {
            "Qout_min": 3 * 15 / 3600,
            # "Qin_max": 126 / 3600,
        },
    }

    m = build_flowsheet(op_limits=op_limits, scenario=None)
    initialize_model(m)

    # Dummy version of fixing value
    print(f"Degrees of freedom before fixing: {degrees_of_freedom(m)}")

    # Param Sweep will fix these two variables. This is the maximum of both flowrate and recovery
    m.fs.ro_train.recovery_vol.fix(0.925)
    m.fs.feed.properties[0].flow_vol_phase["Liq"].fix(0.176)

    print(f"Degrees of freedom after fixing: {degrees_of_freedom(m)}")
    results = optimize(m)

    print(
        f"Stage 3 pump outlet conc (kg/m3): {m.fs.ro_train.stage[3].pump.unit.control_volume.properties_out[0].conc_mass_phase_comp['Liq', 'NaCl']():.1f}"
    )
    print(
        f"Stage 2 pump outlet pressure (psi): {value(pyunits.convert(m.fs.ro_train.stage[2].pump.unit.control_volume.properties_out[0].pressure, to_units=pyunits.psi)):.1f}"
    )
    print(
        f"Stage 3 pump inlet pressure (psi): {value(pyunits.convert(m.fs.ro_train.stage[3].pump.unit.control_volume.properties_in[0].pressure, to_units=pyunits.psi)):.1f}"
    )
    print(
        f"Stage 3 pump outlet pressure (psi): {value(pyunits.convert(m.fs.ro_train.stage[3].pump.unit.control_volume.properties_out[0].pressure, to_units=pyunits.psi)):.1f}"
    )

    print(
        f"Stage 3 feed flow vol (m3/s): {m.fs.ro_train.stage[3].feed.properties[0].flow_vol_phase['Liq']():.3f}"
    )
    print(
        f"Stage 1 RR: {m.fs.ro_train.stage[1].ro.unit.recovery_vol_phase[0, 'Liq']():.3f}"
    )
    print(
        f"Stage 2 RR: {m.fs.ro_train.stage[2].ro.unit.recovery_vol_phase[0, 'Liq']():.3f}"
    )
    print(
        f"Stage 3 RR: {m.fs.ro_train.stage[3].ro.unit.recovery_vol_phase[0, 'Liq']():.3f}"
    )

    assert_optimal_termination(results)
