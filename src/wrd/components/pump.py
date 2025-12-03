import os
import yaml
from pyomo.environ import (
    ConcreteModel,
    Objective,
    Var,
    Param,
    Constraint,
    NonNegativeReals,
    TransformationFactory,
    assert_optimal_termination,
    value,
    units as pyunits,
)
from pyomo.network import Arc

from idaes.core import FlowsheetBlock, UnitModelCostingBlock
from idaes.core.util.initialization import propagate_state
from idaes.models.unit_models import StateJunction
from idaes.core.util.model_statistics import degrees_of_freedom
import idaes.core.util.scaling as iscale
from idaes.core.util.scaling import (
    calculate_scaling_factors,
    set_scaling_factor,
)

from watertap.property_models.NaCl_prop_pack import NaClParameterBlock
from watertap.unit_models.pressure_changer import Pump
from watertap.core.solvers import get_solver


def load_config(config):
    with open(config, "r") as file:
        return yaml.safe_load(file)


def get_config_value(
    config,
    key,
    section,
    subsection=None,
):
    """
    Get a value from the configuration file.
    """

    if section in config:
        if subsection:
            if subsection in config[section]:
                if key in config[section][subsection]:
                    if (
                        isinstance(config[section][subsection][key], dict)
                        and "value" in config[section][subsection][key]
                        and "units" in config[section][subsection][key]
                    ):
                        return config[section][subsection][key]["value"] * getattr(
                            pyunits, config[section][subsection][key]["units"]
                        )
                    return config[section][subsection][key]
                else:
                    raise KeyError(
                        f"Key '{key}' not found in subsection '{subsection}' of section '{section}' of the configuration."
                    )
            else:
                raise KeyError(
                    f"Section '{section}' or subsection '{subsection}' not found in the configuration."
                )
        else:
            if key in config[section]:
                if (
                    isinstance(config[section][key], dict)
                    and "value" in config[section][key]
                    and "units" in config[section][key]
                ):
                    return config[section][key]["value"] * getattr(
                        pyunits, config[section][key]["units"]
                    )
                return config[section][key]
            else:
                raise KeyError(
                    f"Key '{key}' not found in section '{section}' of the configuration."
                )
    else:
        raise KeyError(f"Section '{section}' not found in the configuration.")


def build_system(**kwargs):
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.ro_properties = NaClParameterBlock()
    m.fs.pump_system = FlowsheetBlock(dynamic=False)
    build_wrd_pump(m.fs.pump_system, prop_package=m.fs.ro_properties, **kwargs)
    return m


def build_wrd_pump(blk, stage_num=1, prop_package=None):
    m = blk.model()

    if prop_package is None:
        prop_package = m.fs.ro_properties

    blk.feed_in = StateJunction(property_package=prop_package)
    blk.feed_out = StateJunction(property_package=prop_package)

    # Get the absolute path of the current script
    current_script_path = os.path.abspath(__file__)
    # Get the directory containing the current script
    current_directory = os.path.dirname(current_script_path)
    # Get the parent directory of the current directory (one folder prior)
    parent_directory = os.path.dirname(current_directory)

    config = os.path.join(parent_directory, "meta_data", "wrd_ro_inputs.yaml")

    blk.config_data = load_config(config)
    blk.pump = Pump(property_package=prop_package)

    #Create Variables for simple "surrogate"
     
    blk.pump.efficiency_eq_constant = Param(
        initialize=0.389,
        mutable=True,
        units=pyunits.dimensionless,
        doc="Constant term of Efficiency equation",
    )

    blk.pump.efficiency_eq_linear = Param(
        initialize= -0.535,
        mutable=True,
        units= (pyunits.m**3 / pyunits.s)**-1,
        doc="Linear term of Efficiency equation",
    )

    blk.pump.efficiency_eq_squared = Param(
            initialize= 41.373,
            mutable=True,
            units= (pyunits.m**3 / pyunits.s)**-2,
            doc="Squared term of Efficiency equation",
        )
    
    blk.pump.efficiency_eq_cubed = Param(
            initialize= -138.82,
            mutable=True,
            units=(pyunits.m**3 / pyunits.s)**-3,
            doc="Cubed term of Efficiency equation",
        )
    flow = blk.feed_in.properties[0].flow_vol_phase["Liq"]
    blk.pump.efficiency_surr_eq = Constraint(
        expr = blk.pump.efficiency_pump[0] ==
        blk.pump.efficiency_eq_cubed * flow**3 + blk.pump.efficiency_eq_squared * flow**2
        + blk.pump.efficiency_eq_linear * flow + blk.pump.efficiency_eq_constant,
        doc = "Efficiency surrogate equation"
    )

    # Add Arcs
    blk.feed_in_to_pump = Arc(source=blk.feed_in.outlet, destination=blk.pump.inlet)
    blk.pump_to_feed_out = Arc(source=blk.pump.outlet, destination=blk.feed_out.inlet)
    TransformationFactory("network.expand_arcs").apply_to(blk)


def set_pump_op_conditions(blk, stage_num=1):
    # Configure with input values
    # blk.pump.efficiency_pump.fix(
    #     get_config_value(
    #         blk.config_data, "pump_efficiency", "pumps", f"pump_{stage_num}"
    #     )
    # )
    blk.pump.control_volume.properties_out[0].pressure.fix(
        get_config_value(
            blk.config_data, "pump_outlet_pressure", "pumps", f"pump_{stage_num}"
        )
    )
    # This will be written over by the surrogate model for the pumps


def set_inlet_conditions(blk, Qin=0.154, Cin=0.542, P_in=1):
    """
    Set the operation conditions for the Pump
    """
    Qin = (Qin) * pyunits.m**3 / pyunits.s  # Feed flow rate in m3/s
    Cin = Cin * pyunits.g / pyunits.L  # Feed concentration in g/L
    rho = 1000 * pyunits.kg / pyunits.m**3  # Approximate density of water
    feed_mass_flow_water = Qin * rho
    feed_mass_flow_salt = Cin * Qin

    blk.feed_in.properties[0].flow_mass_phase_comp["Liq", "H2O"].fix(
        feed_mass_flow_water
    )
    blk.feed_in.properties[0].flow_mass_phase_comp["Liq", "NaCl"].fix(
        feed_mass_flow_salt
    )
    blk.feed_in.properties[0].temperature.fix(298.15 * pyunits.K)  # 25 C
    blk.feed_in.properties[0].pressure.fix(P_in * pyunits.bar)
    blk.feed_in.properties[0].flow_vol  # Touching

    # Scaling properties
    m = blk.model()
    m.fs.ro_properties.set_default_scaling(
        "flow_mass_phase_comp", 1e-1, index=("Liq", "H2O")
    )
    m.fs.ro_properties.set_default_scaling(
        "flow_mass_phase_comp", 1e2, index=("Liq", "NaCl")
    )


def add_pump_scaling(blk):
    # Properties
    set_scaling_factor(blk.pump.work_mechanical[0], 1e-3)
    # Isn't there a needed scaling factor for electricity costs?


def initialize_pump(blk):
    blk.feed_in.initialize()
    propagate_state(blk.feed_in_to_pump)
    blk.pump.initialize()
    propagate_state(blk.pump_to_feed_out)
    blk.feed_out.initialize()


def report_pump(blk, w=30):
    title = "Pump Report"
    side = int(((3 * w) - len(title)) / 2) - 1
    header = "=" * side + f" {title} " + "=" * side
    print(f"\n{header}\n")
    print(f'{"Parameter":<{w}s}{"Value":<{w}s}{"Units":<{w}s}')
    print(f"{'-' * (3 * w)}")

    total_flow = blk.feed_in.properties[0].flow_vol
    deltaP = blk.pump.deltaP[0]
    work = blk.pump.work_mechanical[0]
    print(
        f'{f"Total Flow Rate (MGD)":<{w}s}{value(pyunits.convert(total_flow, to_units=pyunits.Mgallons /pyunits.day)):<{w}.3f}{"MGD"}'
    )
    print(f'{f"Total Flow Rate (m3/s)":<{w}s}{value(total_flow):<{w}.3e}{"m3/s"}')
    print(
        f'{f"Total Flow Rate (gpm)":<{w}s}{value(pyunits.convert(total_flow, to_units=pyunits.gallons / pyunits.minute)):<{w}.3f}{"gpm"}'
    )
    print(f'{f"Pressure Change (Pa)":<{w}s}{value(deltaP):<{w}.3e}{"Pa"}')
    print(
        f'{f"Pressure Change (bar)":<{w}s}{value(pyunits.convert(deltaP,to_units=pyunits.bar)):<{w}.3e}{"bar"}'
    )
    print(f'{f"Work Mech. (J)":<{w}s}{value(work):<{w}.3e}{"Joules"}')
    print(
        f'{f"Work Mech. (kW)":<{w}s}{value(pyunits.convert(work, to_units=pyunits.kW)):<{w}.3f}{"kW"}'
    )


def main():
    m = build_system()  # optional input of stage_num
    set_inlet_conditions(m.fs.pump_system, Qin=0.154, Cin=0.542, P_in=1)
    set_pump_op_conditions(m.fs.pump_system)
    add_pump_scaling(m.fs.pump_system)
    calculate_scaling_factors(m)
    initialize_pump(m.fs.pump_system)
    m.fs.obj = Objective(
        expr=m.fs.pump_system.feed_out.properties[0].flow_vol_phase["Liq"]
    )
    solver = get_solver()
    results = solver.solve(m)
    assert_optimal_termination(results)


if __name__ == "__main__":
    m = build_system()  # optional input of stage_num
    print(f"{degrees_of_freedom(m)} degrees of freedom after build")
    set_inlet_conditions(m.fs.pump_system, Qin=0.154, Cin=0.542, P_in=1)
    set_pump_op_conditions(m.fs.pump_system)
    print(f"{degrees_of_freedom(m)} degrees of freedom after setting op conditions")
    add_pump_scaling(m.fs.pump_system)
    calculate_scaling_factors(m)
    initialize_pump(m.fs.pump_system)
    m.fs.obj = Objective(
        expr=m.fs.pump_system.feed_out.properties[0].flow_vol_phase["Liq"]
    )  # There is no D.o.f to optimize with
    solver = get_solver()
    results = solver.solve(m)
    assert_optimal_termination(results)
    # print(f"{iscale.jacobian_cond(m.fs.pump_system):.2e}")
    report_pump(m.fs.pump_system)
