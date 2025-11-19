from pyomo.environ import (
    ConcreteModel,
    Param,
    Objective,
    Expression,
    check_optimal_termination,
    value,
    assert_optimal_termination,
    units as pyunits,
    value,
    TransformationFactory,
)


from idaes.core.util.initialization import propagate_state
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core import FlowsheetBlock
from idaes.models.unit_models import (
    MixingType,
    MomentumMixingType,
    Mixer,
    Separator,
    Product,
    Feed,
    StateJunction,
)

from pyomo.network import Arc
import idaes.core.util.scaling as iscale

from watertap.unit_models.reverse_osmosis_1D import (
    ReverseOsmosis1D,
    PressureChangeType,
    MassTransferCoefficient,
    ConcentrationPolarizationType,
)
from watertap.core.util.model_diagnostics.infeasible import *
from watertap.property_models.NaCl_T_dep_prop_pack import NaClParameterBlock
from watertap.core.solvers import get_solver

from idaes.core.util.scaling import (
    constraint_scaling_transform,
    calculate_scaling_factors,
    set_scaling_factor,
    list_badly_scaled_variables,
    extreme_jacobian_rows,
)

import yaml
import os


solver = get_solver()


def load_config(config):
    with open(config, "r") as file:
        return yaml.safe_load(file)


def get_config_value( # Will need to update this to match new yaml
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


def relax_bounds_for_low_salinity_waters(blk):
    blk.feed_side.cp_modulus.setub(5)
    for e in blk.feed_side.K:
        blk.feed_side.K[e].setub(0.01)
        blk.feed_side.K[e].setlb(1e-7)

    for e in blk.feed_side.cp_modulus:
        blk.feed_side.cp_modulus[e].setlb(1e-5)

    for e in blk.recovery_mass_phase_comp:
        if e[-1] == "NaCl":
            blk.recovery_mass_phase_comp[e].setlb(1e-9)
            blk.recovery_mass_phase_comp[e].setub(1e-1)

    for e in blk.flux_mass_phase_comp:
        if e[-1] == "NaCl":
            blk.flux_mass_phase_comp[e].setlb(1e-9)
            blk.flux_mass_phase_comp[e].setub(1e-1)

    for e in blk.recovery_mass_phase_comp:
        if e[-1] == "H2O":
            blk.recovery_mass_phase_comp[e].setlb(1e-4)
            blk.recovery_mass_phase_comp[e].setub(0.999)

    for e in blk.flux_mass_phase_comp:
        if e[-1] == "H2O":
            blk.flux_mass_phase_comp[e].setlb(1e-5)
            blk.flux_mass_phase_comp[e].setub(0.999)


def build_system(**kwargs):
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.ro_properties = NaClParameterBlock()

    m.fs.ro_system = FlowsheetBlock(dynamic=False)

    build_wrd_ro_system(m.fs.ro_system, prop_package=m.fs.ro_properties, **kwargs)

    return m


def build_wrd_ro_system(blk, prop_package=None):
    """
    Build reverse osmosis system for WRD
    """

    m = blk.model()
    if prop_package is None:
        prop_package = m.fs.ro_properties

    # Feed stream, permeate, and brine 
    blk.feed = StateJunction(property_package=prop_package)
    blk.permeate = Product(property_package=prop_package)
    blk.brine = Product(property_package=prop_package)

    blk.recovery = Var( # Creating variable for RR. This may no longer be needed here, but moved to main flowsheet as an overall recovery from the different stages
        initialize=0.5,
        bounds=(0, 1),
        units=pyunits.dimensionless,
        doc="Overall recovery",
    )

    @blk.Constraint(doc="Overall recovery constraint")
    def eq_recovery(b):
        return (
            b.recovery == b.permeate.properties[0].flow_vol_phase["Liq"] / b.feed.properties[0].flow_vol_phase["Liq"]
        )

    # WRD RO configurations input file. References to all values included in yml file
    # Get the absolute path of the current script
    current_script_path = os.path.abspath(__file__)
    # Get the directory containing the current script
    current_directory = os.path.dirname(current_script_path)
    # Get the parent directory of the current directory (one folder prior)
    parent_directory = os.path.dirname(current_directory)

    config = parent_directory + "/meta_data/wrd_ro_inputs.yaml"
    blk.config_data = load_config(config)

    """
    add_ro_units(train, prop_package=prop_package)
    Add RO units to a single RO train
    """

    blk.add_component(
        "ro", # May want to replace with a name
        ReverseOsmosis1D(
            property_package=prop_package,
            has_pressure_change=True,
            # pressure_change_type=PressureChangeType.calculated, # Why this setting?
            pressure_change_type=PressureChangeType.fixed_per_stage,
            mass_transfer_coefficient=MassTransferCoefficient.calculated,
            concentration_polarization_type=ConcentrationPolarizationType.calculated,
            transformation_scheme="BACKWARD",
            transformation_method="dae.finite_difference",
            module_type="spiral_wound",
            finite_elements=7,
            has_full_reporting=True,
        ),
    )

    """
    add_ro_connections(train)
    Add connections between the units in the RO system
    """
    # Connect permeate mixer to permeate product stream
    blk.RO_to_permeate = Arc(
        source=blk.feed.outlet, destination=blk.ro.inlet
    )

    blk.RO_to_permeate = Arc(
        source=blk.ro.outlet, destination=blk.permeate.inlet
    )
    blk.RO_to_brine = Arc(
        source=blk.ro.retentate, destination=blk.brine.inlet
    )
    TransformationFactory("network.expand_arcs").apply_to(blk)
    print("Degrees of freedom after adding units:", degrees_of_freedom(blk))


def set_inlet_conditions(blk, Qin=0.154, Cin=0.542):
    """
    Set the operation conditions for the RO system
    """
    Qin = (Qin) * pyunits.m**3 / pyunits.s  # Feed flow rate in m3/s
    Cin = Cin * pyunits.g / pyunits.L  # Feed concentration in g/L
    rho = 1000 * pyunits.kg / pyunits.m**3  # Approximate density of water
    feed_mass_flow_water = Qin * rho
    feed_mass_flow_salt = Cin * Qin

    blk.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"].fix(feed_mass_flow_water)
    blk.feed.properties[0].flow_mass_phase_comp["Liq", "NaCl"].fix(feed_mass_flow_salt)
    blk.feed.properties[0].temperature.fix(298.15 * pyunits.K)  # 25 C
    blk.feed.properties[0].pressure.fix(101325 * pyunits.Pa)  # 1 bar


def set_ro_system_op_conditions(blk):
    """
    Set the operation conditions for the RO system
    """
    # Set pump operating conditions
    for t in range(1, blk.number_trains + 1):
        train = blk.find_component(f"train_{t}")
        for s in range(1, (train.number_stages + 1)):
            pump = train.find_component(f"pump{s}")

            pump.control_volume.properties_out[0].pressure.fix(
                get_config_value(
                    blk.config_data, "pump_outlet_pressure", "pumps", f"pump_{s}"
                )
            )

            pump.efficiency_pump.fix(
                get_config_value(
                    blk.config_data, "pump_efficiency", "pumps", f"pump_{s}"
                )
            )

            # Set RO configuration for each stage
            ro_stage = train.find_component(f"ro_stage_{s}")

            ro_stage.A_comp.fix(
                get_config_value(
                    blk.config_data, "A_comp", "reverse_osmosis_1d", f"stage_{s}"
                )
            )
            ro_stage.B_comp.fix(
                get_config_value(
                    blk.config_data, "B_comp", "reverse_osmosis_1d", f"stage_{s}"
                )
            )

            ro_stage.feed_side.channel_height.fix(
                get_config_value(
                    blk.config_data,
                    "channel_height",
                    "reverse_osmosis_1d",
                    f"stage_{s}",
                )
            )
            ro_stage.feed_side.spacer_porosity.fix(
                get_config_value(
                    blk.config_data,
                    "spacer_porosity",
                    "reverse_osmosis_1d",
                    f"stage_{s}",
                )
            )

            ro_stage.feed_side.length.fix(
                get_config_value(
                    blk.config_data,
                    "number_of_elements_per_vessel",
                    "reverse_osmosis_1d",
                    f"stage_{s}",
                )
                * get_config_value(
                    blk.config_data,
                    "element_length",
                    "reverse_osmosis_1d",
                    f"stage_{s}",
                )
            )

            ro_stage.area.setub(1e6)
            ro_stage.width.setub(1e5)

            ro_stage.area.fix(
                get_config_value(
                    blk.config_data,
                    "element_membrane_area",
                    "reverse_osmosis_1d",
                    f"stage_{s}",
                )
                * get_config_value(
                    blk.config_data,
                    "number_of_vessels",
                    "reverse_osmosis_1d",
                    f"stage_{s}",
                )
                * get_config_value(
                    blk.config_data,
                    "number_of_elements_per_vessel",
                    "reverse_osmosis_1d",
                    f"stage_{s}",
                )
            )

            ro_stage.deltaP.fix(
                get_config_value(
                    blk.config_data, "pressure_drop", "reverse_osmosis_1d", f"stage_{s}"
                )
            )

            ro_stage.recovery_mass_phase_comp[0, "Liq", "H2O"].fix(
                get_config_value(
                    blk.config_data,
                    "water_recovery_mass_phase",
                    "reverse_osmosis_1d",
                    f"stage_{s}",
                )
            )


def add_ro_scaling(blk):
    """
    Add scaling to the units in the RO system
    """
    m = blk.model()

    m.fs.ro_properties.set_default_scaling(
        "flow_mass_phase_comp", 1, index=("Liq", "H2O")
    )
    m.fs.ro_properties.set_default_scaling(
        "flow_mass_phase_comp", 1e2, index=("Liq", "NaCl")
    )
    for t in range(1, blk.number_trains + 1):
        train = blk.find_component(f"train_{t}")
        for s in range(1, train.number_stages + 1):
            pump = train.find_component(f"pump{s}")

            set_scaling_factor(pump.control_volume.work, 1e-3)

            ro_stage = train.find_component(f"ro_stage_{s}")

            # Calculate RO scaling factors
            set_scaling_factor(ro_stage.feed_side.length, 1e-1)
            set_scaling_factor(ro_stage.feed_side.width, 1e-3)
            set_scaling_factor(ro_stage.area, 1e-5)
            set_scaling_factor(ro_stage.feed_side.spacer_porosity, 1e-1)
            # set_scaling_factor(ro_stage.feed_side.channel_height, 1e-5)
            for i, x in ro_stage.feed_side.mass_transfer_term.items():
                if i[3] == "NaCl":
                    set_scaling_factor(x, 1e4)
                else:
                    set_scaling_factor(x, 1)
            constraint_scaling_transform(ro_stage.feed_side.eq_dh, 1e-5)
            constraint_scaling_transform(ro_stage.eq_area, 1e-5)
            for i, c in ro_stage.feed_side.eq_K.items():
                set_scaling_factor(c, 1e4)
            # constraint_scaling_transform(ro_stage.feed_side.eq_K, 1e4

    calculate_scaling_factors(blk)


def initialize_ro_system(blk):

    blk.feed.initialize()
    propagate_state(blk.feed_to_feed_splitter)

    blk.feed_splitter.initialize()

    for t in range(1, blk.number_trains + 1):
        splitter_out = blk.feed_splitter.find_component(f"train_{t}_feed")
        train = blk.find_component(f"train_{t}")
        splitter_to_train = train.find_component(f"feed_to_train_{t}")
        propagate_state(splitter_to_train)
        train.feed.initialize()
        propagate_state(train.feed_to_pump1)
        train.pump1.initialize()
        for s in range(1, (train.number_stages + 1)):
            propagate_state(train.find_component(f"pump{s}_to_ro_stage_{s}"))
            stage = train.find_component(f"ro_stage_{s}")
            relax_bounds_for_low_salinity_waters(stage)
            stage.initialize()
            if s != train.number_stages:
                propagate_state(train.find_component(f"ro_stage_{s}_to_pump{s+1}"))
                pump = train.find_component(f"pump{s+1}")
                pump.initialize()
            propagate_state(train.find_component(f"ro_stage_{s}_to_permeate_mixer"))

        train.permeate_mixer.initialize()
        train.permeate.initialize()

        propagate_state(train.last_stage_to_brine)
        train.brine.initialize()

        propagate_state(train.find_component(f"train_{t}_brine_to_mixer"))

    blk.permeate_mixer.initialize()
    propagate_state(blk.permeate_mixer_to_permeate)
    blk.permeate.initialize()

    blk.brine_mixer.initialize()
    propagate_state(blk.brine_mixer_to_brine)
    blk.brine.initialize()


def report_pump(blk, w=30):
    title = "Pump Report"
    side = int(((3 * w) - len(title)) / 2) - 1
    header = "=" * side + f" {title} " + "=" * side
    print(f"\n{header}\n")
    print(f'{"Parameter":<{w}s}{"Value":<{w}s}{"Units":<{w}s}')
    print(f"{'-' * (3 * w)}")

    total_flow = 0
    total_power = 0
    for t in range(1, blk.number_trains + 1):
        train = blk.find_component(f"train_{t}")
        for s in range(1, (train.number_stages + 1)):
            pump = train.find_component(f"pump{s}")
            title = f"Train {t}, Stage {s}"
            side = int(((3 * w) - len(title)) / 2) - 1
            header = "." * side + f" {title} " + "." * side
            if s == 1:
                total_flow += pump.control_volume.properties_out[0].flow_vol
            total_power += pyunits.convert(pump.work_mechanical[0], to_units=pyunits.kW)
            print(f"\n{header}\n")
            print(
                f'{f"Stage {s} Flow In (MGD)":<{w}s}{value(pyunits.convert(pump.control_volume.properties_out[0].flow_vol, to_units=pyunits.Mgallons / pyunits.day)):<{w}.3f}{"MGD"}'
            )
            print(
                f'{f"Stage {s} Flow In (m3/s)":<{w}s}{value(pump.control_volume.properties_out[0].flow_vol):<{w}.3e}{"m3/s"}'
            )
            print(
                f'{f"Stage {s} Flow In (gpm)":<{w}s}{value(pyunits.convert(pump.control_volume.properties_out[0].flow_vol, to_units=pyunits.gallons / pyunits.minute)):<{w}.3f}{"gpm"}'
            )
            print(
                f'{f"Stage {s} Pump Pressure In":<{w}s}{value(pyunits.convert(pump.control_volume.properties_in[0].pressure, to_units=pyunits.psi)):<{w}.1f}{"psi"}'
            )
            print(
                f'{f"Stage {s} Pump Pressure Out":<{w}s}{value(pyunits.convert(pump.control_volume.properties_out[0].pressure, to_units=pyunits.psi)):<{w}.1f}{"psi"}'
            )
            print(
                f'{f"Stage {s} Pump ∆P (psi)":<{w}s}{value(pyunits.convert(pump.control_volume.deltaP[0], to_units=pyunits.psi)):<{w}.1f}{"psi"}'
            )
            print(
                f'{f"Stage {s} Pump ∆P (Pa)":<{w}s}{value(pump.control_volume.deltaP[0]):<{w}.1f}{"Pa"}'
            )
            print(
                f'{f"Stage {s} Pump Work Fluid":<{w}s}{value(pyunits.convert(pump.work_fluid[0], to_units=pyunits.kW)):<{w}.3f}{"kW"}'
            )
            print(
                f'{f"Stage {s} Pump Work Mech.":<{w}s}{value(pyunits.convert(pump.work_mechanical[0], to_units=pyunits.kW)):<{w}.3f}{"kW"}'
            )

    print(f"{'.' * (3 * w)}")
    print(f"{'.' * (3 * w)}")
    print(
        f'{f"Total Flow Rate (MGD)":<{w}s}{value(pyunits.convert(total_flow, to_units=pyunits.Mgallons /pyunits.day)):<{w}.3f}{"MGD"}'
    )
    print(f'{f"Total Flow Rate (m3/s)":<{w}s}{value(total_flow):<{w}.3e}{"m3/s"}')
    print(
        f'{f"Total Flow Rate (gpm)":<{w}s}{value(pyunits.convert(total_flow, to_units=pyunits.gallons / pyunits.minute)):<{w}.3f}{"gpm"}'
    )
    print(
        f'{f"Total Power Consumption":<{w}s}{value(pyunits.convert(total_power, to_units=pyunits.kW)):<{w}.3f}{"kW"}'
    )


if __name__ == "__main__":
    m = build_system(number_trains=4, number_stages=3)
    set_inlet_conditions(m.fs.ro_system, Qin=4*0.154, Cin=0.542)
    set_ro_system_op_conditions(m.fs.ro_system)
    add_ro_scaling(m.fs.ro_system)
    initialize_ro_system(m.fs.ro_system)
    m.fs.obj = Objective(
        expr=m.fs.ro_system.permeate.properties[0].flow_vol_phase["Liq"]
    )
    results = solver.solve(m)
    assert_optimal_termination(results)

    print(f"{iscale.jacobian_cond(m.fs.ro_system):.2e}")
    # m.fs.ro_system.recovery.display()
    m.fs.ro_system.total_power_consumption.display()
    report_pump(m.fs.ro_system, w=40)

    # for t in range(1, m.fs.ro_system.number_trains + 1):
    #     train = m.fs.ro_system.find_component(f"train_{t}")
    #     # print(f"\n--- Train {t} ---")
    #     # train.feed.display()
    #     # train.permeate.display()
    #     # train.brine.display()
    #     for s in range(1, (train.number_stages + 1)):
    # #         if s == train.number_stages:
    #         pump = train.find_component(f"pump{s}")
    #         print(f"\n--- Train {t}, Stage {s} ---")
    #         print(f"Pump {s}- Pressure In: {value(pyunits.convert(pump.control_volume.properties_in[0].pressure, to_units=pyunits.psi))} psi")
    #         print(f"Pump {s}- Pressure Out: {value(pyunits.convert(pump.control_volume.properties_out[0].pressure , to_units=pyunits.psi))} psi")
    #         print(f"Pump {s}- Power: {value(pump.work_mechanical[0] * 1e-3)} kW")
    #             p = pump.control_volume.properties_out[0].pressure.value
    #             pump.control_volume.properties_out[0].pressure.unfix()
    #             pump.control_volume.properties_out[0].pressure.set_value(p * 0.5/0.93)

    #         # stage = train.find_component(f"ro_stage_{s}")
    #         # print(f"\n Stage {s}")
    #         # stage.display()
    # # m.fs.del_component(m.fs.obj
    # m.fs.obj.deactivate()
    # m.fs.obj2 = Objective(
    #     expr=m.fs.ro_system.recovery
    # )
    # print(f"dof = {degrees_of_freedom(m.fs.ro_system)}")
    # m.fs.ro_system.recovery.fix(0.8)
    # print(f"dof = {degrees_of_freedom(m.fs.ro_system)}")
    # results = solver.solve(m)
    # assert_optimal_termination(results)

    # # m.fs.ro_system.display()
