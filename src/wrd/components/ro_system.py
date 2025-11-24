"""
This flowsheet models one train, calling ro_stage to build each train. 
The outputs could be multiplied by 2,3, or 4 to simulate operation of different numbers of trains. 

# Feed
# Feed Pump
# Stage 1 (PRO1)
# Interstage Booster Pump 1
# Stage 2 (PRO2)
# Interstage Booster Pump 2
# Stage Three (TRO)
# Mixer for each stage's permeate
# Permeate
# Brine
"""

from pyomo.environ import (
    ConcreteModel,
    Param,
    Var,
    Constraint,
    NonNegativeReals,
    TransformationFactory,
    RangeSet,
    assert_optimal_termination,
    value,
    units as pyunits,
)
from pyomo.network import Arc

from idaes.core import FlowsheetBlock, UnitModelCostingBlock
from idaes.core.util.initialization import propagate_state
from idaes.models.unit_models import Product, Feed, StateJunction
from idaes.models.unit_models.mixer import (
    Mixer,
    MomentumMixingType,
)
from idaes.core.util.model_statistics import degrees_of_freedom
import idaes.core.util.scaling as iscale

from watertap.property_models.NaCl_prop_pack import NaClParameterBlock

from wrd.components.ro_stage import *
from wrd.components.pump import *
from watertap_contrib.reflo.costing import TreatmentCosting
from watertap_contrib.reflo.flowsheets.KBHDP.utils import solve, calc_scale


def build_system(**kwargs):
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.ro_properties = NaClParameterBlock()
    m.fs.ro_system = FlowsheetBlock(dynamic=False)
    build_wrd_ro_system(m.fs.ro_system, prop_package=m.fs.ro_properties, **kwargs)
    return m


def build_wrd_ro_system(blk, number_of_stages=1, prop_package=None):
    print(f'\n{"=======> BUILDING RO SYSTEM <=======":^60}\n')

    if prop_package is None:
        m = blk.model()
        prop_package = m.fs.RO_properties

    blk.feed = StateJunction(property_package=prop_package)
    blk.product = StateJunction(property_package=prop_package)
    blk.disposal = StateJunction(property_package=prop_package)
    blk.numberOfStages = Param(initialize=number_of_stages)  # Is Param really needed?
    # blk.Stages = RangeSet(blk.numberOfStages)

    # blk.FirstStage = blk.Stages.first()
    # blk.LastStage = blk.Stages.last()
    # blk.NonFinalStages = RangeSet(number_of_stages - 1)
    # blk.stage = FlowsheetBlock(RangeSet(number_of_stages), dynamic=False) # I think I understand this, but opted to do it differently

    blk.permeate_mixer = Mixer(
        property_package=prop_package,
        inlet_list=[f"ro{i+1}_permeate" for i in range(number_of_stages)],
        has_holdup=False,
        energy_mixing_type=MixingType.extensive,
        momentum_mixing_type=MomentumMixingType.minimize,
    )

    for stage in range(1, number_of_stages + 1):
        # Add a pump based on the characteristics of that stage
        blk.add_component(f"pump{stage}", FlowsheetBlock(dynamic=False))
        build_wrd_pump(
            blk.find_component(f"pump{stage}"), stage, prop_package=prop_package
        )

        # Add RO based on that stage
        blk.add_component(f"ro{stage}", FlowsheetBlock(dynamic=False))
        build_wrd_ro(blk.find_component(f"ro{stage}"), stage, prop_package=prop_package)

        # Somewhere need to add the overall srecovery, like in ro_multi_skid

        # Permeate to mixer
        blk.add_component(
            f"stage{stage}_permeate_to_mixer",
            Arc(
                source=blk.find_component(f"ro{stage}").ro.permeate,
                destination=blk.permeate_mixer.find_component(
                    f"ro{stage}_permeate"
                ),  # This might be wrong
            ),
        )
        if stage == 1:
            # Create Arcs for the 1 stage. Only inlet b/c other components not yet created
            blk.feed_to_pump1 = Arc(
                source=blk.feed.outlet,
                destination=blk.find_component(f"pump{stage}").feed_in.inlet,
            )

            blk.pump1_to_ro1 = Arc(
                source=blk.find_component(
                    f"pump{stage}"
                ).feed_out.outlet,  # Pump exits are feed_out
                destination=blk.find_component(f"ro{stage}").feed.inlet,
            )

        elif stage != number_of_stages:  # Not the final stage
            # prev stage to current pump
            blk.add_component(
                f"ro{stage-1}_to_pump{stage}",
                Arc(
                    source=blk.find_component(f"ro{stage-1}").retentate,
                    destination=blk.find_component(f"pump{stage}").feed_in.inlet,
                ),
            )

            # current pump to current stage
            blk.add_component(
                f"pump{stage}_to_ro{stage}",
                Arc(
                    source=blk.find_component(f"pump{stage}").feed_out.outlet,
                    destination=blk.find_component(f"ro{stage}").feed,
                ),
            )
        # else:
        # Last stage stuff outside loop

    # Final Pump to final stage, if there are multiple stages
    if number_of_stages > 1:
        blk.add_component(
            f"ro{number_of_stages-1}_to_pump{number_of_stages}",
            Arc(
                source=blk.find_component(f"ro{stage-1}").retentate,
                destination=blk.find_component(f"pump{stage}").feed_in.inlet,
            ),
        )

    # Final Stage to Brine
    blk.last_stage_retentate_to_disposal = Arc(
        source=blk.find_component(f"ro{number_of_stages}").retentate.outlet,
        destination=blk.disposal.inlet,
    )

    # Permeate Mixer to Product
    blk.permeate_mixer_to_product = Arc(
        source=blk.permeate_mixer.outlet,
        destination=blk.product.inlet,
    )
    TransformationFactory("network.expand_arcs").apply_to(blk)

    blk.recovery = Var(
        initialize=0.5,
        bounds=(0, 1),
        units=pyunits.dimensionless,
        doc="Overall recovery",
    )

    @blk.Constraint(doc="Overall recovery constraint")
    def eq_recovery(b):
        return (
            b.recovery
            == b.product.properties[0].flow_vol_phase["Liq"]
            / b.feed.properties[0].flow_vol_phase["Liq"]
        )

    # Touch key properties
    blk.feed.properties[0].conc_mass_phase_comp
    blk.product.properties[0].conc_mass_phase_comp
    blk.disposal.properties[0].conc_mass_phase_comp


def set_ro_system_op_conditions(blk, number_of_stages=1):
    # Could load and pass config data here?
    for stage in range(1, number_of_stages + 1):
        set_pump_op_conditions(blk.find_component(f"pump{stage}"))
        set_ro_op_conditions(blk.find_component(f"ro{stage}"))


def set_inlet_conditions(blk, Qin=0.154, Cin=0.542, P_in=1):
    """
    Set the operation conditions for the entire RO System
    """
    Qin = (Qin) * pyunits.m**3 / pyunits.s  # Feed flow rate in m3/s
    Cin = Cin * pyunits.g / pyunits.L  # Feed concentration in g/L
    rho = 1000 * pyunits.kg / pyunits.m**3  # Approximate density of water
    feed_mass_flow_water = Qin * rho
    feed_mass_flow_salt = Cin * Qin

    blk.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"].fix(feed_mass_flow_water)
    blk.feed.properties[0].flow_mass_phase_comp["Liq", "NaCl"].fix(feed_mass_flow_salt)
    blk.feed.properties[0].temperature.fix(298.15 * pyunits.K)  # 25 C
    blk.feed.properties[0].pressure.fix(P_in * pyunits.bar)


def add_ro_system_scaling(blk, number_of_stages=1):
    # Properties. Potentially, this could occur in the treatment train?
    m = blk.model()
    m.fs.ro_properties.set_default_scaling(
        "flow_mass_phase_comp", 1e-1, index=("Liq", "H2O")
    )
    m.fs.ro_properties.set_default_scaling(
        "flow_mass_phase_comp", 1e2, index=("Liq", "NaCl")
    )

    for stage in range(1, number_of_stages + 1):
        add_ro_scaling(blk.find_component(f"ro{stage}"))
        add_pump_scaling(blk.find_component(f"pump{stage}"))

    calculate_scaling_factors(m)


def initialize_ro_system(blk, number_of_stages=1):

    blk.feed.initialize()

    for stage in range(1, number_of_stages + 1):
        if stage == 1:
            propagate_state(blk.feed_to_pump1)
        else:
            propagate_state(blk.find_component(f"ro{stage-1}_to_pump{stage}"))

        initialize_pump(blk.find_component(f"pump{stage}"))
        propagate_state(blk.find_component(f"pump{stage}_to_ro{stage}"))
        initialize_ro(blk.find_component(f"ro{stage}"))
        propagate_state(blk.find_component(f"ro{stage}_permeate_to_mixer"))

    # Final stage to products
    blk.permeate_mixer.initialize()
    propagate_state(blk.permeate_mixer_to_product)
    blk.product.initialize()

    propagate_state(blk.last_stage_retentate_to_disposal)
    blk.disposal.initialize()


def report_ro_system(blk, w=30):
    title = "Pump Report"
    side = int(((3 * w) - len(title)) / 2) - 1
    header = "=" * side + f" {title} " + "=" * side
    print(f"\n{header}\n")
    print(f'{"Parameter":<{w}s}{"Value":<{w}s}{"Units":<{w}s}')
    print(f"{'-' * (3 * w)}")

    total_flow = blk.feed_in.properties[0].flow_vol
    overall_recovery = blk.recovery
    deltaP = (
        blk.product.properties[0].pressure - blk.feed.properties[0].pressure
    )  # Probably not the properties to report
    print(
        f'{f"Total Flow Rate (MGD)":<{w}s}{value(pyunits.convert(total_flow, to_units=pyunits.Mgallons /pyunits.day)):<{w}.3f}{"MGD"}'
    )
    print(f'{f"Total Flow Rate (m3/s)":<{w}s}{value(total_flow):<{w}.3e}{"m3/s"}')
    print(
        f'{f"Total Flow Rate (gpm)":<{w}s}{value(pyunits.convert(total_flow, to_units=pyunits.gallons / pyunits.minute)):<{w}.3f}{"gpm"}'
    )
    print(f'{f"Overall Recovery (-)":<{w}s}{value(overall_recovery):<{w}.1e}{"-"}')
    print(f'{f"Pressure Change (Pa)":<{w}s}{value(deltaP):<{w}.3e}{"Pa"}')
    print(
        f'{f"Pressure Change (bar)":<{w}s}{value(pyunits.convert(deltaP,to_units=pyunits.bar)):<{w}.3e}{"bar"}'
    )


if __name__ == "__main__":
    m = build_system(number_of_stages=1)  # should test 1 stage and 3 stage
    print(f"{degrees_of_freedom(m)} degrees of freedom after build")
    set_inlet_conditions(m.fs.ro_system, Qin=0.154, Cin=0.542, P_in=1)
    set_ro_system_op_conditions(m.fs.ro_system, number_of_stages=1)
    print(f"{degrees_of_freedom(m)} degrees of freedom after setting op conditions")
    add_ro_system_scaling(m.fs.ro_system, number_of_stages=1)
    initialize_ro_system(m.fs.ro_system, number_of_stages=1)
    m.fs.obj = Objective(
        expr=m.fs.ro_system.product.properties[0].flow_vol_phase["Liq"]
    )  # There is no D.o.f to optimize with
    solver = get_solver()
    results = solver.solve(m)
    assert_optimal_termination(results)

    report_ro_system(m.fs.ro_system)
