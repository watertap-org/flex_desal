from pyomo.environ import ConcreteModel, TransformationFactory, value, units as pyunits
from pyomo.network import Arc

from idaes.core import FlowsheetBlock
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.initialization import propagate_state
import idaes.core.util.scaling as iscale
from idaes.models.unit_models import Feed, Mixer, Separator, StateJunction
from idaes.models.unit_models.mixer import MomentumMixingType
from idaes.models.unit_models.separator import SplittingType

from watertap.property_models.seawater_prop_pack import SeawaterParameterBlock

from srp.utils.utils import touch_flow_and_conc


def build_system(Qin=11343, Cin=1467, split=0.5, feed_temp=27):

    m = ConcreteModel()

    m.Qin = Qin * pyunits.gallons / pyunits.minute
    m.Cin = Cin * pyunits.mg / pyunits.L
    m.feed_temp = feed_temp
    m.split1 = split
    m.split2 = 1 - split

    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = SeawaterParameterBlock()

    m.fs.permeate = Feed(property_package=m.fs.properties)
    m.fs.wells = Feed(property_package=m.fs.properties)

    m.fs.raw_water_tank = FlowsheetBlock(dynamic=False)

    build_raw_water_tank(m.fs.raw_water_tank, prop_package=m.fs.properties)
    set_system_scaling(m)
    set_system_op_conditions(m)
    set_raw_water_tank_op_conditions(m.fs.raw_water_tank)
    init_raw_water_tank(m.fs.raw_water_tank)

    return m


def build_raw_water_tank(
    blk,
    prop_package=None,
    inlet_list=["from_permeate", "from_wells"],
    outlet_list=["to_cooling_tower", "to_service_and_fire"],
):

    print(f'\n{"=======> BUILDING RAW WATER TANK UNIT <=======":^60}\n')

    if prop_package is None:
        prop_package = m.fs.properties

    assert len(inlet_list) == 2

    blk.mixer = Mixer(
        property_package=prop_package,
        momentum_mixing_type=MomentumMixingType.none,
        inlet_list=inlet_list,
    )

    for inlet in inlet_list:
        blk.add_component(f"{inlet}", StateJunction(property_package=prop_package))

    blk.unit = Separator(
        property_package=prop_package,
        outlet_list=outlet_list,
        split_basis=SplittingType.componentFlow,
    )

    # Connect feeds to mixer
    for inlet in inlet_list:
        inlet_port = blk.mixer.find_component(f"{inlet}")
        blk.add_component(
            f"{inlet}_to_unit",
            Arc(
                source=blk.find_component(f"{inlet}").outlet,
                destination=inlet_port,
            ),
        )

    # Connect mixer to separator
    blk.mixer_to_unit = Arc(
        source=blk.mixer.outlet,
        destination=blk.unit.inlet,
    )

    # for outlet in outlet_list:
    #     outlet_port = blk.unit.find_component(f"{outlet}")
    #     blk.add_component(f)

    blk.product = StateJunction(property_package=prop_package)
    blk.unit_to_product = Arc(
        source=blk.unit.outlet,
        destination=blk.product.inlet,
    )

    touch_flow_and_conc(blk.unit)

    TransformationFactory("network.expand_arcs").apply_to(blk)


def set_system_scaling(m):

    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp",
        1 / value(pyunits.convert(m.Qin, to_units=pyunits.m**3 / pyunits.s) * 1000),
        index=("Liq", "H2O"),
    )
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp",
        1 / value(pyunits.convert(m.Qin, to_units=pyunits.m**3 / pyunits.s)),
        index=("Liq", "TDS"),
    )

    iscale.calculate_scaling_factors(m)


def set_system_op_conditions(m, perm_flow_guess=49):

    if perm_flow_guess is None:
        perm_flow_guess = m.Qin * m.split1
    else:
        perm_flow_guess = perm_flow_guess * pyunits.gallons / pyunits.minute
        well_flow_guess = m.Qin - perm_flow_guess

    m.fs.permeate.properties.calculate_state(
        var_args={
            ("flow_vol_phase", ("Liq")): perm_flow_guess,
            ("conc_mass_phase_comp", ("Liq", "TDS")): m.Cin,
            ("pressure", None): 101325,
            ("temperature", None): 273.15 + m.feed_temp,
        },
        hold_state=True,
    )

    m.fs.wells.properties.calculate_state(
        var_args={
            ("flow_vol_phase", ("Liq")): well_flow_guess,
            ("conc_mass_phase_comp", ("Liq", "TDS")): m.Cin,
            ("pressure", None): 101325,
            ("temperature", None): 273.15 + m.feed_temp,
        },
        hold_state=True,
    )


def set_raw_water_tank_op_conditions(blk):
    blk.unit.outlet.pressure[0].fix(101325)


def init_raw_water_tank(blk):

    for inlet in blk.unit.config.inlet_list:
        fb = blk.find_component(f"{inlet}")
        fb.initialize()
        a = blk.find_component(f"{inlet}_to_unit")
        propagate_state(a)


if __name__ == "__main__":
    m = build_system()
    print(f"dof = {degrees_of_freedom(m)}")
    print(m.fs.raw_water_tank.unit.config.inlet_list)
    # m.fs.raw_water_tank.from_wells.outlet.display()
    # # build_raw_water_tank(m.fs.raw_water_tank, prop_package=m.fs.properties)
