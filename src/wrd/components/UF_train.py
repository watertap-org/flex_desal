from pyomo.environ import (
    ConcreteModel,
    Expression,
    value,
    assert_optimal_termination,
    units as pyunits,
    value,
    Set,
    TransformationFactory,
)
from pyomo.network import Arc

from idaes.core.util.initialization import propagate_state
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core import FlowsheetBlock
from idaes.models.unit_models import (
    Feed,
    Product,
    MomentumMixingType,
    Mixer,
    StateJunction,
)
from idaes.core.util.scaling import calculate_scaling_factors

from watertap.costing import WaterTAPCosting
from watertap.core.util.model_diagnostics.infeasible import *
from watertap.property_models.NaCl_T_dep_prop_pack import NaClParameterBlock
from watertap.core.solvers import get_solver

from wrd.utilities import load_config, get_config_file
from wrd.components.pump import *
from wrd.components.UF_separator import *
from srp.utils import touch_flow_and_conc

__all__ = [
    "build_ro_train",
    "set_ro_train_op_conditions",
    "set_ro_train_scaling",
    "initialize_ro_train",
    "report_ro_train",
    "add_ro_train_costing",
]

solver = get_solver()


def build_system(file="wrd_uf_inputs_8_19_21.yaml"):
    # Will want to combine all inputs into one yaml instead of having separate ones
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = NaClParameterBlock()

    m.fs.uf_train = FlowsheetBlock(dynamic=False)
    build_uf_train(
        m.fs.uf_train, prop_package=m.fs.properties, file=file
    )

    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e-1, index=("Liq", "H2O")  # changed from 1
    )
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e2, index=("Liq", "NaCl")
    )
    return m

def build_uf_train(blk, file="wrd_ro_inputs_8_19_21.yaml", prop_package=None):
    if prop_package is None:
        m = blk.model()
        prop_package = m.fs.properties

    blk.config_data = load_config(get_config_file(file))

    blk.feed = StateJunction(property_package=prop_package)
    touch_flow_and_conc(blk.feed)

    blk.product = StateJunction(property_package=prop_package)
    touch_flow_and_conc(blk.product)

    blk.disposal = StateJunction(property_package=prop_package)
    touch_flow_and_conc(blk.disposal)

    blk.pump = FlowsheetBlock(dynamic=False)
    build_pump(blk.pump,file=file,prop_package=prop_package)

    blk.UF = FlowsheetBlock(dynamic=False)
    build_separator(blk.UF,prop_package=prop_package,outlet_list=["product","disposal"])

    blk.feed_to_pump = Arc(source=blk.feed.outlet, destination=blk.pump.feed.inlet)
    blk.pump_to_UF = Arc(source=blk.pump.product.outlet, destination=blk.UF.feed.inlet)
    blk.UF_to_disposal = Arc(source=blk.UF.disposal.outlet, destination=blk.disposal.inlet)
    blk.UF_to_product = Arc(source=blk.UF.product.outlet, destination=blk.product.inlet)

    TransformationFactory("network.expand_arcs").apply_to(blk)

# def report_ro_train(blk, train_num=None, w=30):

#     if train_num is None:
#         title = "RO Train Report"
#         title2 = f"Overall Train Performance"
#     else:
#         title = f"RO Train {train_num} Report"
#         title2 = f"Overall Train {train_num} Performance"

#     side = int(((3 * w) - len(title)) / 2) - 1
#     header = "=" * side + f" {title} " + "=" * side
#     print(f"\n{header}\n")

#     for i in blk.stages:
#         title = f"Stage {i}"
#         side = int(((3 * w) - len(title)) / 2) - 1
#         header = "_" * side + f" {title} " + "_" * side
#         print(f"\n\n{header}\n")
#         print(
#             f'{f"Stage {i} Feed Flow":<{w}s}{value(pyunits.convert(blk.stage[i].feed.properties[0].flow_vol_phase["Liq"], to_units=pyunits.gallons / pyunits.minute)):<{w}.3f}{"gpm"}'
#         )
#         print(
#             f'{f"Stage {i} Feed Conc.":<{w}s}{value(pyunits.convert(blk.stage[i].feed.properties[0].conc_mass_phase_comp["Liq", "NaCl"], to_units=pyunits.mg / pyunits.L)):<{w}.3f}{"mg/L"}'
#         )
#         report_ro_stage(blk.stage[i], w=w)

#     side = int(((3 * w) - len(title2)) / 2) - 1
#     header = "(" * side + f" {title2} " + ")" * side
#     print(f"\n\n{header}\n")
#     print(f'{"Parameter":<{w}s}{"Value":<{w}s}{"Units":<{w}s}')
#     print(f"{'-' * (3 * w)}")
#     print(
#         f'{f"Total Pump Power":<{w}s}{value(pyunits.convert(blk.total_pump_power, to_units=pyunits.kilowatt)):<{w}.3f}{"kW"}'
#     )

#     print(f'{f"Overall Recovery":<{w}s}{value(blk.recovery_vol)*100:<{w}.3f}{"%"}')
#     print(
#         f'{f"Total Feed Flow":<{w}s}{value(pyunits.convert(blk.feed.properties[0].flow_vol_phase["Liq"], to_units=pyunits.gallons / pyunits.minute)):<{w}.3f}{"gpm"}'
#     )
#     print(
#         f'{f"Inlet Feed Conc":<{w}s}{value(pyunits.convert(blk.feed.properties[0].conc_mass_phase_comp["Liq", "NaCl"], to_units=pyunits.mg / pyunits.L)):<{w}.3f}{"mg/L"}'
#     )
#     print(
#         f'{f"Total Perm Flow":<{w}s}{value(pyunits.convert(blk.product.properties[0].flow_vol_phase["Liq"], to_units=pyunits.gallons / pyunits.minute)):<{w}.3f}{"gpm"}'
#     )
#     print(
#         f'{f"Final Perm Conc":<{w}s}{value(pyunits.convert(blk.product.properties[0].conc_mass_phase_comp["Liq", "NaCl"], to_units=pyunits.mg / pyunits.L)):<{w}.3f}{"mg/L"}'
#     )
#     print(
#         f'{f"Total Brine Flow":<{w}s}{value(pyunits.convert(blk.disposal.properties[0].flow_vol_phase["Liq"], to_units=pyunits.gallons / pyunits.minute)):<{w}.3f}{"gpm"}'
#     )
#     print(
#         f'{f"Final Brine Conc":<{w}s}{value(pyunits.convert(blk.disposal.properties[0].conc_mass_phase_comp["Liq", "NaCl"], to_units=pyunits.mg / pyunits.L)):<{w}.3f}{"mg/L"}'
#     )
#     print()
#     for i, inlet in enumerate(blk.mixer.config.inlet_list, 1):
#         sb = blk.mixer.find_component(f"{inlet}_state")
#         print(
#             f'{f"  Stage {i} Feed Flow":<{w}s}{value(pyunits.convert(blk.stage[i].feed.properties[0].flow_vol_phase["Liq"], to_units=pyunits.gallons / pyunits.minute)):<{w}.3f}{"gpm"}'
#         )
#         print(
#             f'{f"  Stage {i} Feed Conc":<{w}s}{value(pyunits.convert(blk.stage[i].feed.properties[0].conc_mass_phase_comp["Liq", "NaCl"], to_units=pyunits.mg / pyunits.L)):<{w}.3f}{"mg/L"}'
#         )
#         print(
#             f'{f"  Stage {i} Perm Flow":<{w}s}{value(pyunits.convert(sb[0].flow_vol_phase["Liq"], to_units=pyunits.gallons / pyunits.minute)):<{w}.3f}{"gpm"}'
#         )
#         print(
#             f'{f"  Stage {i} Perm Conc":<{w}s}{value(pyunits.convert(sb[0].conc_mass_phase_comp["Liq", "NaCl"], to_units=pyunits.mg / pyunits.L)):<{w}.3f}{"mg/L"}'
#         )



def main(
    Qin=2637,
    Cin=0.528,
    Tin=302,
    Pin=101325,
    file="wrd_ro_inputs_8_19_21.yaml",
    add_costing=True,
):

    m = build_system(file=file)
    # set_ro_train_scaling(m.fs.ro_train)
    # calculate_scaling_factors(m)
    # set_inlet_conditions(m, Qin=Qin, Cin=Cin, Tin=Tin, Pin=Pin)
    # set_ro_train_op_conditions(m.fs.ro_train)

    # initialize_system(m)
    # assert degrees_of_freedom(m) == 0
    # results = solver.solve(m, tee=True)
    # assert_optimal_termination(results)

    # if add_costing:
    #     m.fs.costing = WaterTAPCosting()
    #     add_ro_train_costing(m.fs.ro_train)
    #     m.fs.costing.cost_process()
    #     m.fs.costing.add_LCOW(m.fs.product.properties[0].flow_vol_phase["Liq"])
    #     m.fs.costing.add_specific_energy_consumption(
    #         m.fs.product.properties[0].flow_vol_phase["Liq"],
    #         name="SEC",
    #     )
    #     m.fs.costing.initialize()

    #     assert degrees_of_freedom(m) == 0
    #     results = solver.solve(m, tee=True)
    #     assert_optimal_termination(results)

    # report_ro_train(m.fs.ro_train, w=30)

    return m

if __name__ == "__main__":
    m = main()
