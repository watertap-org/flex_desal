# Flowsheet that will include pumps and RO unit models for each stage

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

from wrd.components.ro import *
from wrd.components.pump import *
from watertap_contrib.reflo.costing import TreatmentCosting
from watertap_contrib.reflo.flowsheets.KBHDP.utils import solve, calc_scale


def build_ro_system(blk, number_of_stages=1, prop_package=None):
    print(f'\n{"=======> BUILDING RO SYSTEM <=======":^60}\n')

    if prop_package is None:
        m = blk.model()
        prop_package = m.fs.RO_properties

    blk.feed = StateJunction(property_package=prop_package)
    blk.product = StateJunction(property_package=prop_package)
    blk.disposal = StateJunction(property_package=prop_package)
    blk.numberOfStages = Param(initialize=number_of_stages)  # Is Param really needed?
    blk.Stages = RangeSet(blk.numberOfStages)
    # blk.booster_pumps = False

    blk.FirstStage = blk.Stages.first()
    blk.LastStage = blk.Stages.last()
    blk.NonFinalStages = RangeSet(number_of_stages - 1)

    blk.permeate_mixer = Mixer(
        property_package=prop_package,
        inlet_list=[f"ro{i+1}_permeate" for i in range(number_of_stages)],
        has_holdup=False,
        energy_mixing_type=MixingType.extensive,
        momentum_mixing_type=MomentumMixingType.minimize,
    )

    # blk.stage = FlowsheetBlock(RangeSet(number_of_stages), dynamic=False)

    for stage in range(1, number_of_stages + 1):
        # Add a pump based on the characteristics of that stage
        blk.add_component(f"pump{stage}",FlowsheetBlock(dynamic=False))
        build_wrd_pump(
            blk.find_component(f"pump{stage}"), stage, prop_package=prop_package
        ) 

        # Add RO based on that stage
        blk.add_component(f"ro{stage}",FlowsheetBlock(dynamic=False))
        build_wrd_ro(
            blk.find_component(f"ro{stage}"), stage, prop_package=prop_package
        ) 

        # Permeate to mixer
        blk.add_component(
            f"stage{stage}_permeate_to_mixer",
            Arc(
                source=blk.find_component(f"ro{stage}").ro.permeate.outlet, 
                destination=blk.permeate_mixer.find_component(f"ro{stage}_permeate"), # This might be wrong
            ),
        )
        if (
            stage in blk.FirstStage
        ):  # Create Arcs for the 1 stage. Only inlet b/c other components not yet created
            blk.feed_to_pump1 = Arc(
                source=blk.feed.outlet,
                destination=blk.find_component(f"pump{stage}").feed_in.inlet,
            )

            blk.pump1_to_ro1 = Arc(
                source=blk.find_component(
                    f"pump{stage}"
                ).feed_out.outlet,  # Pump exits are feed_out
                destination=blk.find_component(f"pump{stage}").feed.inlet,
            )

        elif stage in blk.NonFinalStages:
            # prev stage to current pump
            blk.add_component(
                f"ro{stage-1}_to_pump{stage}",
                Arc(
                    source=blk.find_component(f"ro{stage-1}").retentate.outlet,
                    destination=blk.find_component(f"pump{stage}").feed_in.inlet,
                ),
            )

            # current pump to current stage
            blk.add_component(
                f"ro{stage}_to_pump{stage}",
                Arc(
                    source=blk.find_component(f"pump{stage}").feed_out.outlet,
                    destination=blk.find_component(f"ro{stage}").feed.inlet,
                ),
            )
        # else:
        # Last stage stuff outside loop

    # Final Pump to final stage, if there are multiple stages
    if number_of_stages > 1:
        blk.add_component(
            f"ro{number_of_stages-1}_to_pump{number_of_stages}",
            Arc(
                source=blk.find_component(f"ro{stage-1}").retentate.outlet,
                destination=blk.find_component(f"pump{stage}").feed_in.inlet,
            ),
        )

    # Final Stage to Brine
    blk.last_stage_retentate_to_ro_retentate = Arc(
        source=blk.find_component(
            f"ro{number_of_stages}"
        ).retentate.outlet,  # Might not be named that...
        destination=blk.disposal.inlet,
    )

    # Permeate Mixer to Product
    blk.permeate_mixer_to_product = Arc(
        source=blk.permeate_mixer.outlet,
        destination=blk.product.inlet,
    )

    # blk.stage_retentate_to_next_stage = Arc( # May be worth understanding this syntax!
    #     blk.NonFinalStages,
    #     rule=lambda blk, n: {
    #         "source": blk.stage[n].retentate.outlet,
    #         "destination": blk.stage[n + 1].feed.inlet,
    #     },
    # )

    # Touch key properties
    blk.feed.properties[0].conc_mass_phase_comp
    blk.product.properties[0].conc_mass_phase_comp
    blk.disposal.properties[0].conc_mass_phase_comp

def scale_ro_system(blk,number_of_stages):
    # Properties. Potentially, this could occur in the treatment train?
    m = blk.model()
    m.fs.ro_properties.set_default_scaling(
        "flow_mass_phase_comp", 1, index=("Liq", "H2O")  # 1e-2 ????
    )
    m.fs.ro_properties.set_default_scaling(
        "flow_mass_phase_comp", 1e2, index=("Liq", "NaCl")
    )

    for stage in range(1,number_of_stages+1):
        add_ro_scaling(blk.find_component(f"ro{stage}"))
        add_pump_scaling(blk.find_component(f"pump{stage}"))
