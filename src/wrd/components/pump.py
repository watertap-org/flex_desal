import os
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
from watertap.unit_models.pressure_changer import Pump
from wrd.components.ro import load_config, get_config_value

from idaes.core.util.scaling import (
    constraint_scaling_transform,
    calculate_scaling_factors,
    set_scaling_factor,
    list_badly_scaled_variables,
    extreme_jacobian_rows,
)

def build_system(**kwargs): # For testing
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.ro_properties = NaClParameterBlock()
    m.fs.pump_system = FlowsheetBlock(dynamic=False)
    build_wrd_pump(m.fs.pump_system, prop_package=m.fs.ro_properties, **kwargs)
    return m


def build_wrd_pump(blk,stage_num=1,prop_pack=None):
    m = blk.model()
    
    if prop_package is None:
        prop_package = m.fs.ro_properties
    
    blk.feed_in = StateJunction(property_package=prop_pack)
    blk.feed_out  = StateJunction(property_package=prop_pack)

    # Get the absolute path of the current script       # Consider moving config to the ro_system, then passing as input
    current_script_path = os.path.abspath(__file__)
    # Get the directory containing the current script
    current_directory = os.path.dirname(current_script_path)
    # Get the parent directory of the current directory (one folder prior)
    parent_directory = os.path.dirname(current_directory)

    config = parent_directory + "/meta_data/wrd_ro_system_inputs.yaml" # Should change ro back and delete the other yaml file (ro_inputs)
    blk.config_data = load_config(config)
    blk.pump = Pump(property_package=prop_pack)    

    # Add Arcs
    blk.feed_in_to_pump = Arc(source=blk.feed_in.outlet,destination=blk.pump.inlet)
    blk.pump_to_feed_out = Arc(source=blk.pump.outlet,destination=blk.feed_out.outlet)
    TransformationFactory("network.expand_arcs").apply_to(blk)  
   # print("Degrees of freedom after adding units:", degrees_of_freedom(blk))

def set_pump_op_conditions(blk,stage_num=1,prop_pack):
    #Configure with input values
    blk.pump.efficiency_pump.fix(
    get_config_value(blk.config_data, "pump_efficiency", "pumps", f"pump_{stage_num}")
    )
    blk.pump.control_volume.properties_out[0].pressure.fix(
    get_config_value(
        blk.config_data, "pump_outlet_pressure", "pumps", f"pump_{i}"
    )
)

def add_pump_scaling(blk):
    set_scaling_factor(blk.pump.work_mechanical[0], 1e-3) # Not sure what value to use here yet
    # Isn't there a needed scaling factor for electricity costs? Where is that scaled?


