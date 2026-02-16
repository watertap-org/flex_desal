from idaes.apps.grid_integration import OperationModel
from pyomo.environ import (
    Constraint,
    NonNegativeReals,
    Param,
    RangeSet,
    Var,
    exp,
    units as pyunits,
)
from pricetaker.flowsheets import params as um_params
from pricetaker.flowsheets.unit_models import _add_required_variables


def ro_skid_operation_model(blk, params: um_params.WRD_ROParams):
    """
    Builds operation model for an RO skid

    Parameters
    ----------
    blk : OperationModel
        IDAES OperationModel instance

    params : object
        Input parameters needed for the model
    """
    _add_required_variables(blk)
    blk.coeffs = Param(["a", "b"], initialize=params.surrogate_coeffs)

    blk.operational_limits = Constraint(
        expr=blk.feed_flowrate == blk.op_mode * params.nominal_flowrate
    )

    if params.surrogate_type == "constant_energy_intensity":
        blk.calculate_energy_intensity = Constraint(
            expr=blk.energy_intensity
            == (blk.coeffs["a"] + blk.coeffs["b"] * blk.feed_flowrate),
            doc="Calculates the specific energy requirement",
        )
    else:
        raise ValueError("Unrecognized surrogate type")


def wrd_reverse_osmosis_operation_model(blk, params: um_params.WRD_ROParams):
    """
    Builds operation model for the reverse osmosis unit

    Parameters
    ----------
    blk : OperationModel
        IDAES OperationModel instance

    params : object
        Input parameters needed for the model
    """
    # Declare required variables
    _add_required_variables(blk)
    blk.inlet_flowrate = Var(within=NonNegativeReals, units=pyunits.m**3 / pyunits.hr)
    # Defining a slack variable for flowrate that is not accounted
    # for by the sum of RO intake pumps
    blk.leftover_flow = Var(within=NonNegativeReals, units=pyunits.m**3 / pyunits.hr)

    # Build RO skid models
    blk.set_ro_skids = RangeSet(params.num_ro_skids)
    blk.ro_skid = OperationModel(
        blk.set_ro_skids,
        model_func=ro_skid_operation_model,
        model_args={"params": params},
        minimum_up_time=params.minimum_uptime,
        minimum_down_time=params.minimum_downtime,
    )

    # Remove overall mass balance and power consumption calculation
    blk.del_component(blk.recovery)
    blk.del_component(blk.energy_intensity)
    blk.del_component(blk.mass_balance)
    blk.del_component(blk.calculate_product_flowrate)
    blk.del_component(blk.calculate_power_consumption)

    # Declare required constraints
    blk.calculate_leftover_flow = Constraint(
        expr=blk.feed_flowrate == blk.inlet_flowrate + blk.leftover_flow,
        doc="Calculates leftover flowrate",
    )
    blk.feed_mass_balance = Constraint(
        expr=blk.inlet_flowrate
        == sum(blk.ro_skid[i].feed_flowrate for i in blk.set_ro_skids),
        doc="Mass balance at the feed",
    )
    blk.product_mass_balance = Constraint(
        expr=blk.product_flowrate
        == sum(blk.ro_skid[i].product_flowrate for i in blk.set_ro_skids),
        doc="Mass balance on permeate side",
    )
    blk.reject_mass_balance = Constraint(
        expr=blk.reject_flowrate
        == sum(blk.ro_skid[i].reject_flowrate for i in blk.set_ro_skids),
        doc="Mass balance on brine side",
    )
    blk.calculate_power_consumption = Constraint(
        expr=blk.power_consumption
        == sum(blk.ro_skid[i].power_consumption for i in blk.set_ro_skids),
        doc="Calculates the total power requirement for RO",
    )

    # symmetry breaking for >1 skid. Skids can only operate if the previous skid is on
    @blk.Constraint(blk.set_ro_skids)
    def symmetry_breaking_cuts(b, index):
        if index == 1:
            return Constraint.Skip
        return b.ro_skid[index].op_mode <= b.ro_skid[index - 1].op_mode

    # Ensure that the operation of minimum number of skids is identical
    blk.set_min_operating_skids = RangeSet(2, params.minimum_operating_skids)

    @blk.Constraint(blk.set_min_operating_skids)
    def minimum_ro_skids_startup(b, index):
        return b.ro_skid[index].startup == b.ro_skid[1].startup

    @blk.Constraint(blk.set_min_operating_skids)
    def minimum_ro_skids_op_mode(b, index):
        return b.ro_skid[index].op_mode == b.ro_skid[1].op_mode

    @blk.Constraint(blk.set_min_operating_skids)
    def minimum_ro_skids_shutdown(b, index):
        return b.ro_skid[index].shutdown == b.ro_skid[1].shutdown

    # Update bounds on recovery and energy intensity for all skids
    ei_lb, ei_ub = params.get_energy_intensity_bounds()
    for skid in blk.set_ro_skids:
        blk.ro_skid[skid].feed_flowrate.setlb(params.minimum_flowrate)
        blk.ro_skid[skid].feed_flowrate.setub(params.maximum_flowrate)
        blk.ro_skid[skid].energy_intensity.setlb(ei_lb)
        blk.ro_skid[skid].energy_intensity.setub(ei_ub)


# Currently implementing UF same way as RO.
# However, this will increase the decision variables significantly, increasing solve time.


def uf_pump_operation_model(blk, params: um_params.WRD_UFParams):
    """
    Builds operation model for a UF pump

    Parameters
    ----------
    blk : OperationModel
        IDAES OperationModel instance

    params : object
        Input parameters needed for the model
    """
    _add_required_variables(blk)
    blk.coeffs = Param(["a", "b"], initialize=params.surrogate_coeffs)

    blk.operational_limits = Constraint(
        expr=blk.feed_flowrate == blk.op_mode * params.nominal_flowrate
    )

    if params.surrogate_type == "constant_energy_intensity":
        blk.calculate_energy_intensity = Constraint(
            expr=blk.energy_intensity
            == (blk.coeffs["a"] + blk.coeffs["b"] * blk.feed_flowrate),
            doc="Calculates the specific energy requirement",
        )
    else:
        raise ValueError("Unrecognized surrogate type")


def wrd_uf_operation_model(blk, params: um_params.WRD_UFParams):
    """
    Builds operation model for UF unit in WRD case

    Parameters
    ----------
    blk : OperationModel
        IDAES OperationModel instance

    params : object
        Input parameters needed for the model
    """
    # Declare required variables
    _add_required_variables(blk)
    blk.inlet_flowrate = Var(within=NonNegativeReals, units=pyunits.m**3 / pyunits.hr)
    # Defining a slack variable for flowrate that is not accounted
    # for by the sum of RO intake pumps
    blk.leftover_flow = Var(within=NonNegativeReals, units=pyunits.m**3 / pyunits.hr)

    # Build RO skid models
    blk.set_uf_pumps = RangeSet(params.num_uf_pumps)
    blk.uf_pumps = OperationModel(
        blk.set_uf_pumps,
        model_func=uf_pump_operation_model,
        model_args={"params": params},
        minimum_up_time=params.minimum_uptime,
        minimum_down_time=params.minimum_downtime,
    )

    # Remove overall mass balance and power consumption calculation
    blk.del_component(blk.recovery)
    blk.del_component(blk.energy_intensity)
    blk.del_component(blk.mass_balance)
    blk.del_component(blk.calculate_product_flowrate)
    blk.del_component(blk.calculate_power_consumption)

    # Declare required constraints
    blk.calculate_leftover_flow = Constraint(
        expr=blk.feed_flowrate == blk.inlet_flowrate + blk.leftover_flow,
        doc="Calculates leftover flowrate",
    )
    blk.feed_mass_balance = Constraint(
        expr=blk.inlet_flowrate
        == sum(blk.uf_pumps[i].feed_flowrate for i in blk.set_uf_pumps),
        doc="Mass balance at the feed",
    )
    blk.product_mass_balance = Constraint(
        expr=blk.product_flowrate
        == sum(blk.uf_pumps[i].product_flowrate for i in blk.set_uf_pumps),
        doc="Mass balance on permeate side",
    )
    blk.reject_mass_balance = Constraint(
        expr=blk.reject_flowrate
        == sum(blk.uf_pumps[i].reject_flowrate for i in blk.set_uf_pumps),
        doc="Mass balance on brine side",
    )
    blk.calculate_power_consumption = Constraint(
        expr=blk.power_consumption
        == sum(blk.uf_pumps[i].power_consumption for i in blk.set_uf_pumps),
        doc="Calculates the total power requirement for RO",
    )

    # symmetry breaking for >1 skid. Skids can only operate if the previous skid is on
    # This is not true in practice, but which exact pump is on shouldn't matter
    @blk.Constraint(blk.set_uf_pumps)
    def symmetry_breaking_cuts(b, index):
        if index == 1:
            return Constraint.Skip
        return b.uf_pumps[index].op_mode <= b.uf_pumps[index - 1].op_mode

    # Ensure that the operation of minimum number of skids is identical
    blk.set_min_operating_pumps = RangeSet(2, params.minimum_operating_pumps)

    @blk.Constraint(blk.set_min_operating_pumps)
    def minimum_uf_pumps_startup(b, index):
        return b.uf_pumps[index].startup == b.uf_pumps[1].startup

    @blk.Constraint(blk.set_min_operating_pumps)
    def minimum_uf_pumps_op_mode(b, index):
        return b.uf_pumps[index].op_mode == b.uf_pumps[1].op_mode

    @blk.Constraint(blk.set_min_operating_pumps)
    def minimum_uf_pumps_shutdown(b, index):
        return b.uf_pumps[index].shutdown == b.uf_pumps[1].shutdown

    # Update bounds on recovery and energy intensity for all skids
    ei_lb, ei_ub = params.get_energy_intensity_bounds()
    for pump in blk.set_uf_pumps:
        blk.uf_pumps[pump].feed_flowrate.setlb(params.minimum_flowrate)
        blk.uf_pumps[pump].feed_flowrate.setub(params.maximum_flowrate)
        blk.uf_pumps[pump].energy_intensity.setlb(ei_lb)
        blk.uf_pumps[pump].energy_intensity.setub(ei_ub)
