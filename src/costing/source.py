import pyomo.environ as pyo
from watertap.costing.util import register_costing_parameter_block
from watertap_contrib.reflo.costing.util import (
    make_capital_cost_var,
    make_variable_operating_cost_var,
)


def build_source_cost_param_block(blk):

    costing = blk.parent_block()
    blk.unit_cost = pyo.Var(
        initialize=0.15,
        units=costing.base_currency / pyo.units.m**3,
        doc="Source cost per cubic meter",
    )


@register_costing_parameter_block(
    build_rule=build_source_cost_param_block, parameter_block_name="source"
)
def cost_source(blk):

    make_capital_cost_var(blk)
    make_variable_operating_cost_var(blk)
    blk.costing_package.add_cost_factor(blk, None)

    blk.capital_cost_constraint = pyo.Constraint(
        expr=blk.capital_cost == 0 * blk.costing_package.base_currency
    )

    blk.base_period_flow = pyo.units.convert(
        blk.unit_model.properties[0].flow_vol_phase["Liq"]
        * blk.costing_package.utilization_factor,
        to_units=pyo.units.m**3 / blk.costing_package.base_period,
    )

    blk.variable_operating_cost_constraint = pyo.Constraint(
        expr=blk.variable_operating_cost
        == pyo.units.convert(
            blk.costing_package.source.unit_cost * blk.base_period_flow,
            to_units=blk.costing_package.base_currency
            / blk.costing_package.base_period,
        )
    )
