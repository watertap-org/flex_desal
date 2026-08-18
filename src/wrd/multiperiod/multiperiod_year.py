# Model on minute frequency
# Ramp up is associated with period of off-spec water production

# TODO:
# Build basic multiperiod flowsheet
# List of variables that vary operation: Number of RO skids online, pricing chart for 24h
# What should the monthly target be to meet 10,000 AFY?
# Check which components have enough data for pump efficiency

import numpy as np
import pandas as pd
import logging
import matplotlib.pyplot as plt

# Pyomo imports
from pyomo.environ import ConcreteModel, Var, Param, units as pyunits, Objective
from watertap.core.util.model_diagnostics import *
from idaes.core import FlowsheetBlock

# IDAES imports
from idaes.apps.grid_integration.multiperiod.multiperiod import MultiPeriodModel
from idaes.core.solvers.get_solver import get_solver
import idaes.logger as idaeslog

# Based on rates in 2021 from GRIP Cost Tracker
elec_price_invoice_1 = np.array(
    [0.15, 0.16, 0.17, 0.16, 0.16, 0.29, 0.2, 0.21, 0.21, 0.17, 0.16, 0.24]
)

elec_price_invoice_2 = np.array(
    [0.14, 0.14, 0.14, 0.14, 0.14, 0.25, 0.18, 0.2, 0.21, 0.16, 0.15, 0.26]
)

elec_price = elec_price_invoice_1 + elec_price_invoice_2


def build_wrd_flowsheet(m=None, elec_price=0.1):
    if m is None:
        m = ConcreteModel()

    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.energy_per_MG = Param(
        initialize=3,
        units=pyunits.MWh / pyunits.megagallons,
        doc="Energy required per MG of water treated",
    )

    m.fs.acc_prod = Var(
        initialize=0,
        bounds=(0, None),
        units=pyunits.megagallons,
        doc="Accumulate water produced in MG",
    )

    m.fs.pre_acc_prod = Var(
        initialize=0,
        bounds=(0, None),
        units=pyunits.megagallons,
        doc="Accumulate water produced in MG from previous step",
    )

    m.fs.acc_energy = Var(
        initialize=0,
        bounds=(0, None),
        units=pyunits.MWh,
        doc="Accumulate water produced in MG",
    )

    m.fs.pre_acc_energy = Var(
        initialize=0,
        bounds=(0, None),
        units=pyunits.MWh,
        doc="Accumulate water produced in MG from previous step",
    )

    m.fs.water_prod = Var(
        initialize=0,
        bounds=(0, 400),
        units=pyunits.megagallons,
        doc="Water produced in a month in MG",
    )

    m.fs.grid_cost = Var(
        initialize=0,
        bounds=(0, None),
        units=pyunits.dimensionless,
        doc="Grid cost for the month",
    )

    @m.Constraint(doc="Constraint to accumulate water production")
    def eq_acc_water_prod(b):
        return b.fs.acc_prod == b.fs.pre_acc_prod + b.fs.water_prod

    @m.Constraint(doc="Constraint to calculate total energy consumption")
    def eq_acc_energy(b):
        return (
            b.fs.acc_energy
            == b.fs.pre_acc_energy + b.fs.water_prod * b.fs.energy_per_MG
        )

    @m.Constraint(doc="Grid cost")
    def eq_grid_cost(b):
        return b.fs.grid_cost == elec_price * b.fs.water_prod * b.fs.energy_per_MG

    return m


def get_wrd_variable_pairs(t1, t2):
    # Connect the accumulated water produced
    return [
        (t1.fs.acc_prod, t2.fs.pre_acc_prod),
        (t1.fs.acc_energy, t2.fs.pre_acc_energy),
    ]


def unfix_dof(m):
    m.fs.water_prod.unfix()


def create_wrd_mp(
    n_time_points=12,
    elec_price=elec_price,
):
    """
    This function creates a multi-period flowsheet object for each month for the WRD plant. This object contains
    a pyomo model with a block for each time instance.

    Args:
        n_time_points: Number of time blocks to create

    Returns:
        Object containing multi-period vagmd batch flowsheet model
    """
    mp = MultiPeriodModel(
        n_time_points=n_time_points,
        process_model_func=build_wrd_flowsheet,
        linking_variable_func=get_wrd_variable_pairs,
        initialization_func=None,
        unfix_dof_func=unfix_dof,
        outlvl=logging.WARNING,
    )

    """
    Specify the initialization conditions of each period
    """

    flowsheet_options = {
        t: {
            "elec_price": elec_price[t],
        }
        for t in range(n_time_points)
    }

    mp.build_multi_period_model(
        model_data_kwargs=flowsheet_options,
        flowsheet_options=None,
        initialization_options=None,
        unfix_dof_options=None,
    )

    # winter_months = [10, 11]
    # @mp.Constraint(winter_months, doc="November & December has limited access to spreading grounds")
    # def eq_nov_production(b, month):
    #     return b.blocks[month].process.fs.water_prod == 0*pyunits.megagallons

    @mp.Expression(doc="Total cost")
    def total_cost(b):
        return sum([b.blocks[i].process.fs.grid_cost for i in range(n_time_points)])

    @mp.Constraint(doc="Total production")
    def total_production(b):
        return (
            sum([b.blocks[i].process.fs.water_prod for i in range(n_time_points)])
            >= 3258.51 * pyunits.megagallons
        )

    # Set objective
    mp.obj = Objective(expr=mp.total_cost)

    return mp


if __name__ == "__main__":

    n_time_points = 12

    mp = create_wrd_mp(
        n_time_points=n_time_points,
        elec_price=elec_price,
    )
    solver = get_solver()
    results = solver.solve(mp)

    prod = [mp.blocks[i].process.fs.water_prod() for i in range(n_time_points)]
    energy = [
        mp.blocks[i].process.fs.water_prod() * mp.blocks[i].process.fs.energy_per_MG()
        for i in range(n_time_points)
    ]
    demand_peak = [
        mp.blocks[i].process.fs.water_prod() * 5.6 for i in range(n_time_points)
    ]

    months = np.linspace(0, 11, 12)

    fig, ax = plt.subplots(figsize=(5, 5))
    ax.plot(months, prod, label="Water Production", color="blue", marker="o")
    ax.set_xlim(0, 11)

    ax.set_ylim(0, 500)
    ax1 = ax.twinx()

    ax1.plot(months, elec_price, label="Electricity Price", color="black", marker="o")
    ax1.set_ylim(0, 2)

    ax2 = ax.twinx()
    ax2.plot(months, energy, label="Energy Consumption", color="orange", marker="o")
    ax2.spines["right"].set_position(("outward", 55))
    ax2.set_ylim(0, 2500)

    ax2.plot(months, demand_peak, label="Demand Peak", color="green", marker="o")

    handle, label = ax.get_legend_handles_labels()
    handle1, label1 = ax1.get_legend_handles_labels()
    handle2, label2 = ax2.get_legend_handles_labels()

    handles = handle + handle1 + handle2
    labels = label + label1 + label2

    plt.legend(handles=handles, labels=labels)

    ax1.set_ylabel("Electricity Price (2021 $/kWh)")
    ax.set_ylabel("Water production (MG)")
    ax2.set_ylabel("Energy Consumption(MWh)")

    month_names = [
        "Jan",
        "Feb",
        "Mar",
        "Apr",
        "May",
        "Jun",
        "Jul",
        "Aug",
        "Sep",
        "Oct",
        "Nov",
        "Dec",
    ]
    ax.set_xticklabels(month_names)
    ax.xaxis.set_major_locator(plt.MaxNLocator(12))
    ax.set_xlabel("Month")

    fig.tight_layout()
    plt.show()

    print("Total production in MG (Target 3258.1):", mp.total_production())
    print("Total energy cost:", mp.total_cost())
