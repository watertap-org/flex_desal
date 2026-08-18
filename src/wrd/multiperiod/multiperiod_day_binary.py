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
from pyomo.environ import (
    ConcreteModel,
    Var,
    Param,
    units as pyunits,
    Objective,
    Binary,
)
import matplotlib.dates as mdates

from idaes.core import FlowsheetBlock
from watertap.core.util.model_diagnostics import *
from idaes.core.util.model_statistics import *

# IDAES imports
from idaes.apps.grid_integration.multiperiod.multiperiod import MultiPeriodModel
from idaes.core.solvers.get_solver import get_solver
import idaes.logger as idaeslog
from pyomo.environ import SolverFactory

# Based on rates in 2021 from GRIP Cost Tracker
# Based on rates in 2021 from GRIP Cost Tracker
elec_price_invoice_1 = [
    0.15,
    0.16,
    0.17,
    0.16,
    0.16,
    0.29,
    0.2,
    0.21,
    0.21,
    0.17,
    0.16,
    0.24,
]

elec_price_invoice_2 = [
    0.14,
    0.14,
    0.14,
    0.14,
    0.14,
    0.25,
    0.18,
    0.2,
    0.21,
    0.16,
    0.15,
    0.26,
]

elec_price = elec_price_invoice_1 + elec_price_invoice_2

elec_price = elec_price_invoice_1 + elec_price_invoice_2


def build_elec_price_summer():
    # Delivery Pricing $/kWh
    on_peak_del = 0.01971
    mid_peak_del = 0.01971
    off_peak_del = 0.01971
    super_off_peak_del = 0

    # Generation Pricing $/kWh
    on_peak_gen = 0.09934
    mid_peak_gen = 0.0891
    off_peak_gen = 0.05782
    super_off_peak_gen = 0

    elec_price = np.ones(25)

    # off peak 12 AM - 4 PM
    elec_price[0:16] = off_peak_del + off_peak_gen
    # one peak 4 PM - 9 PM
    elec_price[16:22] = on_peak_del + mid_peak_gen
    # off peak 9 PM - 12 AM
    elec_price[22:25] = off_peak_del + off_peak_gen

    return elec_price


def build_elec_price_winter():
    # Delivery Pricing
    on_peak_del = 0
    mid_peak_del = 0.0239
    off_peak_del = 0.0228
    super_off_peak_del = 0.022137

    # Generation Pricing $/kWh
    on_peak_gen = 0
    mid_peak_gen = 0.07663
    off_peak_gen = 0.06397
    super_off_peak_gen = 0.04026

    elec_price = np.ones(25)

    # Off peak 12 AM - 8 AM
    elec_price[0:8] = off_peak_del + off_peak_gen
    # Super off peak 8 AM - 4 PM
    elec_price[8:16] = super_off_peak_del + super_off_peak_gen
    # Mid peak 4 PM - 9 PM
    elec_price[16:21] = mid_peak_del + mid_peak_gen
    # Off peak peak 9 PM - 12 AM
    elec_price[21:25] = off_peak_del + off_peak_gen

    return elec_price


def build_wrd_flowsheet(m=None, elec_price=0.1):
    if m is None:
        m = ConcreteModel()

    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.energy_per_MG = Param(
        initialize=10,
        units=pyunits.kWh / pyunits.megagallons,
        doc="Energy required per MG",
    )

    m.fs.flow_100 = Var(initialize=0, domain=Binary, doc="binary?")

    m.fs.flow_50 = Var(initialize=1, domain=Binary, doc="binary?")

    m.fs.acc_prod = Var(
        initialize=0,
        bounds=(0, None),
        units=pyunits.megagallons,
        doc="Accumulate water produces in MG",
    )

    m.fs.pre_acc_prod = Var(
        initialize=0,
        bounds=(0, None),
        units=pyunits.megagallons,
        doc="Accumulate water produces in MG from previous step",
    )

    m.fs.acc_energy = Var(
        initialize=0,
        bounds=(0, None),
        units=pyunits.kWh,
        doc="Accumulate water produces in MG",
    )

    m.fs.pre_acc_energy = Var(
        initialize=0,
        bounds=(0, None),
        units=pyunits.kWh,
        doc="Accumulate water produces in MG from previous step",
    )

    m.fs.water_prod_100 = Var(
        initialize=14.8 / 24,
        bounds=(14.8 / 24, 14.8 / 24),
        units=pyunits.megagallons,
        doc="Water produced in a hour in MG",
    )

    m.fs.water_prod_50 = Var(
        initialize=14.8 / 24 / 2,
        bounds=(14.8 / 24 / 2, 14.8 / 24 / 2),
        units=pyunits.megagallons,
        doc="Water produced in a hour in MG",
    )

    m.fs.grid_cost = Var(
        initialize=0,
        bounds=(0, None),
        units=pyunits.dimensionless,
        doc="Grid cost for the month",
    )

    @m.Constraint(doc="Constraint to accumulate water production")
    def eq_acc_water_prod(b):
        return (
            b.fs.acc_prod
            == b.fs.pre_acc_prod
            + b.fs.water_prod_100 * b.fs.flow_100
            + b.fs.water_prod_50 * (b.fs.flow_50)
        )

    @m.Constraint(doc="Constraint to calculate total energy consumption")
    def eq_acc_energy(b):
        return (
            b.fs.acc_energy
            == b.fs.pre_acc_energy
            + (
                b.fs.water_prod_100 * b.fs.flow_100
                + b.fs.water_prod_50 * (b.fs.flow_50)
            )
            * b.fs.energy_per_MG
        )

    @m.Constraint(doc="Grid cost")
    def eq_grid_cost(b):
        return (
            b.fs.grid_cost
            == elec_price
            * (
                b.fs.water_prod_100 * b.fs.flow_100
                + b.fs.water_prod_50 * (b.fs.flow_50)
            )
            * b.fs.energy_per_MG
        )

    @m.Constraint(doc="AND gate")
    def eq_and(b):
        return b.fs.flow_50 * b.fs.flow_100 == 0

    @m.Constraint(doc="OR gate")
    def eq_or(b):
        return b.fs.flow_50 + b.fs.flow_100 == 1

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

    @mp.Expression(doc="Total cost")
    def total_cost(b):
        return sum([b.blocks[i].process.fs.grid_cost for i in range(n_time_points)])

    # @mp.Constraint(doc="Total production")
    # def total_production(b):
    #     return sum([b.blocks[i].process.fs.water_prod for i in range(n_time_points)]) >=10*pyunits.megagallons

    @mp.Constraint(doc="Total production")
    def total_production(b):
        return (
            sum(
                [
                    (
                        b.blocks[i].process.fs.water_prod_100
                        * b.blocks[i].process.fs.flow_100
                        + b.blocks[i].process.fs.water_prod_50
                        * (b.blocks[i].process.fs.flow_50)
                    )
                    for i in range(n_time_points)
                ]
            )
            >= 5 * pyunits.megagallons
        )

    # Set objective
    mp.obj = Objective(expr=mp.total_cost)

    return mp


if __name__ == "__main__":

    n_time_points = 25

    season = "winter"
    season = "summer"

    if season == "winter":
        elec_price = build_elec_price_winter()
    else:
        elec_price = build_elec_price_summer()
    print(elec_price)

    mp = create_wrd_mp(
        n_time_points=n_time_points,
        elec_price=elec_price,
    )
    # solver = get_solver()
    solver = SolverFactory("mindtpy")
    results = solver.solve(mp)

    flow_100 = [mp.blocks[i].process.fs.flow_100() for i in range(n_time_points)]
    flow_50 = [mp.blocks[i].process.fs.flow_50() for i in range(n_time_points)]
    prod = [
        mp.blocks[i].process.fs.water_prod_100() * mp.blocks[i].process.fs.flow_100()
        + mp.blocks[i].process.fs.water_prod_50() * (mp.blocks[i].process.fs.flow_50())
        for i in range(n_time_points)
    ]
    energy = [
        prod[i] * mp.blocks[i].process.fs.energy_per_MG() for i in range(n_time_points)
    ]

    # print(degrees_of_freedom(mp))
    time = np.linspace(0, n_time_points - 1, n_time_points)

    fig, ax = plt.subplots(figsize=(5, 5))
    ax.plot(time, flow_100, label="Binary 100%", color="blue")
    ax.plot(time, flow_50, label="Binary 50%", color="blue", linestyle="--")
    ax.plot(time, prod, label="Water Production (MG)", color="blue", marker="o")
    ax.set_xlim(0, 24)

    ax.set_ylim(0, 1.5)
    ax1 = ax.twinx()

    ax1.plot(time, elec_price, label="Electricity Price", color="black", marker="o")
    ax1.set_ylim(0, 0.5)

    ax2 = ax.twinx()
    ax2.plot(time, energy, label="Energy Consumption", color="orange", marker="o")
    ax2.spines["right"].set_position(("outward", 55))
    ax2.set_ylim(0, 100)

    if season == "summer":
        ax2.axvspan(0, 15, facecolor="lemonchiffon", alpha=0.3, label="Off Peak")
        ax2.axvspan(15, 21, facecolor="gold", alpha=0.3, label="On Peak")
        ax2.axvspan(21, 24, facecolor="lemonchiffon", alpha=0.3)
    elif season == "winter":
        ax2.axvspan(0, 8, facecolor="khaki", alpha=0.3, label="Off Peak")
        ax2.axvspan(8, 16, facecolor="lemonchiffon", alpha=0.3, label="Super Off Peak")
        ax2.axvspan(16, 21, facecolor="gold", alpha=0.3, label="Mid Peak")
        ax2.axvspan(21, 24, facecolor="khaki", alpha=0.3, label="Off Peak")

    ax.axhline(y=14.8 / 24, label="Maximum Production Capacity (MG per hour)")

    handle, label = ax.get_legend_handles_labels()
    handle1, label1 = ax1.get_legend_handles_labels()
    handle2, label2 = ax2.get_legend_handles_labels()

    handles = handle + handle1 + handle2
    labels = label + label1 + label2

    plt.legend(handles=handles, labels=labels, loc="upper left")

    ax.set_ylabel("Water production (MG)")
    ax1.set_ylabel("Electricity Price (2021 $/kWh)")
    ax2.set_ylabel("Energy Consumption (kWh)")

    # ax.set_xticklabels(month_names)
    ax.xaxis.set_major_locator(plt.MaxNLocator(25))
    ax.set_xlabel("Hours")

    plt.title(season)
    fig.tight_layout()
    plt.show()

    print("Total production in MG (Target 10):", mp.total_production())
    print("Total energy cost:", mp.total_cost())
