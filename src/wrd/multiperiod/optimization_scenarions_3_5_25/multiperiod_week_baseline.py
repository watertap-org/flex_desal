import numpy as np
import os
import pandas as pd
import logging
import matplotlib.pyplot as plt
from pathlib import Path

# Pyomo imports
from pyomo.environ import ConcreteModel, Var, Param, units as pyunits, Objective
from pyomo.util.check_units import assert_units_consistent
import matplotlib.dates as mdates

from idaes.core import FlowsheetBlock
from pyomo.environ import Var, Binary, Constraint, Objective, Expression, value
from watertap.core.util.model_diagnostics import *
from idaes.core.util.model_statistics import *

# IDAES imports
from idaes.apps.grid_integration.multiperiod.multiperiod import MultiPeriodModel
from idaes.core.solvers.get_solver import get_solver
import idaes.logger as idaeslog
from pyomo.environ import SolverFactory

if hasattr(pyunits, "USD_2021"):
    CURRENCY_UNIT = pyunits.USD_2021
elif hasattr(pyunits, "USD"):
    CURRENCY_UNIT = pyunits.USD
else:
    pyunits.load_definitions_from_strings(["USD = [currency]"])
    CURRENCY_UNIT = pyunits.USD

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


def build_elec_price_summer(n):
    # Delivery Pricing $/kWh
    on_peak_del = 0.01885
    mid_peak_del = 0.01766
    off_peak_del = 0.01741
    super_off_peak_del = 0

    # Generation Pricing $/kWh
    on_peak_gen = 0.13361
    mid_peak_gen = 0.12228  # MID PEAK ONLY OCCURS ON WEEKENDS
    off_peak_gen = 0.08419
    super_off_peak_gen = 0

    elec_price = np.ones(24)

    # off peak 12 AM - 4 PM
    elec_price[0:16] = off_peak_del + off_peak_gen
    # on peak 4 PM - 9 PM
    elec_price[16:21] = on_peak_del + on_peak_gen
    # off peak 9 PM - 12 AM
    elec_price[21:24] = off_peak_del + off_peak_gen

    # Repeat for the 7 days
    if value(n) > 24:
        elec_price = np.tile(elec_price, int(value(n) / 24))

    return elec_price


def build_elec_price_winter(n):
    # 2023
    # Delivery Pricing
    on_peak_del = 0
    mid_peak_del = 0.01927
    off_peak_del = 0.01811
    super_off_peak_del = 0.01745

    # Generation Pricing $/kWh
    on_peak_gen = 0
    mid_peak_gen = 0.09639
    off_peak_gen = 0.09695
    super_off_peak_gen = 0.05329

    elec_price = np.ones(24)

    # Off peak 12 AM - 8 AM
    elec_price[0:8] = off_peak_del + off_peak_gen
    # Super off peak 8 AM - 4 PM
    elec_price[8:16] = super_off_peak_del + super_off_peak_gen
    # Mid peak 4 PM - 9 PM
    elec_price[16:21] = mid_peak_del + mid_peak_gen
    # Off peak peak 9 PM - 12 AM
    elec_price[21:24] = off_peak_del + off_peak_gen

    if value(n) > 24:
        elec_price = np.tile(elec_price, int(value(n) / 24))
    return elec_price


def build_wrd_flowsheet(m=None, elec_price=0.1, peak_hours=list(range(16, 21))):

    if m is None:
        m = ConcreteModel()

    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.time_step = Param(
        initialize=1,
        mutable=True,
        units=pyunits.h,
        doc="Duration of each multiperiod time block",
    )

    m.fs.electricity_price = Param(
        initialize=elec_price,
        mutable=True,
        units=CURRENCY_UNIT / pyunits.kWh,
        doc="Electricity price for the current time block",
    )

    total_plant_production_capacity = (
        53150 / 24 * pyunits.m**3 / pyunits.h
    )  # m3 per hour
    train_production_capacity = (
        total_plant_production_capacity / 4
    )  # m3 per hour per train
    max_ro_train_power = (
        0.7421 * train_production_capacity / (pyunits.m**3 / pyunits.hr) - 169.73
    ) * pyunits.kW
    max_uf_power = (
        0.199 * total_plant_production_capacity / (pyunits.m**3 / pyunits.hr) - 30.0
    ) * pyunits.kW
    max_treatment_power = 4 * max_ro_train_power + max_uf_power

    m.fs.total_water_production = Var(
        initialize=total_plant_production_capacity,
        bounds=(0, total_plant_production_capacity),
        units=pyunits.m**3 / pyunits.h,
        doc="Water produced in a hour in m3",
    )

    m.fs.water_production_ro_train_1 = Var(
        initialize=train_production_capacity,
        bounds=(0, train_production_capacity),
        units=pyunits.m**3 / pyunits.h,
        doc="Volume of water treated by RO train 1",
    )

    m.fs.water_production_ro_train_2 = Var(
        initialize=train_production_capacity,
        bounds=(0, train_production_capacity),
        units=pyunits.m**3 / pyunits.h,
        doc="Volume of water treated by RO train 2",
    )

    m.fs.water_production_ro_train_3 = Var(
        initialize=train_production_capacity,
        bounds=(0, train_production_capacity),
        units=pyunits.m**3 / pyunits.h,
        doc="Volume of water treated by RO train 3",
    )

    m.fs.water_production_ro_train_4 = Var(
        initialize=train_production_capacity,
        bounds=(0, train_production_capacity),
        units=pyunits.m**3 / pyunits.h,
        doc="Volume of water treated by RO train 4",
    )

    # Create binary variables to indicate if train is on or off
    m.fs.train_1_on = Var(
        initialize=1,
        domain=Binary,
        doc="Binary variable indicating if RO train 1 is on",
    )

    m.fs.train_2_on = Var(
        initialize=1,
        domain=Binary,
        doc="Binary variable indicating if RO train 2 is on",
    )

    m.fs.train_3_on = Var(
        initialize=1,
        domain=Binary,
        doc="Binary variable indicating if RO train 3 is on",
    )

    m.fs.train_4_on = Var(
        initialize=1,
        domain=Binary,
        doc="Binary variable indicating if RO train 4 is on",
    )

    # Constraint connecting binary variable to the flowrate
    @m.Constraint(doc="Upper bound for flow depends on the binary variable")
    def eq_train_1_ub(b):
        return (
            b.fs.water_production_ro_train_1
            == train_production_capacity * b.fs.train_1_on
        )

    @m.Constraint(doc="Upper bound for flow depends on the binary variable")
    def eq_train_2_ub(b):
        return (
            b.fs.water_production_ro_train_2
            == train_production_capacity * b.fs.train_2_on
        )

    @m.Constraint(doc="Upper bound for flow depends on the binary variable")
    def eq_train_3_ub(b):
        return (
            b.fs.water_production_ro_train_3
            == train_production_capacity * b.fs.train_3_on
        )

    @m.Constraint(doc="Upper bound for flow depends on the binary variable")
    def eq_train_4_ub(b):
        return (
            b.fs.water_production_ro_train_4
            == train_production_capacity * b.fs.train_4_on
        )

    # Constraint to connect total water production to sum of RO train production
    @m.Constraint(doc="Total water production is sum of RO train production")
    def eq_total_water_production(b):
        return (
            b.fs.total_water_production
            == b.fs.water_production_ro_train_1
            + b.fs.water_production_ro_train_2
            + b.fs.water_production_ro_train_3
            + b.fs.water_production_ro_train_4
        )

    m.fs.treatment_energy_rate = Var(
        initialize=0,
        bounds=(0, max_treatment_power),
        units=pyunits.kWh / pyunits.h,
        doc="Total treatment energy required per hour",
    )

    def calculate_ro_power(flow, train_on):
        # Linear train power model (kW): 0 when train is off, fitted line when on
        return (
            0.7421 * train_production_capacity / (pyunits.m**3 / pyunits.hr) - 169.73
        ) * pyunits.kW

    def calculate_uf_power(flow, uf_on):
        # Linear UF power model (kW)
        # If train_1_on == 0, then the whole system is off and no power use from UF
        return (0.199 * flow / (pyunits.m**3 / pyunits.hr) - 30.0 * uf_on) * pyunits.kW

    @m.Constraint(doc="Calculate total treatment energy rate")
    def eq_treatment_energy_rate(b):
        return b.fs.treatment_energy_rate == (
            calculate_ro_power(b.fs.water_production_ro_train_1, b.fs.train_1_on)
            + calculate_ro_power(b.fs.water_production_ro_train_2, b.fs.train_2_on)
            + calculate_ro_power(b.fs.water_production_ro_train_3, b.fs.train_3_on)
            + calculate_ro_power(b.fs.water_production_ro_train_4, b.fs.train_4_on)
            + calculate_uf_power(b.fs.total_water_production, b.fs.train_1_on)
        )

    m.fs.acc_production = Var(
        initialize=0,
        bounds=(0, None),
        units=pyunits.m**3,
        doc="Accumulate water produces in m3",
    )

    m.fs.pre_acc_production = Var(
        initialize=0,
        bounds=(0, None),
        units=pyunits.m**3,
        doc="Accumulate water produced in m3 from previous step",
    )

    m.fs.acc_energy = Var(
        initialize=0,
        bounds=(0, None),
        units=pyunits.kWh,
        doc="Accumulate energy consumption in kWh",
    )

    m.fs.pre_acc_energy = Var(
        initialize=0,
        bounds=(0, None),
        units=pyunits.kWh,
        doc="Accumulate energy consumption in kWh from previous step",
    )

    m.fs.grid_cost = Var(
        initialize=0,
        bounds=(0, None),
        units=CURRENCY_UNIT,
        doc="Electricity cost for each time step",
    )

    @m.Constraint(doc="Constraint to accumulate water production")
    def eq_acc_water_prod(b):
        return (
            b.fs.acc_production
            == b.fs.pre_acc_production + b.fs.total_water_production * b.fs.time_step
        )

    @m.Constraint(doc="Constraint to calculate total energy consumption")
    def eq_acc_energy(b):
        return (
            b.fs.acc_energy
            == b.fs.pre_acc_energy + b.fs.treatment_energy_rate * b.fs.time_step
        )

    @m.Constraint(doc="Grid cost")
    def eq_grid_cost(b):
        return (
            b.fs.grid_cost
            == b.fs.electricity_price * b.fs.treatment_energy_rate * b.fs.time_step
        )

    return m


def get_wrd_variable_pairs(t1, t2):
    # Connect the accumulated water produced
    return [
        (t1.fs.acc_production, t2.fs.pre_acc_production),
        (t1.fs.acc_energy, t2.fs.pre_acc_energy),
    ]


def unfix_dof(m):
    # Train 1 and 2 are always on, so we only vary the fraction of water treated by train 3 and 4
    m.fs.water_production_ro_train_3.unfix()
    m.fs.water_production_ro_train_4.unfix()
    m.fs.train_3_on.unfix()
    m.fs.train_4_on.unfix()
    m.fs.water_production_ro_train_1.unfix()
    m.fs.water_production_ro_train_2.unfix()
    m.fs.train_1_on.unfix()
    m.fs.train_2_on.unfix()
    return None


def initialize_mp(m):
    print("Initializing multi-period model...")
    # Check if first time step
    max_train_flow = 53150 / 24 / 4  # m3/hr
    m.fs.water_production_ro_train_1.fix(max_train_flow)
    m.fs.water_production_ro_train_2.fix(max_train_flow)
    m.fs.water_production_ro_train_3.fix(max_train_flow * 1)
    m.fs.water_production_ro_train_4.fix(max_train_flow * 1)

    m.fs.train_1_on.fix(1)
    m.fs.train_2_on.fix(1)
    m.fs.train_3_on.fix(1)
    m.fs.train_4_on.fix(1)


def create_wrd_mp(
    n_days=1,
    n_time_points=24,
    elec_price=elec_price,
    daily_production_target=12 * pyunits.m**3 / pyunits.day,
    total_water_production_target=12 * pyunits.m**3 / pyunits.day,
    peak_hours=list(range(16, 21)),  # 4 PM to 9 PM
    demand_charges={"fixed_demand_price": 5, "on_peak_demand_price": 15},
):
    """
    This function creates a multi-period flowsheet object for each month for the WRD plant. This object contains
    a pyomo model with a block for each time instance.

    Args:
        n_time_points: Number of time blocks to create

    Returns:
        Object containing multi-period vagmd batch flowsheet model
    """
    m = ConcreteModel()

    m.fs = FlowsheetBlock(dynamic=False)

    total_plant_production_capacity = 53150 / 24 * pyunits.m**3 / pyunits.h
    train_production_capacity = total_plant_production_capacity / 4
    max_ro_train_power = (
        0.7421 * train_production_capacity / (pyunits.m**3 / pyunits.hr) - 169.73
    ) * pyunits.kW
    max_uf_power = (
        0.199 * total_plant_production_capacity / (pyunits.m**3 / pyunits.hr) - 30.0
    ) * pyunits.kW
    max_treatment_power = 4 * max_ro_train_power + max_uf_power

    m.fs.mp = MultiPeriodModel(
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
            # "peak_hours": list(range(16, 21)),  # 4 PM to 9 PM # Not sure that's the best way to implement
        }
        for t in range(n_time_points)
    }

    unfix_dof_options = {t: {} for t in range(n_time_points)}

    m.fs.mp.build_multi_period_model(
        model_data_kwargs=flowsheet_options,
        flowsheet_options=flowsheet_options,
        initialization_options=None,
        unfix_dof_options=None,
    )

    for t in range(n_time_points):
        initialize_mp(m.fs.mp.blocks[t].process)
        unfix_dof(m.fs.mp.blocks[t].process)

    m.fs.mp.blocks[0].process.fs.pre_acc_production.fix(0)
    m.fs.mp.blocks[0].process.fs.pre_acc_energy.fix(0)

    split_idx = int(0.75 * n_time_points)
    time_points_75 = range(split_idx)
    time_points_25 = range(split_idx, n_time_points)

    @m.Constraint(
        time_points_75, doc="Force all four RO trains ON for first 75% of time periods"
    )
    def eq_all_trains_on_first_75_percent(b, i):
        return (
            b.fs.mp.blocks[i].process.fs.train_1_on
            + b.fs.mp.blocks[i].process.fs.train_2_on
            + b.fs.mp.blocks[i].process.fs.train_3_on
            + b.fs.mp.blocks[i].process.fs.train_4_on
            == 4
        )

    @m.Constraint(
        time_points_25, doc="Force all four RO trains OFF for last 25% of time periods"
    )
    def eq_all_trains_off_last_25_percent(b, i):
        return (
            b.fs.mp.blocks[i].process.fs.train_1_on
            + b.fs.mp.blocks[i].process.fs.train_2_on
            + b.fs.mp.blocks[i].process.fs.train_3_on
            + b.fs.mp.blocks[i].process.fs.train_4_on
            == 0
        )

    # Take values as inputs?
    m.fs.fixed_demand_price = Param(
        initialize=demand_charges["fixed_demand_price"],
        mutable=True,
        units=CURRENCY_UNIT / pyunits.kW,
        doc="Demand charge associate with time of peak power",
    )

    m.fs.on_peak_demand_price = Param(
        initialize=demand_charges["on_peak_demand_price"],
        mutable=True,
        units=CURRENCY_UNIT / pyunits.kW,
        doc="Demand charge associated with greatest power value within the on-peak hours",
    )

    m.fs.highest_demand = Var(
        initialize=0,
        bounds=(0, max_treatment_power),
        units=pyunits.kW,
        doc="Demand during highest period of energy use",
    )

    m.fs.highest_on_peak_demand = Var(
        initialize=0,
        bounds=(0, max_treatment_power),
        units=pyunits.kW,
        doc="Demand during highest period of energy use within the peak hours",
    )

    m.fs.fixed_demand_charge = Expression(
        expr=m.fs.fixed_demand_price * m.fs.highest_demand,
        doc="Demand charge for highest period of energy use",
    )

    m.fs.on_peak_demand_charge = Expression(
        expr=m.fs.on_peak_demand_price * m.fs.highest_on_peak_demand,
        doc="Demand charge for highest period of energy use within the peak hours",
    )

    @m.Constraint(
        range(n_time_points),
        doc="Upper bound highest demand by each period energy rate",
    )
    def eq_highest_demand(b, i):
        # this could also be formulated using max(), but idk how that interacts with the solver
        return b.fs.highest_demand >= b.fs.mp.blocks[i].process.fs.treatment_energy_rate

    m.fs.peak_hours = [h for h in range(n_time_points) if (h % 24) in peak_hours]

    @m.Constraint(
        m.fs.peak_hours,
        doc="Upper bound on highest on peak demand by each period energy rate during peak hours",
    )
    def eq_highest_on_peak_demand(b, i):
        return (
            b.fs.highest_on_peak_demand
            >= b.fs.mp.blocks[i].process.fs.treatment_energy_rate
        )

    @m.Expression(doc="Total cost")
    def total_cost(b):
        return (
            (30 / n_days)
            * sum(
                [b.fs.mp.blocks[i].process.fs.grid_cost for i in range(n_time_points)]
            )
            + b.fs.fixed_demand_charge
            + b.fs.on_peak_demand_charge
        )

    @m.Constraint(doc="Total production")
    def total_production(b):
        return sum(
            [
                b.fs.mp.blocks[i].process.fs.total_water_production
                * b.fs.mp.blocks[i].process.fs.time_step
                for i in range(n_time_points)
            ]
        ) >= pyunits.convert(total_water_production_target, to_units=pyunits.m**3)

    # Set objective
    m.fs.obj = Objective(expr=m.total_cost)

    return m


def plot_function(m, n_time_points, season):

    time = np.linspace(0, n_time_points - 1, n_time_points)

    fig, (ax, ax_trains) = plt.subplots(2, 1, figsize=(10, 8))

    # First subplot: Total production, electricity price, and energy
    ax.plot(
        time + 0.5, prod, label="Water Production (m3/hr)", color="blue", marker="o"
    )
    ax.set_ylim(0, 3000)
    ax1 = ax.twinx()

    ax1.plot(
        time + 0.5, elec_price, label="Electricity Price", color="black", marker="o"
    )
    ax1.set_ylim(0, 0.4)

    ax2 = ax.twinx()
    ax2.plot(time + 0.5, energy, label="Energy Consumption", color="orange", marker="o")
    ax2.spines["right"].set_position(("outward", 55))
    ax2.set_ylim(0, 2500)

    if season == "summer":
        for i in range(int(n_time_points / 24)):
            ax2.axvspan(
                24 * i,
                24 * i + 16,
                facecolor="lemonchiffon",
                alpha=0.3,
                label="Off Peak" if i == 0 else "_nolegend_",
                zorder=0,
            )
            ax2.axvspan(
                24 * i + 16,
                24 * i + 21,
                facecolor="gold",
                alpha=0.3,
                label="On Peak" if i == 0 else "_nolegend_",
                zorder=0,
            )
            ax2.axvspan(
                24 * i + 21,
                24 * i + 24,
                facecolor="lemonchiffon",
                alpha=0.3,
                label="_nolegend_",
                zorder=0,
            )
    elif season == "winter":
        for i in range(int(n_time_points / 24)):
            ax2.axvspan(
                24 * i,
                8 + 24 * i,
                facecolor="khaki",
                alpha=0.3,
                label="Off Peak" if i == 0 else "_nolegend_",
                zorder=0,
            )
            ax2.axvspan(
                8 + 24 * i,
                16 + 24 * i,
                facecolor="lemonchiffon",
                alpha=0.3,
                label="Super Off Peak" if i == 0 else "_nolegend_",
                zorder=0,
            )
            ax2.axvspan(
                16 + 24 * i,
                21 + 24 * i,
                facecolor="gold",
                alpha=0.3,
                label="Mid Peak" if i == 0 else "_nolegend_",
                zorder=0,
            )
            ax2.axvspan(
                21 + 24 * i,
                24 + 24 * i,
                facecolor="khaki",
                alpha=0.3,
                label="_nolegend_",
                zorder=0,
            )

    ax.axhline(y=53150 / 24, label="Maximum Production Capacity (m3/h)")

    handle, label = ax.get_legend_handles_labels()
    handle1, label1 = ax1.get_legend_handles_labels()
    handle2, label2 = ax2.get_legend_handles_labels()

    handles = handle + handle1 + handle2
    labels = label + label1 + label2

    leg = ax2.legend(
        handles=handles,
        labels=labels,
        loc="upper left",
        framealpha=1.0,
        ncol=2,
        fontsize=11,
    )
    leg.set_zorder(1000)
    leg.get_frame().set_facecolor("white")  # optional, keeps it clean
    ax.set_ylabel("Water production (m3/h)", fontsize=12)
    ax1.set_ylabel("Electricity Price (2021 $/kWh)", fontsize=11)
    ax2.set_ylabel("Energy Consumption (kWh)", fontsize=11)
    ax.xaxis.set_major_locator(plt.MaxNLocator(24))
    ax.set_xlabel("Hours", fontsize=12)
    ax.set_title(season + ": Plant Shutdown Scenario", fontsize=14, fontweight="bold")
    # Tick labels (all axes)
    for a in (ax, ax1, ax2, ax_trains):
        a.tick_params(axis="both", labelsize=11)

    # Extract RO train flow rates and convert to % of max flow
    max_train_flow = 53150 / 24 / 4  # m3/hr
    train_1_flows = [
        m.fs.mp.blocks[i].process.fs.water_production_ro_train_1()
        * m.fs.mp.blocks[i].process.fs.train_1_on()
        / max_train_flow
        * 100
        for i in range(n_time_points)
    ]
    train_2_flows = [
        m.fs.mp.blocks[i].process.fs.water_production_ro_train_2()
        * m.fs.mp.blocks[i].process.fs.train_2_on()
        / max_train_flow
        * 100
        for i in range(n_time_points)
    ]
    train_3_flows = [
        m.fs.mp.blocks[i].process.fs.water_production_ro_train_3()
        * m.fs.mp.blocks[i].process.fs.train_3_on()
        / max_train_flow
        * 100
        for i in range(n_time_points)
    ]
    train_4_flows = [
        m.fs.mp.blocks[i].process.fs.water_production_ro_train_4()
        * m.fs.mp.blocks[i].process.fs.train_4_on()
        / max_train_flow
        * 100
        for i in range(n_time_points)
    ]
    train_flows = [train_1_flows, train_2_flows, train_3_flows, train_4_flows]

    # Second subplot: RO train flow rates as % of max flow
    ax_trains.plot(
        time + 0.5, train_flows[0], label="RO Train 1", marker="o", linewidth=2
    )
    ax_trains.plot(
        time + 0.5, train_flows[1], label="RO Train 2", marker="s", linewidth=2
    )
    ax_trains.plot(
        time + 0.5, train_flows[2], label="RO Train 3", marker="^", linewidth=2
    )
    ax_trains.plot(
        time + 0.5, train_flows[3], label="RO Train 4", marker="d", linewidth=2
    )

    ax_trains.set_ylabel("Flow Rate (% of Max)", fontsize=12)
    ax_trains.set_xlabel("Hours", fontsize=12)
    ax_trains.set_title("RO Train Flow Rates as % of Maximum", fontsize=14)
    ax_trains.set_ylim(0, 110)
    ax_trains.axhline(
        y=100,
        color="red",
        linestyle="--",
        linewidth=1,
        alpha=0.5,
        label="Max Capacity",
        zorder=0,
    )
    ax_trains.legend(loc="lower left", fontsize=11)
    ax_trains.grid(True, alpha=0.3)
    ax_trains.xaxis.set_major_locator(plt.MaxNLocator(24))
    ax_trains.tick_params(axis="both", labelsize=11)
    fig.tight_layout()
    output_path = Path(__file__).resolve().parent / f"{Path(__file__).stem}_wrd_{season}_plant_shutdown.png"
    fig.savefig(output_path, dpi=600)
    plt.show()


def plot_grid_cost_over_time(m, n_time_points, season=None):
    time = np.arange(n_time_points)
    grid_cost = [
        value(
            pyunits.convert(
                m.fs.mp.blocks[i].process.fs.grid_cost,
                to_units=CURRENCY_UNIT,
            )
        )
        for i in range(n_time_points)
    ]

    fig, ax = plt.subplots(figsize=(10, 4))
    ax.plot(time + 0.5, grid_cost, marker="o", linewidth=1.5, color="tab:green")
    ax.set_xlabel("Hours")
    ax.set_ylabel("Grid Cost (2021 $/hr)")
    title = "Grid Cost Over Time"
    if season is not None:
        title = f"{title} - {season.capitalize()}"
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    ax.xaxis.set_major_locator(plt.MaxNLocator(24))
    fig.tight_layout()
    plt.show()


def print_unfixed_vars(model):
    print("Unfixed variables contributing to degrees of freedom:")
    for v in model.component_data_objects(ctype=Var, descend_into=True):
        if not v.fixed:
            print(f"  {v.name}")


if __name__ == "__main__":
    n_days = 7
    n_time_points = 24 * n_days
    daily_production_target = 0 * pyunits.m**3 / pyunits.day
    total_water_production_target = (
        0.74 * 53150 * pyunits.m**3 / pyunits.day * n_days * pyunits.day
    )  # 74 to give a bit of wiggle room

    season = "summer"

    if season == "winter":
        elec_price = build_elec_price_winter(n=n_time_points)
        peak_hours = list(range(16, 21))
        demand_charges = {
            "fixed_demand_price": 19.62,
            "on_peak_demand_price": 7.99 + 2.55,
        }  # March 2023
    else:
        elec_price = build_elec_price_summer(n=n_time_points)
        peak_hours = list(range(16, 21))  # But only the weekdays actually!
        demand_charges = {
            "fixed_demand_price": 19.94,
            "on_peak_demand_price": 22.10 + 14.68,
        }  # June 2023

    m = create_wrd_mp(
        n_days=n_days,
        n_time_points=n_time_points,
        elec_price=elec_price,
        daily_production_target=daily_production_target,
        total_water_production_target=total_water_production_target,
        peak_hours=peak_hours,
        demand_charges=demand_charges,
    )
    assert_units_consistent(m)
    # print_unfixed_vars(m)

    # solver = get_solver()
    # solver = SolverFactory('mindtpy')
    # results = solver.solve(m)
    os.environ["PATH"] = (
        r"C:\Users\rchurchi\AppData\Local\anaconda3\pkgs\glpk-4.65-h17947e8_4\Library\bin"
        + os.pathsep
        + os.environ.get("PATH", "")
    )

    # dt = DiagnosticsToolbox(m)

    # solver = SolverFactory("mindtpy")
    # results = solver.solve(
    #     m,
    #     strategy="OA",
    #     mip_solver="glpk",
    #     nlp_solver="ipopt",
    #     tee=True,
    # )

    solver = SolverFactory("glpk")
    # Practical stopping criteria for long-horizon MILPs
    solver.options["mipgap"] = 0.005
    solver.options["tmlim"] = 180
    results = solver.solve(m, tee=True)

    prod = [
        m.fs.mp.blocks[i].process.fs.total_water_production()
        for i in range(n_time_points)
    ]
    energy = [
        value(
            pyunits.convert(
                m.fs.mp.blocks[i].process.fs.treatment_energy_rate,
                to_units=pyunits.kWh / pyunits.h,
            )
        )
        for i in range(n_time_points)
    ]

    print(degrees_of_freedom(m))

    print("-" * 10, "Ouputs over Period", "-" * 10)
    print("Total production in m3:", m.total_production())
    print("Total target water production in m3:", total_water_production_target())
    print(
        "Total energy consumption in kWh:",
        m.fs.mp.blocks[n_time_points - 1].process.fs.acc_energy(),
    )

    print("-" * 10, "Monthly Costs", "-" * 10)
    print("Fixed demand charge:", m.fs.fixed_demand_charge(), "2021 $")
    print("On-peak demand charge:", m.fs.on_peak_demand_charge(), "2021 $")
    print(
        "Consumption charge:",
        (28 / n_days)
        * sum(m.fs.mp.blocks[i].process.fs.grid_cost() for i in range(n_time_points)),
        "2021 $",
    )
    print("Total electricity cost for month:", m.total_cost(), "2021 $")

    plot_function(m, n_time_points, season)
    # plot_grid_cost_over_time(m, n_time_points, season)
