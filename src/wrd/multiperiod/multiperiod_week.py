import numpy as np
import os
import pandas as pd
import logging
import matplotlib.pyplot as plt

# Pyomo imports
from pyomo.environ import ConcreteModel, Var, Param, units as pyunits, Objective
from pyomo.util.check_units import assert_units_consistent
import matplotlib.dates as mdates

from idaes.core import FlowsheetBlock
from pyomo.environ import Var, Binary, Constraint, Objective, value
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
    on_peak_del = 0.01971
    mid_peak_del = 0.01971
    off_peak_del = 0.01971
    super_off_peak_del = 0

    # Generation Pricing $/kWh
    on_peak_gen = 0.09934
    mid_peak_gen = 0.0891
    off_peak_gen = 0.05782
    super_off_peak_gen = 0

    elec_price = np.ones(24)

    # off peak 12 AM - 4 PM
    elec_price[0:16] = off_peak_del + off_peak_gen
    # on peak 4 PM - 9 PM
    elec_price[16:21] = on_peak_del + mid_peak_gen
    # off peak 9 PM - 12 AM
    elec_price[21:24] = off_peak_del + off_peak_gen

    # Repeat for the 7 days
    if value(n) > 24:
        elec_price = np.tile(elec_price, int(value(n) / 24))

    return elec_price


def build_elec_price_winter(n):
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


def build_wrd_flowsheet(
    m=None,
    elec_price=0.1,
):

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

    total_plant_production_capacity = 53150 / 24  # m3 per hour
    train_production_capacity = (
        total_plant_production_capacity / 4
    )  # m3 per hour per train

    m.fs.total_water_production = Var(
        initialize=total_plant_production_capacity,
        bounds=(total_plant_production_capacity / 2, total_plant_production_capacity),
        units=pyunits.m**3 / pyunits.h,
        doc="Water produced in a hour in m3",
    )

    m.fs.water_production_ro_train_1 = Var(
        initialize=train_production_capacity,
        bounds=(train_production_capacity * 0.85, train_production_capacity),
        units=pyunits.m**3 / pyunits.h,
        doc="Volume of water treated by RO train 1",
    )

    m.fs.water_production_ro_train_2 = Var(
        initialize=train_production_capacity,
        bounds=(train_production_capacity * 0.85, train_production_capacity),
        units=pyunits.m**3 / pyunits.h,
        doc="Volume of water treated by RO train 2",
    )

    m.fs.water_production_ro_train_3 = Var(
        initialize=train_production_capacity,
        bounds=(train_production_capacity * 0.85, train_production_capacity),
        units=pyunits.m**3 / pyunits.h,
        doc="Volume of water treated by RO train 3",
    )

    m.fs.water_production_ro_train_4 = Var(
        initialize=train_production_capacity,
        bounds=(train_production_capacity * 0.85, train_production_capacity),
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

    # Constraint to connect total water production to sum of RO train production
    @m.Constraint(doc="Total water production is sum of RO train production")
    def eq_total_water_production(b):
        return (
            b.fs.total_water_production
            == b.fs.water_production_ro_train_1 * b.fs.train_1_on
            + b.fs.water_production_ro_train_2 * b.fs.train_2_on
            + b.fs.water_production_ro_train_3 * b.fs.train_3_on
            + b.fs.water_production_ro_train_4 * b.fs.train_4_on
        )

    # Function to calculate energy consumption per m3 of water treated.
    def calculate_energy_intensity(flow):
        return (
            (5e-06 * flow / (pyunits.m**3 / pyunits.hr) + 1.1)
            * pyunits.kWh
            / pyunits.m**3
        )
        # return (2E-06*flow**2 - 0.0026*flow + 1.1007)

    m.fs.treatment_energy_per_m3 = Var(
        initialize=0.48,
        units=pyunits.kWh / pyunits.m**3,
        doc="Total treatment energy required per m3",
    )

    # Constraint to calculate treatment energy per m3 based on flow through each train
    @m.Constraint(doc="Calculate treatment energy per m3")
    def eq_treatment_energy_per_m3(b):
        return (
            b.fs.treatment_energy_per_m3
            == (
                calculate_energy_intensity(b.fs.water_production_ro_train_1)
                * b.fs.water_production_ro_train_1
                * b.fs.train_1_on
                + calculate_energy_intensity(b.fs.water_production_ro_train_2)
                * b.fs.water_production_ro_train_2
                * b.fs.train_2_on
                + calculate_energy_intensity(b.fs.water_production_ro_train_3)
                * b.fs.water_production_ro_train_3
                * b.fs.train_3_on
                + calculate_energy_intensity(b.fs.water_production_ro_train_4)
                * b.fs.water_production_ro_train_4
                * b.fs.train_4_on
            )
            / b.fs.total_water_production
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
            == b.fs.pre_acc_energy
            + b.fs.total_water_production
            * b.fs.treatment_energy_per_m3
            * b.fs.time_step
        )

    @m.Constraint(doc="Grid cost")
    def eq_grid_cost(b):
        return (
            b.fs.grid_cost
            == b.fs.electricity_price
            * b.fs.total_water_production
            * b.fs.treatment_energy_per_m3
            * b.fs.time_step
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
    return None


def initialize_mp(m):
    print("Initializing multi-period model...")
    # Check if first time step
    max_train_flow = 53150 / 24 / 4  # m3/hr
    m.fs.water_production_ro_train_1.fix(max_train_flow)
    m.fs.water_production_ro_train_2.fix(max_train_flow)
    m.fs.water_production_ro_train_3.fix(max_train_flow * 0.85)
    m.fs.water_production_ro_train_4.fix(max_train_flow * 0.85)

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

    # Constraints to connect non-working hours water production to be constant
    # Non working hours
    non_working_hours_morning = [0, 1, 2, 3, 4, 5, 6, 7, 8]
    non_working_hours_evening = [17, 18, 19, 20, 21, 22, 23]

    m.fs.non_working_hours_morning = [
        h for h in range(n_time_points) if (h % 24) in non_working_hours_morning
    ]
    m.fs.non_working_hours_evening = [
        h for h in range(n_time_points) if (h % 24) in non_working_hours_evening
    ]

    # Constraints to keep production constant during non-working hours (skip at transitions)
    @m.Constraint(
        m.fs.non_working_hours_morning,
        doc="Production constant during non-working hours",
    )
    def eq_non_working_hrs_morning(b, h):
        if h == 0:  # or (h % 24)==0:
            return Constraint.Skip
        return (
            b.fs.mp.blocks[h].process.fs.total_water_production
            == b.fs.mp.blocks[h - 1].process.fs.total_water_production
        )

    @m.Constraint(
        m.fs.non_working_hours_evening,
        doc="Production constant during non-working hours",
    )
    def eq_non_working_hrs_evening(b, h):
        # if h == 16 or (h % 24)==16:
        #     return Constraint.Skip
        return (
            b.fs.mp.blocks[h].process.fs.total_water_production
            == b.fs.mp.blocks[h - 1].process.fs.total_water_production
        )

    # # Constraints to keep trains on/off state constant during non-working hours
    # @m.Constraint(non_working_hours, doc="Train 1 on/off state constant during non-working hours")
    # def eq_train_1_on_nwh(b, h):
    #     if h == 0 or (h % 24) in [9, 16]:
    #         return Constraint.Skip
    #     return b.fs.mp.blocks[h].process.fs.train_1_on == b.fs.mp.blocks[h - 1].process.fs.train_1_on

    # @m.Constraint(non_working_hours, doc="Train 2 on/off state constant during non-working hours")
    # def eq_train_2_on_nwh(b, h):
    #     if h == 0 or (h % 24) in [9, 16]:
    #         return Constraint.Skip
    #     return b.fs.mp.blocks[h].process.fs.train_2_on == b.fs.mp.blocks[h - 1].process.fs.train_2_on

    @m.Constraint(
        m.fs.non_working_hours_morning,
        doc="Train 3 on/off state constant during non-working hours",
    )
    def eq_train_3_on_nwh_morning(b, h):
        if h == 0:  # This should ensure the midnight constraint
            return Constraint.Skip
        return (
            b.fs.mp.blocks[h].process.fs.train_3_on
            == b.fs.mp.blocks[h - 1].process.fs.train_3_on
        )

    @m.Constraint(
        m.fs.non_working_hours_evening,
        doc="Train 3 on/off state constant during non-working hours",
    )
    def eq_train_3_on_nwh_evening(b, h):
        if h == 0:
            return Constraint.Skip
        return (
            b.fs.mp.blocks[h].process.fs.train_3_on
            == b.fs.mp.blocks[h - 1].process.fs.train_3_on
        )

    # Not sure if the below is really needed due to the eq_train_3_over_train_4 constraint.
    # @m.Constraint(m.fs.non_working_hours_morning, doc="Train 4 on/off state constant during non-working hours")
    # def eq_train_4_on_nwh_morning(b, h):
    #     if h == 0: #or (h % 24)==0: --> this causes issues but idk why
    #         return Constraint.Skip
    #     return b.fs.mp.blocks[h].process.fs.train_4_on == b.fs.mp.blocks[h - 1].process.fs.train_4_on

    # @m.Constraint(m.fs.non_working_hours_evening, doc="Train 4 on/off state constant during non-working hours")
    # def eq_train_4_on_nwh_evening(b, h):
    #     if h == 0:
    #         return Constraint.Skip
    #     return b.fs.mp.blocks[h].process.fs.train_4_on == b.fs.mp.blocks[h - 1].process.fs.train_4_on

    # Ensure train 4 turns off first before train 3
    @m.Constraint(
        range(n_time_points),
        doc="Train 3 will always have the higher production than train 4",
    )
    def eq_train_3_over_train_4(b, i):
        return (
            b.fs.mp.blocks[i].process.fs.water_production_ro_train_3
            >= b.fs.mp.blocks[i].process.fs.water_production_ro_train_4
        )

    @m.Constraint(range(n_time_points), doc="Train 3 must be on for train 4 to be on")
    def eq_train_3_on_train_4_on(b, i):
        return (
            b.fs.mp.blocks[i].process.fs.train_3_on
            >= b.fs.mp.blocks[i].process.fs.train_4_on
        )

    # # Constraint disallowing off state for only one hour
    # @m.Constraint(range(n_time_points), doc="Trains cannot be on for only one hour")
    # def eq_no_one_hour_on_train3(b, i):
    #     if i < 2:
    #         return Constraint.Skip
    #     return b.fs.mp.blocks[i].process.fs.train_3_on * b.fs.mp.blocks[i - 1].process.fs.train_3_on >= b.fs.mp.blocks[i].process.fs.train_3_on - b.fs.mp.blocks[i-2].process.fs.train_3_on

    # @m.Constraint(range(n_time_points), doc="Trains cannot be on for only one hour")
    # def eq_no_one_hour_on_train4(b, i):
    #     if i < 2:
    #         return Constraint.Skip
    #     return b.fs.mp.blocks[i].process.fs.train_4_on * b.fs.mp.blocks[i - 1].process.fs.train_4_on >= b.fs.mp.blocks[i].process.fs.train_4_on - b.fs.mp.blocks[i-2].process.fs.train_4_on

    @m.Constraint(
        range(1, n_time_points - 1), doc="Disallow on-off-on pattern for train 3"
    )
    def eq_no_on_off_on_train3(b, i):
        return (
            b.fs.mp.blocks[i - 1].process.fs.train_3_on
            + b.fs.mp.blocks[i + 1].process.fs.train_3_on
            - b.fs.mp.blocks[i].process.fs.train_3_on
            <= 1
        )

    @m.Constraint(
        range(1, n_time_points - 1), doc="Disallow on-off-on pattern for train 4"
    )
    def eq_no_on_off_on_train4(b, i):
        return (
            b.fs.mp.blocks[i - 1].process.fs.train_4_on
            + b.fs.mp.blocks[i + 1].process.fs.train_4_on
            - b.fs.mp.blocks[i].process.fs.train_4_on
            <= 1
        )

    # @m.Constraint(range(n_time_points), doc="Production should not change at midnight")
    # def eq_midnight_constraint_on_off(b, h):
    #     if h == 0 or (h % 24) != 0:
    #         return Constraint.Skip
    #     return (
    #         b.fs.mp.blocks[h].process.fs.train_3_on
    #         == b.fs.mp.blocks[h - 1].process.fs.train_3_on
    #     )

    @m.Expression(doc="Total cost")
    def total_cost(b):
        return sum(
            [b.fs.mp.blocks[i].process.fs.grid_cost for i in range(n_time_points)]
        )

    # Add constraint for every 24 hours of water production
    m.fs.days = range(int(value(n_days)))

    # @m.Constraint(m.fs.days, doc="Daily production target for each 24-hour period")
    # def daily_production_constraint(b, day):
    #     start_idx = day * 24
    #     end_idx = min((day + 1) * 24, n_time_points)
    #     return (
    #         sum([b.fs.mp.blocks[i].process.fs.total_water_production for i in range(start_idx, end_idx)])
    #         >= pyunits.convert(daily_production_target * pyunits.days, to_units=pyunits.m**3)
    #     )

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
    ax.plot(time + 0.5, prod, label="Water Production (MG)", color="blue", marker="o")
    ax.set_ylim(0, 1)
    ax1 = ax.twinx()

    ax1.plot(
        time + 0.5, elec_price, label="Electricity Price", color="black", marker="o"
    )
    ax1.set_ylim(0, 0.5)

    ax2 = ax.twinx()
    ax2.plot(time + 0.5, energy, label="Energy Consumption", color="orange", marker="o")
    ax2.spines["right"].set_position(("outward", 55))
    # ax2.set_ylim(0, 2000)

    if season == "summer":
        ax2.axvspan(0, 16, facecolor="lemonchiffon", alpha=0.3, label="Off Peak")
        ax2.axvspan(16, 21, facecolor="gold", alpha=0.3, label="On Peak")
        ax2.axvspan(21, 24, facecolor="lemonchiffon", alpha=0.3)
    elif season == "winter":
        ax2.axvspan(0, 8, facecolor="khaki", alpha=0.3, label="Off Peak")
        ax2.axvspan(8, 16, facecolor="lemonchiffon", alpha=0.3, label="Super Off Peak")
        ax2.axvspan(16, 21, facecolor="gold", alpha=0.3, label="Mid Peak")
        ax2.axvspan(21, 24, facecolor="khaki", alpha=0.3)

    ax.axhline(y=53150 / 24, label="Maximum Production Capacity (m3/h)")

    handle, label = ax.get_legend_handles_labels()
    handle1, label1 = ax1.get_legend_handles_labels()
    handle2, label2 = ax2.get_legend_handles_labels()

    handles = handle + handle1 + handle2
    labels = label + label1 + label2

    ax.legend(handles=handles, labels=labels, loc="upper left")

    ax.set_ylabel("Water production (m3/h)")
    ax1.set_ylabel("Electricity Price (2021 $/kWh)")
    ax2.set_ylabel("Energy Consumption (kWh)")
    ax.xaxis.set_major_locator(plt.MaxNLocator(24))
    ax.set_xlabel("Hours")
    ax.set_title(season)

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

    ax_trains.set_ylabel("Flow Rate (% of Max)")
    ax_trains.set_xlabel("Hours")
    ax_trains.set_title("RO Train Flow Rates as % of Maximum")
    ax_trains.set_ylim(0, 110)
    ax_trains.axhline(
        y=100, color="red", linestyle="--", linewidth=1, alpha=0.5, label="Max Capacity"
    )
    ax_trains.legend(loc="best")
    ax_trains.grid(True, alpha=0.3)
    ax_trains.xaxis.set_major_locator(plt.MaxNLocator(24))
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
        0.75 * 53150 * pyunits.m**3 / pyunits.day * n_days * pyunits.day
    )

    season = "winter"
    # season = 'summer'

    if season == "winter":
        elec_price = build_elec_price_winter(n=n_time_points)
    else:
        elec_price = build_elec_price_summer(n=n_time_points)

    m = create_wrd_mp(
        n_days=n_days,
        n_time_points=n_time_points,
        elec_price=elec_price,
        daily_production_target=daily_production_target,
        total_water_production_target=total_water_production_target,
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

    solver = SolverFactory("mindtpy")
    results = solver.solve(
        m,
        strategy="OA",
        mip_solver="glpk",
        nlp_solver="ipopt",
        tee=True,
    )
    prod = [
        m.fs.mp.blocks[i].process.fs.total_water_production()
        for i in range(n_time_points)
    ]
    energy = [
        value(
            pyunits.convert(
                m.fs.mp.blocks[i].process.fs.total_water_production,
                to_units=pyunits.m**3 / pyunits.h,
            )
            * m.fs.mp.blocks[i].process.fs.treatment_energy_per_m3
        )
        for i in range(n_time_points)
    ]

    print(degrees_of_freedom(m))

    print("Total production in m3:", m.total_production())
    print("Total target water production in m3:", total_water_production_target())
    print("Total energy cost:", m.total_cost(), "2021 $")

    plot_function(m, n_time_points, season)
