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
from pyomo.opt import SolverStatus, TerminationCondition

"""
Should allow an input of a csv file of the flowrates and binary variables for each train. 
"""

if hasattr(pyunits, "USD_2021"):
    CURRENCY_UNIT = pyunits.USD_2021
elif hasattr(pyunits, "USD"):
    CURRENCY_UNIT = pyunits.USD
else:
    pyunits.load_definitions_from_strings(["USD = [currency]"])
    CURRENCY_UNIT = pyunits.USD


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

TOTAL_PLANT_PRODUCTION_CAPACITY = 53150 / 24  # m3/hr
MAX_TRAIN_FLOW = TOTAL_PLANT_PRODUCTION_CAPACITY / 4  # m3/hr
TRAIN_IDS = (1, 2, 3, 4)


def _non_overwriting_path(path):
    """Return a unique file path by appending _1, _2, ... when needed."""
    if not path.exists():
        return path

    stem = path.stem
    suffix = path.suffix
    parent = path.parent
    index = 1

    while True:
        candidate = parent / f"{stem}_{index}{suffix}"
        if not candidate.exists():
            return candidate
        index += 1


def load_train_schedule(csv_path, n_time_points=None):
    """
    Load schedule CSV for simulation.

    Required columns:
      - train_1_on ... train_4_on
      - train_1_flow_pct ... train_4_flow_pct

    Optional column:
      - time_index (sorted if present)
    """
    schedule = pd.read_csv(csv_path)

    if "time_index" in schedule.columns:
        schedule = schedule.sort_values("time_index").reset_index(drop=True)

    required_cols = [
        *[f"train_{train_id}_on" for train_id in TRAIN_IDS],
        *[f"train_{train_id}_flow_pct" for train_id in TRAIN_IDS],
    ]
    missing_cols = [col for col in required_cols if col not in schedule.columns]
    if missing_cols:
        raise ValueError(
            f"Schedule CSV is missing required columns: {missing_cols}. "
            f"Expected columns include {required_cols}."
        )

    if n_time_points is not None and len(schedule) != n_time_points:
        raise ValueError(
            f"Schedule has {len(schedule)} rows but n_time_points={n_time_points}. "
            "Provide one row per time point."
        )

    for train_id in TRAIN_IDS:
        on_col = f"train_{train_id}_on"
        pct_col = f"train_{train_id}_flow_pct"

        on_values = schedule[on_col]
        invalid_on = ~on_values.isin([0, 1])
        if invalid_on.any():
            bad_idx = invalid_on[invalid_on].index.tolist()[:5]
            raise ValueError(
                f"{on_col} must contain only 0/1 values. Invalid rows (first 5): {bad_idx}."
            )

        pct_values = pd.to_numeric(schedule[pct_col], errors="coerce")
        if pct_values.isna().any():
            bad_idx = pct_values[pct_values.isna()].index.tolist()[:5]
            raise ValueError(
                f"{pct_col} must be numeric. Invalid rows (first 5): {bad_idx}."
            )
        if ((pct_values < 0) | (pct_values > 100)).any():
            bad_idx = pct_values[
                ((pct_values < 0) | (pct_values > 100))
            ].index.tolist()[:5]
            raise ValueError(
                f"{pct_col} must be within [0, 100]. Invalid rows (first 5): {bad_idx}."
            )

        bad_off_rows = schedule[
            (schedule[on_col] == 0) & (pct_values > 0)
        ].index.tolist()
        if bad_off_rows:
            raise ValueError(
                f"{pct_col} must be 0 when {on_col}=0. Invalid rows (first 5): {bad_off_rows[:5]}."
            )

        bad_on_rows = schedule[
            (schedule[on_col] == 1) & (pct_values < 85)
        ].index.tolist()
        if bad_on_rows:
            raise ValueError(
                f"For this model, on-state trains must run at >=85% flow. "
                f"Invalid rows (first 5): {bad_on_rows[:5]}."
            )

        schedule[on_col] = schedule[on_col].astype(int)
        schedule[pct_col] = pct_values.astype(float)

    return schedule


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

    weekday_elec_price = np.ones(24)
    weekend_elec_price = np.ones(24)

    # off peak 12 AM - 4 PM
    weekday_elec_price[0:16] = off_peak_del + off_peak_gen
    weekend_elec_price[0:16] = off_peak_del + off_peak_gen
    # on peak 4 PM - 9 PM
    weekday_elec_price[16:21] = on_peak_del + on_peak_gen
    weekend_elec_price[16:21] = mid_peak_del + mid_peak_gen
    # off peak 9 PM - 12 AM
    weekday_elec_price[21:24] = off_peak_del + off_peak_gen
    weekend_elec_price[21:24] = off_peak_del + off_peak_gen

    total_hours = int(value(n))
    if total_hours <= 0:
        return np.array([]), []

    # Build repeating weekly pattern assuming the horizon starts on a weekday.
    daily_profiles = [
        weekday_elec_price,
        weekday_elec_price,
        weekday_elec_price,
        weekday_elec_price,
        weekday_elec_price,
        weekend_elec_price,
        weekend_elec_price,
    ]

    full_days, rem_hours = divmod(total_hours, 24)

    day_blocks = []
    for day_idx in range(full_days):
        day_blocks.append(daily_profiles[day_idx % 7])

    if day_blocks:
        elec_price = np.concatenate(day_blocks)
    else:
        elec_price = np.array([])

    if rem_hours > 0:
        next_day_profile = daily_profiles[full_days % 7]
        elec_price = np.concatenate([elec_price, next_day_profile[:rem_hours]])

    # Absolute hourly indices that fall in weekday on-peak window (16:00-21:00).
    peak_hours = []
    for h in range(total_hours):
        day_of_week = (h // 24) % 7
        hour_of_day = h % 24
        if day_of_week < 5 and 16 <= hour_of_day <= 20:
            peak_hours.append(h)

    return elec_price, peak_hours


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
        peak_hours = [h for h in range(int(value(n))) if (h % 24) in range(16, 21)]

    return elec_price, peak_hours


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

    total_plant_production_capacity = (
        TOTAL_PLANT_PRODUCTION_CAPACITY * pyunits.m**3 / pyunits.h
    )  # m3 per hour
    train_production_capacity = (
        total_plant_production_capacity / 4
    )  # m3 per hour per train

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

    m.fs.uf_on = Var(
        initialize=1,
        domain=Binary,
        doc="Binary variable indicating if UF is on (1 if any RO train is on)",
    )

    @m.Constraint(TRAIN_IDS, doc="UF must be on if any RO train is on")
    def eq_uf_on_lb(b, train_id):
        return b.fs.uf_on >= getattr(b.fs, f"train_{train_id}_on")

    @m.Constraint(doc="UF can only be on when at least one RO train is on")
    def eq_uf_on_ub(b):
        return b.fs.uf_on <= sum(
            getattr(b.fs, f"train_{train_id}_on") for train_id in TRAIN_IDS
        )

    # Constraint defining the flowrates based on input file
    @m.Constraint(doc="Upper bound for flow depends on the binary variable")
    def eq_train_1_ub(b):
        return (
            b.fs.water_production_ro_train_1
            <= train_production_capacity * b.fs.train_1_on
        )

    @m.Constraint(doc="Lower bound for flow depends on the binary variable")
    def eq_train_1_lb(b):
        return (
            b.fs.water_production_ro_train_1
            >= train_production_capacity * 0.85 * b.fs.train_1_on
        )

    @m.Constraint(doc="Upper bound for flow depends on the binary variable")
    def eq_train_2_ub(b):
        return (
            b.fs.water_production_ro_train_2
            <= train_production_capacity * b.fs.train_2_on
        )

    @m.Constraint(doc="Lower bound for flow depends on the binary variable")
    def eq_train_2_lb(b):
        return (
            b.fs.water_production_ro_train_2
            >= train_production_capacity * 0.85 * b.fs.train_2_on
        )

    @m.Constraint(doc="Upper bound for flow depends on the binary variable")
    def eq_train_3_ub(b):
        return (
            b.fs.water_production_ro_train_3
            <= train_production_capacity * b.fs.train_3_on
        )

    @m.Constraint(doc="Lower bound for flow depends on the binary variable")
    def eq_train_3_lb(b):
        return (
            b.fs.water_production_ro_train_3
            >= train_production_capacity * 0.85 * b.fs.train_3_on
        )

    @m.Constraint(doc="Upper bound for flow depends on the binary variable")
    def eq_train_4_ub(b):
        return (
            b.fs.water_production_ro_train_4
            <= train_production_capacity * b.fs.train_4_on
        )

    @m.Constraint(doc="Lower bound for flow depends on the binary variable")
    def eq_train_4_lb(b):
        return (
            b.fs.water_production_ro_train_4
            >= train_production_capacity * 0.85 * b.fs.train_4_on
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
        bounds=(0, 3000),
        units=pyunits.kWh / pyunits.h,
        doc="Total treatment energy required per hour",
    )

    def calculate_ro_power(flow, train_on):
        # Linear train power model (kW): 0 when train is off, fitted line when on
        return (
            0.5207 * flow / (pyunits.m**3 / pyunits.hr) - 95.74 * train_on
        ) * pyunits.kW

    def calculate_uf_power(flow, uf_on):
        # Linear UF SYSTEM power model (kW)
        # If train_1_on == 0, then the whole system is off and no power use from UF
        return (0.199 * flow / (pyunits.m**3 / pyunits.hr) - 27.4 * uf_on) * pyunits.kW

    def calculate_UVAOP_power(flow):
        # Linear UVAOP power model (kW)
        return (0.101 * flow / (pyunits.m**3 / pyunits.hr)) * pyunits.kW

    m.fs.ro_train_1_energy_rate = Var(
        initialize=0,
        bounds=(0, 1000),
        units=pyunits.kWh / pyunits.h,
        doc="Energy rate for RO train 1",
    )
    m.fs.ro_train_2_energy_rate = Var(
        initialize=0,
        bounds=(0, 1000),
        units=pyunits.kWh / pyunits.h,
        doc="Energy rate for RO train 2",
    )
    m.fs.ro_train_3_energy_rate = Var(
        initialize=0,
        bounds=(0, 1000),
        units=pyunits.kWh / pyunits.h,
        doc="Energy rate for RO train 3",
    )
    m.fs.ro_train_4_energy_rate = Var(
        initialize=0,
        bounds=(0, 1000),
        units=pyunits.kWh / pyunits.h,
        doc="Energy rate for RO train 4",
    )
    m.fs.uf_energy_rate = Var(
        initialize=0,
        bounds=(0, 1000),
        units=pyunits.kWh / pyunits.h,
        doc="Energy rate for UF system",
    )
    m.fs.uvaop_energy_rate = Var(
        initialize=0,
        bounds=(0, 1000),
        units=pyunits.kWh / pyunits.h,
        doc="Energy rate for UV/AOP system",
    )

    @m.Constraint(doc="RO train 1 energy rate")
    def eq_ro_train_1_energy_rate(b):
        return b.fs.ro_train_1_energy_rate == calculate_ro_power(
            b.fs.water_production_ro_train_1, b.fs.train_1_on
        )

    @m.Constraint(doc="RO train 2 energy rate")
    def eq_ro_train_2_energy_rate(b):
        return b.fs.ro_train_2_energy_rate == calculate_ro_power(
            b.fs.water_production_ro_train_2, b.fs.train_2_on
        )

    @m.Constraint(doc="RO train 3 energy rate")
    def eq_ro_train_3_energy_rate(b):
        return b.fs.ro_train_3_energy_rate == calculate_ro_power(
            b.fs.water_production_ro_train_3, b.fs.train_3_on
        )

    @m.Constraint(doc="RO train 4 energy rate")
    def eq_ro_train_4_energy_rate(b):
        return b.fs.ro_train_4_energy_rate == calculate_ro_power(
            b.fs.water_production_ro_train_4, b.fs.train_4_on
        )

    @m.Constraint(doc="UF system energy rate")
    def eq_uf_energy_rate(b):
        return b.fs.uf_energy_rate == calculate_uf_power(
            b.fs.total_water_production, b.fs.uf_on
        )

    @m.Constraint(doc="UV/AOP energy rate")
    def eq_uvaop_energy_rate(b):
        return b.fs.uvaop_energy_rate == calculate_UVAOP_power(
            b.fs.total_water_production
        )

    @m.Constraint(doc="Total treatment energy rate is sum of component energy rates")
    def eq_treatment_energy_rate_sum(b):
        return b.fs.treatment_energy_rate == (
            b.fs.ro_train_1_energy_rate
            + b.fs.ro_train_2_energy_rate
            + b.fs.ro_train_3_energy_rate
            + b.fs.ro_train_4_energy_rate
            + b.fs.uf_energy_rate
            + b.fs.uvaop_energy_rate
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


# def unfix_dof(m):
#     # Train 1 and 2 are always on, so we only vary the fraction of water treated by train 3 and 4
#     m.fs.water_production_ro_train_3.unfix()
#     m.fs.water_production_ro_train_4.unfix()
#     m.fs.train_3_on.unfix()
#     m.fs.train_4_on.unfix()
#     return None


def initialize_mp(m):
    print("Initializing multi-period model...")
    # Check if first time step
    max_train_flow = MAX_TRAIN_FLOW
    m.fs.water_production_ro_train_1.fix(max_train_flow)
    m.fs.water_production_ro_train_2.fix(max_train_flow)
    m.fs.water_production_ro_train_3.fix(max_train_flow * 0.85)
    m.fs.water_production_ro_train_4.fix(max_train_flow * 0.85)

    m.fs.train_1_on.fix(1)
    m.fs.train_2_on.fix(1)
    m.fs.train_3_on.fix(1)
    m.fs.train_4_on.fix(1)


def apply_train_schedule(process_model, schedule_row):
    for train_id in TRAIN_IDS:
        on_col = f"train_{train_id}_on"
        pct_col = f"train_{train_id}_flow_pct"

        train_on = int(schedule_row[on_col])
        flow_pct = float(schedule_row[pct_col])
        flow_value = MAX_TRAIN_FLOW * flow_pct / 100.0

        getattr(process_model.fs, f"train_{train_id}_on").fix(train_on)
        getattr(process_model.fs, f"water_production_ro_train_{train_id}").fix(
            flow_value
        )


def create_wrd_mp(
    n_days=1,
    n_time_points=24,
    elec_price=elec_price,
    train_schedule=None,
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

    m.fs.mp = MultiPeriodModel(
        n_time_points=n_time_points,
        process_model_func=build_wrd_flowsheet,
        linking_variable_func=get_wrd_variable_pairs,
        initialization_func=None,
        unfix_dof_func=None,
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

    m.fs.mp.build_multi_period_model(
        model_data_kwargs=flowsheet_options,
        flowsheet_options=flowsheet_options,
        initialization_options=None,
        unfix_dof_options=None,
    )

    if train_schedule is None:
        raise ValueError("train_schedule must be provided for simulation mode.")

    if len(train_schedule) != n_time_points:
        raise ValueError(
            f"train_schedule has {len(train_schedule)} rows but n_time_points={n_time_points}."
        )

    for t in range(n_time_points):
        initialize_mp(m.fs.mp.blocks[t].process)
        apply_train_schedule(m.fs.mp.blocks[t].process, train_schedule.iloc[t])

    m.fs.mp.blocks[0].process.fs.pre_acc_production.fix(0)
    m.fs.mp.blocks[0].process.fs.pre_acc_energy.fix(0)

    m.fs.fixed_demand_price = Param(
        initialize=demand_charges["fixed_demand_price"],
        mutable=True,
        units=CURRENCY_UNIT / pyunits.kW,
        doc="Demand charge associated with highest period of energy use",
    )

    m.fs.on_peak_demand_price = Param(
        initialize=demand_charges["on_peak_demand_price"],
        mutable=True,
        units=CURRENCY_UNIT / pyunits.kW,
        doc="Demand charge associated with greatest power value within the on-peak hours",
    )

    m.fs.highest_demand = Var(
        initialize=1000,
        bounds=(0, None),
        units=pyunits.kW,
        doc="Demand during highest period of energy use",
    )

    m.fs.highest_on_peak_demand = Var(
        initialize=1000,
        bounds=(0, None),
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
        return b.fs.highest_demand >= b.fs.mp.blocks[i].process.fs.treatment_energy_rate

    m.fs.peak_hours = peak_hours

    @m.Constraint(
        m.fs.peak_hours,
        doc="Upper bound highest on-peak demand by each period energy rate during peak hours",
    )
    def eq_highest_on_peak_demand(b, i):
        return (
            b.fs.highest_on_peak_demand
            >= b.fs.mp.blocks[i].process.fs.treatment_energy_rate
        )

    # @m.Constraint(range(n_time_points), doc="Production should not change at midnight")
    # def eq_midnight_constraint_on_off(b, h):
    #     if h == 0 or (h % 24) != 0:
    #         return Constraint.Skip
    #     return (
    #         b.fs.mp.blocks[h].process.fs.train_3_on
    #         == b.fs.mp.blocks[h - 1].process.fs.train_3_on
    #     )
    # Adding working hours, but not including the constraints because it's just a simulation.
    non_working_hours_morning = [0, 1, 2, 3, 4, 5, 6, 7, 8]
    non_working_hours_evening = [18, 19, 20, 21, 22, 23]

    # m.fs.non_working_hours_morning = [h for h in range(n_time_points) if (h % 24) in non_working_hours_morning]
    # m.fs.non_working_hours_evening = [h for h in range(n_time_points) if (h % 24) in non_working_hours_evening]

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

    @m.Expression(doc="Total production")
    def total_production(b):
        return sum(
            [
                b.fs.mp.blocks[i].process.fs.total_water_production
                * b.fs.mp.blocks[i].process.fs.time_step
                for i in range(n_time_points)
            ]
        )

    # Set objective
    m.fs.obj = Objective(expr=m.total_cost)

    return m


def plot_function(m, n_time_points, season):
    time = np.linspace(0, n_time_points - 1, n_time_points)
    fig = plt.figure(figsize=(8, 8))
    gs = fig.add_gridspec(2, 1, height_ratios=[1, 1])
    ax_energy = fig.add_subplot(gs[0])
    ax_trains = fig.add_subplot(gs[1], sharex=ax_energy)
    ax_energy.set_facecolor("#f5f5f5")
    ax_trains.set_facecolor("#f5f5f5")

    peak_hours = set(getattr(m.fs, "peak_hours", []))
    if peak_hours:
        peak_legend_added = False
        for i in range(n_time_points):
            if i in peak_hours:
                span_label = "Peak Hours" if not peak_legend_added else None
                ax_energy.axvspan(
                    i,
                    i + 1,
                    color="grey",
                    alpha=0.2,
                    linewidth=0,
                    zorder=-1,
                    hatch="///",
                    label=span_label,
                )
                ax_trains.axvspan(
                    i,
                    i + 1,
                    color="grey",
                    alpha=0.2,
                    linewidth=0,
                    zorder=-1,
                    hatch="///",
                    label=span_label,
                )
                peak_legend_added = True

    total_energy = [
        value(
            pyunits.convert(
                m.fs.mp.blocks[i].process.fs.treatment_energy_rate,
                to_units=pyunits.kWh / pyunits.h,
            )
        )
        for i in range(n_time_points)
    ]
    ro1_energy = [
        value(
            pyunits.convert(
                m.fs.mp.blocks[i].process.fs.ro_train_1_energy_rate,
                to_units=pyunits.kWh / pyunits.h,
            )
        )
        for i in range(n_time_points)
    ]
    ro2_energy = [
        value(
            pyunits.convert(
                m.fs.mp.blocks[i].process.fs.ro_train_2_energy_rate,
                to_units=pyunits.kWh / pyunits.h,
            )
        )
        for i in range(n_time_points)
    ]
    ro3_energy = [
        value(
            pyunits.convert(
                m.fs.mp.blocks[i].process.fs.ro_train_3_energy_rate,
                to_units=pyunits.kWh / pyunits.h,
            )
        )
        for i in range(n_time_points)
    ]
    ro4_energy = [
        value(
            pyunits.convert(
                m.fs.mp.blocks[i].process.fs.ro_train_4_energy_rate,
                to_units=pyunits.kWh / pyunits.h,
            )
        )
        for i in range(n_time_points)
    ]
    uf_energy = [
        value(
            pyunits.convert(
                m.fs.mp.blocks[i].process.fs.uf_energy_rate,
                to_units=pyunits.kWh / pyunits.h,
            )
        )
        for i in range(n_time_points)
    ]
    post_treatment_energy = [
        value(
            pyunits.convert(
                m.fs.mp.blocks[i].process.fs.uvaop_energy_rate,
                to_units=pyunits.kWh / pyunits.h,
            )
        )
        for i in range(n_time_points)
    ]

    ax_energy.stackplot(
        time + 0.5,
        ro1_energy,
        ro2_energy,
        ro3_energy,
        ro4_energy,
        uf_energy,
        post_treatment_energy,
        labels=[
            "RO Train 1",
            "RO Train 2",
            "RO Train 3",
            "RO Train 4",
            "UF System",
            "Post-Treatment",
        ],
        alpha=0.5,
    )

    ax_energy.plot(
        time + 0.5,
        total_energy,
        label="Total Power",
        color="black",
        linestyle="--",
        linewidth=2,
    )
    ax_energy.set_ylim(0, 2500)
    ax_energy.set_ylabel("kW", fontsize=12)
    ax_energy.grid(False)

    ax_price = ax_energy.twinx()
    elec_price = [
        value(
            pyunits.convert(
                m.fs.mp.blocks[i].process.fs.electricity_price,
                to_units=CURRENCY_UNIT / pyunits.kWh,
            )
        )
        for i in range(n_time_points)
    ]
    ax_price.plot(
        time + 0.5,
        np.asarray(elec_price) * 100,
        label="Elec. Price",
        color="orange",
        linewidth=2,
    )
    ax_price.set_ylabel("cent/kWh", fontsize=12)
    ax_price.set_ylim(0, max(np.asarray(elec_price) * 100) + 3)

    handle1, label1 = ax_energy.get_legend_handles_labels()
    handle2, label2 = ax_price.get_legend_handles_labels()
    handles = handle2 + handle1
    labels = label2 + label1
    leg1 = ax_price.legend(
        handles,
        labels,
        loc="lower left",
        bbox_to_anchor=(0.0, 1, 1, 1),
        framealpha=1.0,
        ncol=4,
        fontsize=12,
        mode="expand",
    )
    leg1.set_zorder(1000)
    leg1.get_frame().set_facecolor("white")
    ax_energy.xaxis.set_major_locator(plt.MaxNLocator(24))

    prod = [
        value(
            pyunits.convert(
                m.fs.mp.blocks[i].process.fs.total_water_production,
                to_units=pyunits.m**3 / pyunits.h,
            )
        )
        for i in range(n_time_points)
    ]
    ax_trains.plot(
        time + 0.5,
        prod,
        label="Water Production",
        color="black",
        linestyle="--",
        linewidth=2,
        alpha=0.75,
    )
    ax_trains.set_ylim(0, 2500)
    ax_trains.axhline(
        y=TOTAL_PLANT_PRODUCTION_CAPACITY,
        color="blue",
        linestyle=":",
        linewidth=2,
        alpha=0.75,
        label="Max Production",
        zorder=0,
    )
    ax_trains.set_ylabel("m$^3$/h", fontsize=12)
    ax_trains.set_xlabel("Hours", fontsize=12)
    ax_trains.xaxis.set_major_locator(plt.MaxNLocator(24))
    ax_trains.grid(False)

    train_1_flows = [
        value(
            pyunits.convert(
                m.fs.mp.blocks[i].process.fs.water_production_ro_train_1,
                to_units=pyunits.m**3 / pyunits.h,
            )
        )
        for i in range(n_time_points)
    ]
    train_2_flows = [
        value(
            pyunits.convert(
                m.fs.mp.blocks[i].process.fs.water_production_ro_train_2,
                to_units=pyunits.m**3 / pyunits.h,
            )
        )
        for i in range(n_time_points)
    ]
    train_3_flows = [
        value(
            pyunits.convert(
                m.fs.mp.blocks[i].process.fs.water_production_ro_train_3,
                to_units=pyunits.m**3 / pyunits.h,
            )
        )
        for i in range(n_time_points)
    ]
    train_4_flows = [
        value(
            pyunits.convert(
                m.fs.mp.blocks[i].process.fs.water_production_ro_train_4,
                to_units=pyunits.m**3 / pyunits.h,
            )
        )
        for i in range(n_time_points)
    ]

    ax_trains.stackplot(
        time + 0.5,
        train_1_flows,
        train_2_flows,
        train_3_flows,
        train_4_flows,
        labels=["RO Train 1", "RO Train 2", "RO Train 3", "RO Train 4"],
        alpha=0.5,
    )

    handle_t, label_t = ax_trains.get_legend_handles_labels()
    leg3 = ax_trains.legend(
        handle_t,
        label_t,
        loc="lower left",
        bbox_to_anchor=(0.0, 1, 1, 1),
        framealpha=1.0,
        ncol=3,
        fontsize=12,
        mode="expand",
    )
    leg3.get_frame().set_facecolor("white")

    for a in (ax_energy, ax_trains):
        a.set_xlim(0, n_time_points)
        a.xaxis.set_major_locator(plt.MaxNLocator(24))

    for a in (ax_energy, ax_price, ax_trains):
        a.tick_params(axis="both", labelsize=14)
    for label in ax_trains.get_xticklabels():
        label.set_rotation(45)
        label.set_ha("center")
    for label in ax_energy.get_xticklabels():
        label.set_rotation(45)
        label.set_ha("center")

    fig.tight_layout()
    output_path = _non_overwriting_path(
        Path(__file__).resolve().parent
        / f"{Path(__file__).stem}_wrd_{season}_result.png"
    )
    fig.savefig(output_path, dpi=600)
    plt.show()


def plot_function_top_w_real_energy_use(m, n_time_points, season):

    time = np.linspace(0, n_time_points - 1, n_time_points)

    sim_energy_profile = [
        value(
            pyunits.convert(
                m.fs.mp.blocks[i].process.fs.treatment_energy_rate,
                to_units=pyunits.kWh / pyunits.h,
            )
        )
        for i in range(n_time_points)
    ]
    electricity_cost = [
        value(
            pyunits.convert(
                m.fs.mp.blocks[i].process.fs.electricity_price,
                to_units=CURRENCY_UNIT / pyunits.kWh,
            )
        )
        for i in range(n_time_points)
    ]

    act_energy_profile = pd.read_csv(r"C:\Users\rchurchi\flex_desal\src\wrd\multiperiod\figures\Aug_21_kW_hourly.csv")["total_energy_kW"].to_list()

    fig, ax = plt.subplots(1, 1, figsize=(12, 6))

    sim_energy_line = ax.plot(
        time + 0.5,
        sim_energy_profile,
        label= "Modeled Energy Consumption (kWh)",
        color="orange",
        marker="o",
    )

    act_energy_line = ax.plot(
        time + 0.5,
        act_energy_profile,
        label= "Energy Consumption Data (kWh)",
        color="blue",
        marker="s",
    )

    ax.set_ylim(0, 2500)
    ax.set_ylabel("Energy Consumption (kWh)", fontsize=16)
    ax.set_xlabel("Hours", fontsize=16)
    ax.set_title("Energy Consumption and Cost - August 2021", fontsize=18, fontweight="bold")
    ax.grid(False)
    ax.xaxis.set_major_locator(plt.MaxNLocator(24))

    ax_grid = ax.twinx()
    electricity_cost_line = ax_grid.plot(
        time + 0.5,
        electricity_cost,
        label="Electricity Cost ($/kWh)",
        color="black",
        linestyle="-",
        linewidth=2,
    )
    ax_grid.set_ylabel("Electricity Cost ($/kWh)", fontsize=16)
    ax_grid.set_ylim(0, 0.17)

    ax_grid.legend(
        handles=[electricity_cost_line[0], sim_energy_line[0], act_energy_line[0]],
        loc="lower left",
        framealpha=1.0,
        ncol=1,
        fontsize=11,
    )

    for a in (ax, ax_grid):
        a.set_xlim(0, n_time_points)
        a.tick_params(axis="both", labelsize=11)

    fig.tight_layout()
    output_path = _non_overwriting_path(
        Path(__file__).resolve().parent
        / f"{Path(__file__).stem}_comparison_to_plant_data.png"
    )
    fig.savefig(output_path, dpi=600)
    plt.show()


def plot_function_top_and_num_trains(m, n_time_points, season):
    time = np.linspace(0, n_time_points - 1, n_time_points)

    fig, (ax, ax_trains) = plt.subplots(
        2,
        1,
        figsize=(12, 12),
        gridspec_kw={"height_ratios": [1, 1]},
    )

    # First subplot: Energy consumption (left) + Grid cost (right)
    electricity_cost = [
        value(
            pyunits.convert(
                m.fs.mp.blocks[i].process.fs.electricity_price,
                to_units=CURRENCY_UNIT / pyunits.kWh,
            )
        )
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

    ax.plot(
        time + 0.5, energy, label="Energy Consumption (kWh)", color="orange", marker="o"
    )
    ax.set_ylim(0, 2500)
    ax.set_ylabel("Energy Consumption (kWh)", fontsize=14)
    ax.set_xlabel("Hours", fontsize=16)
    ax.set_title("Flexible Operations Scenario", fontsize=14, fontweight="bold")
    ax.grid(False)

    ax_grid = ax.twinx()
    ax_grid.plot(
        time + 0.5,
        electricity_cost,
        label="Electricity Cost ($/kWh)",
        color="black",
        linestyle="-",
        linewidth=2,
    )
    ax_grid.set_ylabel("Electricity Cost ($/kWh)", fontsize=14)
    ax_grid.set_ylim(0, 0.17)

    handle1, label1 = ax.get_legend_handles_labels()
    handle2, label2 = ax_grid.get_legend_handles_labels()
    handles = [handle2[0], handle1[0]]
    labels = [label2[0], label1[0]]
    leg = ax_grid.legend(
        handles, labels, loc="lower left", framealpha=1.0, ncol=2, fontsize=14
    )
    leg.set_zorder(1000)
    leg.get_frame().set_facecolor("white")

    ax.xaxis.set_major_locator(plt.MaxNLocator(24))

    # Second subplot: Water production (left) + Number of trains in operation (right)
    prod = [
        m.fs.mp.blocks[i].process.fs.total_water_production()
        for i in range(n_time_points)
    ]
    num_trains = [
        value(m.fs.mp.blocks[i].process.fs.train_1_on())
        + value(m.fs.mp.blocks[i].process.fs.train_2_on())
        + value(m.fs.mp.blocks[i].process.fs.train_3_on())
        + value(m.fs.mp.blocks[i].process.fs.train_4_on())
        for i in range(n_time_points)
    ]

    ax_trains.plot(
        time + 0.5,
        prod,
        label="Water Production (m$^3$/h)",
        color="blue",
        marker="o",
        linewidth=2,
    )
    ax_trains.set_ylim(0, 2250)
    ax_trains.set_ylabel("Water Production (m$^3$/h)", fontsize=14)
    ax_trains.set_xlabel("Hours", fontsize=14)
    ax_trains.grid(False)

    ax_trains2 = ax_trains.twinx()
    ax_trains2.plot(
        time + 0.5, num_trains, label="RO Trains in Operation", color="red", marker="s", linewidth=2
    )
    ax_trains2.set_ylim(0, 5)
    ax_trains2.set_ylabel("RO Trains in Operation", fontsize=14)
    ax_trains2.set_yticks([0, 1, 2, 3, 4])

    handle_prod, label_prod = ax_trains.get_legend_handles_labels()
    handle_trains, label_trains = ax_trains2.get_legend_handles_labels()
    leg2 = ax_trains2.legend(
        handle_prod + handle_trains,
        label_prod + label_trains,
        loc="lower left",
        fontsize=14,
        framealpha=1.0,
        ncol=2,
    )
    leg2.get_frame().set_facecolor("white")

    ax_trains.xaxis.set_major_locator(plt.MaxNLocator(24))

    for a in (ax, ax_trains):
        a.set_xlim(0, n_time_points)
        a.xaxis.set_major_locator(plt.MaxNLocator(24))
        a.tick_params(axis="both", labelsize=14)

    ax_grid.tick_params(axis="both", labelsize=14)
    ax_trains2.tick_params(axis="both", labelsize=14)

    fig.tight_layout()
    output_path = _non_overwriting_path(
        Path(__file__).resolve().parent
        / f"{Path(__file__).stem}_wrd_{season}_top_and_num_trains.png"
    )
    fig.savefig(output_path, dpi=600)
    plt.show()
    

def print_unfixed_vars(model):
    print("Unfixed variables contributing to degrees of freedom:")
    for v in model.component_data_objects(ctype=Var, descend_into=True):
        if not v.fixed:
            print(f"  {v.name}")


if __name__ == "__main__":
    schedule_csv = (
        r"src\wrd\multiperiod\simulation\summer_modular_sim.csv"
    )
    season = "summer"
    tee = True

    train_schedule = load_train_schedule(schedule_csv)
    n_time_points = len(train_schedule)
    if n_time_points == 0:
        raise ValueError("Schedule CSV must contain at least one row.")
    n_days = n_time_points / 24

    if season == "winter":
        elec_price, peak_hours = build_elec_price_winter(n=n_time_points)
        demand_charges = {
            "fixed_demand_price": 19.62,
            "on_peak_demand_price": 7.99 + 2.55,
        }  # March 2023
    else:
        elec_price, peak_hours = build_elec_price_summer(n=n_time_points)
        demand_charges = {
            "fixed_demand_price": 19.94,
            "on_peak_demand_price": 22.10 + 14.68,
        }  # June 2023

    m = create_wrd_mp(
        n_days=n_days,
        n_time_points=n_time_points,
        elec_price=elec_price,
        train_schedule=train_schedule,
        peak_hours=peak_hours,
        demand_charges=demand_charges,
    )
    assert_units_consistent(m)
    # print_unfixed_vars(m)
    # Add a water demand value
    aug_total_water = 1324527
    m.total_demand = Param(
        initialize=aug_total_water * 0.98,
        mutable=True,
        units=pyunits.m**3,
        doc="Total water demand in m3",
    )

    @m.Constraint(doc="Total production must meet demand")
    def eq_total_production(b):
        return b.total_production >= b.total_demand

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
    results = solver.solve(m, tee=tee)

    solver_status = results.solver.status
    termination = results.solver.termination_condition
    if not (
        solver_status == SolverStatus.ok
        and termination in (TerminationCondition.optimal, TerminationCondition.feasible)
    ):
        raise RuntimeError(
            "Simulation solve failed. "
            f"Solver status={solver_status}, termination={termination}. "
            "This often indicates an inconsistent schedule row "
            "(e.g., train_on=1 with flow_pct<85)."
        )

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
    # plot_function_top_and_num_trains(m, n_time_points, season) # Outdated formating
    # plot_function_top_w_real_energy_use(m, n_time_points, season)