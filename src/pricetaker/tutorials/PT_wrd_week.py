##########
#
# Pricetaker Implementation for WRD plant using surrogate models for UF and RO energy consumption.
# Short = Only testing one or two days, and no rain events or demand response events.
#
##########

import warnings
import logging

warnings.filterwarnings("ignore", message=".*implicit domain of 'Any'.*")
logging.getLogger("pyomo").setLevel(logging.ERROR)

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from datetime import datetime
from pathlib import Path

import pyomo.environ as pyo
from pyomo.environ import SolverFactory, value

from watertap.flowsheets.flex_desal import wrd_ro_flowsheet as fs
from watertap.flowsheets.flex_desal import utils
from watertap.flowsheets.flex_desal.params import FlexDesalParams
from watertap.core.solvers import get_solver

from idaes.core.util.model_diagnostics import DiagnosticsToolbox
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.apps.grid_integration import PriceTakerModel


def plot_function(m, n_time_points, output_stem, peak_hours=None):
    time = np.linspace(0, n_time_points - 1, n_time_points)
    fig = plt.figure(figsize=(8, 8))
    gs = fig.add_gridspec(2, 1, height_ratios=[1, 1])
    ax_energy = fig.add_subplot(gs[0])
    ax_trains = fig.add_subplot(gs[1], sharex=ax_energy)
    ax_energy.set_facecolor("#f5f5f5")
    ax_trains.set_facecolor("#f5f5f5")

    if peak_hours is not None:
        peak_legend_added = False
        for i, is_peak in enumerate(peak_hours):
            if is_peak:
                # Shade full hourly intervals where variable demand charges apply.
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

    # First subplot: Stacked energy consumption by major equipment
    total_energy = [
        pyo.value(v[None])
        for v in m.period[:, :].net_power_consumption.extract_values()
    ]

    ro1_energy = [
        v[None]
        for v in m.period[:, :]
        .reverse_osmosis.ro_skid[1]
        .power_consumption.extract_values()
    ]
    ro2_energy = [
        v[None]
        for v in m.period[:, :]
        .reverse_osmosis.ro_skid[2]
        .power_consumption.extract_values()
    ]
    ro3_energy = [
        v[None]
        for v in m.period[:, :]
        .reverse_osmosis.ro_skid[3]
        .power_consumption.extract_values()
    ]
    ro4_energy = [
        v[None]
        for v in m.period[:, :]
        .reverse_osmosis.ro_skid[4]
        .power_consumption.extract_values()
    ]

    uf1_energy = [
        v[None]
        for v in m.period[:, :]
        .pretreatment.uf_pumps[1]
        .power_consumption.extract_values()
    ]
    uf2_energy = [
        v[None]
        for v in m.period[:, :]
        .pretreatment.uf_pumps[2]
        .power_consumption.extract_values()
    ]
    uf3_energy = [
        v[None]
        for v in m.period[:, :]
        .pretreatment.uf_pumps[3]
        .power_consumption.extract_values()
    ]

    other_energy = np.array(total_energy) - (
        np.array(ro1_energy)
        + np.array(ro2_energy)
        + np.array(ro3_energy)
        + np.array(ro4_energy)
        + np.array(uf1_energy)
        + np.array(uf2_energy)
        + np.array(uf3_energy)
    )
    # Clip tiny negatives from solver tolerances so stackplot remains well-defined.
    other_energy = np.maximum(other_energy, 0.0)

    ax_energy.stackplot(
        time + 0.5,
        ro1_energy,
        ro2_energy,
        ro3_energy,
        ro4_energy,
        uf1_energy,
        uf2_energy,
        uf3_energy,
        other_energy,
        labels=[
            "RO Train 1",
            "RO Train 2",
            "RO Train 3",
            "RO Train 4",
            "UF Pump 1",
            "UF Pump 2",
            "UF Pump 3",
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
    elec_price = np.asarray(m._config.lmp_data, dtype=float)
    if elec_price.size != n_time_points:
        if elec_price.size % n_time_points == 0:
            elec_price = elec_price.reshape(n_time_points, -1).mean(axis=1)
        else:
            elec_price = elec_price[:n_time_points]
    ax_price.plot(
        time + 0.5,
        elec_price * 100,
        label="Elec. Price",
        color="orange",
        linewidth=2,
    )
    ax_price.set_ylabel("¢/kWh", fontsize=12)
    ax_price.set_ylim(0, max(elec_price * 100) + 3)

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

    # Second subplot: Water production and RO train flow rates
    prod = [
        v[None] for v in m.period[:, :].posttreatment.product_flowrate.extract_values()
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
        y=602 * 4,
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

    # Extract RO train flow rates (m3/hr) for stacked plotting
    train_1_flows = [
        v[None]
        for v in m.period[:, :]
        .reverse_osmosis.ro_skid[1]
        .product_flowrate.extract_values()
    ]
    train_2_flows = [
        v[None]
        for v in m.period[:, :]
        .reverse_osmosis.ro_skid[2]
        .product_flowrate.extract_values()
    ]
    train_3_flows = [
        v[None]
        for v in m.period[:, :]
        .reverse_osmosis.ro_skid[3]
        .product_flowrate.extract_values()
    ]
    train_4_flows = [
        v[None]
        for v in m.period[:, :]
        .reverse_osmosis.ro_skid[4]
        .product_flowrate.extract_values()
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

    # Set consistent x-axis limits and formatting
    for a in (ax_energy, ax_trains):
        a.set_xlim(0, n_time_points)
        a.xaxis.set_major_locator(plt.MaxNLocator(24))

    # Tick labels for all axes
    for a in (ax_energy, ax_price, ax_trains):
        a.tick_params(axis="both", labelsize=14)
    for label in ax_trains.get_xticklabels():
        label.set_rotation(45)
        label.set_ha("center")
    for label in ax_energy.get_xticklabels():
        label.set_rotation(45)
        label.set_ha("center")

    fig.tight_layout()
    fig.savefig(f"{output_stem}.png", dpi=600)
    plt.show()


def _fix_nominal_flowrates(m):
    # TODO: This function is not working correctly because the optimal solution gives minimum and maximum flowrates that are not equal to the nominal flowrate. This is because the surrogate model is not perfectly accurate and the optimal solution is not exactly at the nominal flowrate.
    m.params.wrd_ro.minimum_flowrate = m.params.wrd_ro.nominal_flowrate
    m.params.wrd_ro.maximum_flowrate = m.params.wrd_ro.nominal_flowrate


def _restrict_flexible_trains(m, num_flexible_trains):
    ro_skids = sorted(list(m.period[1, 1].reverse_osmosis.set_ro_skids))
    n_ro_skids = len(ro_skids)

    if num_flexible_trains < 0 or num_flexible_trains > n_ro_skids:
        raise ValueError(
            "Invalid num_flexible_trains "
            f"'{num_flexible_trains}'. Valid range is [0, {n_ro_skids}]."
        )
    if num_flexible_trains == 0:
        # Fix the 4th train to be off
        for p in m.period:
            ro_skid = m.period[p].reverse_osmosis.ro_skid[4]
            ro_skid.feed_flowrate.fix(0)
            ro_skid.recovery.fix(m.params.wrd_ro.nominal_recovery)
            ro_skid.op_mode.fix(0)
        # Fix all the other skids to the nominal flowrate and recovery
        num_flexible_trains = 1

    non_flexible_skids = ro_skids[: n_ro_skids - num_flexible_trains]

    # This would be restricting just the modularity of the trains
    # for p in m.period:
    #     for skid in non_flexible_skids:
    #         ro_skid = m.period[p].reverse_osmosis.ro_skid[skid]
    #         ro_skid.startup.fix(0)
    #         ro_skid.shutdown.fix(0)

    for p in m.period:
        for skid in non_flexible_skids:
            ro_skid = m.period[p].reverse_osmosis.ro_skid[skid]
            ro_skid.feed_flowrate.fix(605.3)
            ro_skid.recovery.fix(m.params.wrd_ro.nominal_recovery)


def _fix_operations_for_first_four_days(m, peak_hours=None):
    """Fix all RO trains to expected behavior for first four days. This could be some part of an initialization strat. to improve solve times."""
    for d, p in m.period:
        if p <= 4 * 24:  # Assuming hourly time steps
            if p <= 2:
                # Avoiding constraint that plant has to be on at first (and therefore last) time step.
                pass
            elif peak_hours is not None and peak_hours[p]:
                # Full shutdown during peak hours. Could also consider just shutting down two RO skids during peak hours
                # This is too strong to impose on model. Turning off during peak hours should be an output of the opt., not prescribed.
                # m.period[d, p].reverse_osmosis.ro_skid[1].op_mode.fix(0)
                m.period[d, p].reverse_osmosis.ro_skid[4].op_mode.fix(
                    0
                )  # 4th skid off during peak hours. If 0 flex skids, forces this train off. But that should be ok for cases we are looking at.
                pass
            else:
                # Just ensure plant is on during the non-peak hours
                m.period[d, p].reverse_osmosis.ro_skid[1].op_mode.fix(
                    1
                )  # Plant must be on


def _begin_and_end_constraint(m):
    """Force RO train 1 op_mode to match between first and last timesteps."""
    period_points = list(m.period.index_set())
    if not period_points:
        return

    first_point = period_points[0]
    last_point = period_points[-1]

    @m.Constraint()
    def match_train_1_at_start_and_end(blk):
        return (
            blk.period[first_point].reverse_osmosis.ro_skid[1].op_mode
            == blk.period[last_point].reverse_osmosis.ro_skid[1].op_mode
        )


def main(season, flex_type, num_flexible_trains=4):
    season_map = {
        "summer": "price_signals/summer_week.csv",
        "winter": "price_signals/winter_week.csv",
    }
    season_key = season.lower()
    if season_key not in season_map:
        raise ValueError(
            f"Invalid season '{season}'. Valid options are: {sorted(season_map)}"
        )

    flex_type_key = flex_type.lower()
    valid_flex_types = {"rr", "flow", "both", "no_flex"}
    if flex_type_key not in valid_flex_types:
        raise ValueError(
            "Invalid flex_type "
            f"'{flex_type}'. Valid options are: {sorted(valid_flex_types)}"
        )

    selected_price_signal_stem = Path(season_map[season_key]).stem
    output_suffix = f"{season_key}_{flex_type_key}_no_brine_start_up_cost"
    if selected_price_signal_stem.upper().endswith("RTP"):
        output_suffix = f"{output_suffix}_RTP"
    if selected_price_signal_stem.upper().endswith("TOU_8"):
        output_suffix = f"{output_suffix}_TOU_8"
    if selected_price_signal_stem.upper().endswith("CPP"):
        output_suffix = f"{output_suffix}_CPP"
    if selected_price_signal_stem.upper().endswith("DR"):
        output_suffix = f"{output_suffix}_DR"

    # Get the directory where this script is located
    script_dir = Path(__file__).parent
    # Load price data
    price_data = pd.read_csv(script_dir / season_map[season_key])
    price_data["Energy Rate"] = (
        price_data["electric_energy_on_peak"]
        + price_data["electric_energy_mid_peak"]
        + price_data["electric_energy_off_peak"]
        + price_data["electric_energy_super_off_peak"]
    )
    price_data["Fixed Demand Rate"] = price_data["electric_demand_fixed"]
    price_data["Var Demand Rate"] = price_data["electric_demand_peak"]
    price_data["Customer Cost"] = price_data["electric_customer_fixed_charge"]
    price_data["Demand_Response_Price"] = price_data["electric_demand_response_price"]
    price_data["Emissions Intensity"] = 0
    peak_hours = price_data["Var Demand Rate"].to_numpy() != 0

    # Load PV data
    pv_kW = price_data["solar_output_kW"]
    pv_capacity = max(pv_kW)
    pv_capacity_factors = pv_kW / pv_capacity

    m = PriceTakerModel()
    # Find start and end datetimes and time step  from the price data
    price_datetimes = pd.to_datetime(price_data["DateTime"])
    data_start = price_datetimes.iloc[0]
    data_next_time = price_datetimes.iloc[1]
    timestep_hours = (data_next_time - data_start).total_seconds() / 3600
    start_date = data_start.strftime("%Y-%m-%d %H:%M:%S")
    end_date = price_datetimes.iloc[-1].strftime("%Y-%m-%d %H:%M:%S")

    # Instantiate an object containing the model parameters
    m.params = FlexDesalParams(
        start_date=start_date,
        end_date=end_date,
        annual_production_AF=12000,
        # * 1.0863,  # This is to compare against the October plant data
        # * 1.2602,  # This is to compare against the Aug plant data
        timestep_hours=timestep_hours,
        include_onsite_solar=False,
        onsite_capacity=pv_capacity,
        nonworking_hours=list(range(0, 8))
        + list(
            range(18, 24)
        ),  # 6pm-8am are nonworking hours (assuming time index starts at 0 for 12am-1am)
        # rainy_days=1,  # This will reduce the maxumim value for annual_production AF
        CAPEX_yr=6498300,  # For WRD, this assumes a 30 yr lifetime
        include_demand_response=True,
        max_daily_shutdowns=1,  # I'd like to change to one a day
    )
    m.baseline_power = 1102  # kW
    m.params.intake.update(
        {
            "energy_intensity": 0,
            "nominal_flowrate": 2500,
            "feed_cost": 0.16,
            "chemical_cost": 0.0332,
        }
    )  # m3/hr

    m.params.wrd_uf.update(
        {
            "minimum_downtime": 2,
            "startup_delay": 2,
            "minimum_flowrate": 344,  # m3/hr
            "nominal_flowrate": 900,
            "maximum_flowrate": 989,
            "surrogate_type": "quadratic_energy_intensity",
            "surrogate_a": 2.71e-1,
            "surrogate_b": -3.32e-4,
            "surrogate_c": 2.39e-7,
            "nominal_recovery": 0.96,
            "num_uf_pumps": 3,
        }
    )

    m.params.wrd_ro.update(
        {
            "startup_delay": 2,  # hours
            "minimum_downtime": 2,  # hours
            "minimum_flowrate": 520,  # m3/hr
            "nominal_flowrate": 602,
            "maximum_flowrate": 635,
            "allow_variable_recovery": flex_type_key not in {"flow", "no_flex"},
            "surrogate_type": "PySMO_polyfit",
            "surrogate_file": script_dir / "ro_SEC_poly_fit_order_1.json",
            "minimum_recovery": 0.88,
            "nominal_recovery": 0.925,
            "maximum_recovery": 0.925,
            "num_ro_skids": 4,
            "replacement_types": ["membranes", "motors"],
            "replacement_costs": [
                500 * 4 * (72 + 30 + 15),
                125000 * 4,
            ],  # $ per replacement
            "replacement_lifetimes": [5, 17.5],  # years
            "replacement_max_flex_penalty": [
                0.1,
                0.1,
            ],  # Reduction in lifetime if shutdowns occur twice a day
        }
    )

    m.params.posttreatment.update(
        {
            "energy_intensity": 0.101,
            "chemical_cost": 0.0310,
        }
    )  # kWh/m3 #$/m3

    m.params.brinedischarge.update({"brine_cost": 0.43, "energy_intensity": 0})

    # Append LMP data to the model
    m.append_lmp_data(lmp_data=price_data["Energy Rate"])

    m.build_multiperiod_model(
        flowsheet_func=fs.build_desal_flowsheet,
        flowsheet_options={"params": m.params},
    )

    _restrict_flexible_trains(m, num_flexible_trains=num_flexible_trains)

    _begin_and_end_constraint(m)

    # if season_key == "summer":
    #     _fix_operations_for_first_four_days(m, peak_hours=peak_hours)

    # Update the time-varying parameters other than the LMP, such as
    # demand costs and emissions intensity. LMP value is updated by default
    m.update_operation_params(
        {
            "fixed_demand_rate": price_data["Fixed Demand Rate"],
            "variable_demand_rate": price_data["Var Demand Rate"],
            "emissions_intensity": price_data["Emissions Intensity"],
            "customer_cost": price_data["Customer Cost"],
            "demand_response_price": price_data["Demand_Response_Price"],
        }
    )
    if m.params.include_onsite_solar:
        m.update_operation_params(
            {"power_generation.capacity_factor": pv_capacity_factors}
        )

    # Add demand cost and fixed cost calculation constraints
    fs.add_demand_and_fixed_costs(m)

    # Add the startup delay constraints
    fs.add_delayed_startup_constraints(m)
    fs.add_delayed_shutdown_constraints(m)
    # fs.repeat_weekdays(m)

    m.total_water_production = pyo.Expression(
        expr=m.params.timestep_hours
        * sum(m.period[:, :].posttreatment.product_flowrate)
    )
    m.total_energy_cost = pyo.Expression(expr=sum(m.period[:, :].energy_cost))

    # Demand costs are automatically normalized by number of months. So for a sample week, it multiplies by 7/31.
    m.total_demand_cost = pyo.Expression(
        expr=m.fixed_demand_cost + m.variable_demand_cost
    )
    m.total_customer_cost = pyo.Expression(
        expr=sum(m.period[:, :].customer_cost) * m.params.num_months
    )

    fs.add_flow_costs(m)  # Flow costs = Feed, Brine, and Chemicals
    fs.add_useful_expressions(m)
    # This adds the total_demand_response_revenue, which only represents one of the available SCE DR options.

    m.total_op_cost = pyo.Expression(
        expr=m.total_energy_cost
        + m.total_demand_cost
        + m.total_customer_cost
        - m.total_demand_response_revenue
        + m.total_feed_cost
        + m.total_brine_cost
        + m.total_chemical_cost
    )
    # add CAPEX as a fixed cost to calculate LCOW
    m.fixed_cost = pyo.Expression(expr=m.params.CAPEX_yr * m.params.num_months / 12)
    m.total_cost = pyo.Expression(expr=m.total_op_cost + m.fixed_cost)

    m.LCOW = pyo.Expression(expr=m.total_cost / m.total_water_production)  # $/m3

    fs.constrain_water_production(m)

    # If water recovery is static, it must be fixed
    if not m.params.wrd_ro.allow_variable_recovery:
        utils.wrd_fix_ro_recovery(
            m,
            ro_recovery=m.params.wrd_ro.nominal_recovery,
        )
    # Always want to fix the UF recovery
    utils.wrd_fix_uf_recovery(
        m,
        uf_recovery=m.params.wrd_uf.nominal_recovery,
    )

    if flex_type_key == "rr":
        _fix_nominal_flowrates(m)

    # Could cause feasibility issues b/c this is a slack variable essentially.
    # m.fix_operation_var("reverse_osmosis.leftover_flow", 0)

    # Flowrates not fixed, but shouldn't randomly fluctuate either.
    fs.add_flow_changes_penalty_binary(m)

    # restricts number of shutdowns per 24 hours period, mainly to reduce solution space
    # fs.add_maximum_shutdowns(m)

    # fs.add_working_hours_constraint(m)

    fs.add_rain_shutdowns(m)

    # This does not include the replacement costs atm because they don't drive the optimization. Also I removed the flexibility penalty
    m.obj = pyo.Objective(
        expr=1e-4
        * (
            m.total_energy_cost
            + m.total_demand_cost
            + m.total_customer_cost
            - m.total_demand_response_revenue
            + m.total_feed_cost
            + m.total_brine_cost
            + m.total_chemical_cost
            + m.flow_changes_penalty
        ),
        sense=pyo.minimize,
    )

    # Only to find the baseline power for this water production
    if flex_type_key == "no_flex" and num_flexible_trains == 0:
        m.enforce_steady_state = pyo.Constraint(expr=m.flow_changes_penalty == 0)

    ##### ADDING FOR TESTING ####
    # _shutdown_times = {42, 43, 44, 45, 46}
    # m.enforce_one_plant_shutdown = pyo.Constraint(
    #     expr=sum(
    #         m.period[1, t].reverse_osmosis.ro_skid[1].op_mode for t in _shutdown_times
    #     )
    #     == 0
    # )
    # m.enforce_one_plant_on = pyo.ConstraintList()
    # for _d, _t in m.period.index_set():
    #     if _d == 1 and _t not in _shutdown_times:
    #         m.enforce_one_plant_on.add(
    #             m.period[1, _t].reverse_osmosis.ro_skid[1].op_mode == 1
    #         )
    #### END TESTING CONSTRAINTS ####

    print(degrees_of_freedom(m))

    # dt = DiagnosticsToolbox(m)
    # dt.report_structural_issues()

    # IPOPT
    solver = get_solver()

    # mip_gap = 0.01
    # solver = pyo.SolverFactory("gurobi_direct_minlp")
    # solver.options["MIPGap"] = mip_gap  # 1.0 %
    # solver.options["MIPGapAbs"] = (
    #     0.1  # $1,000 (b/c objective function is scaled down by 1e-4)
    # )
    # solver.options["MIPFocus"] = 1
    results = solver.solve(m, tee=True)

    print(f"m.flow_changes_penalty(): {m.flow_changes_penalty()}")
    print(f"Total operational cost: {m.total_op_cost():.2f}")

    pyo.assert_optimal_termination(results)

    # Baseline power is a function of the target water production, but needs to be calculated by running this model!
    # The OPEX value does not include the replacement costs... so I guess they aren't being included in the LVOF
    if season_key == "winter":
        if selected_price_signal_stem.upper().endswith("TOU_8"):
            baseline_OPEX = 115031
        elif selected_price_signal_stem.upper().endswith("RTP"):
            baseline_OPEX = 116862
        else:
            baseline_OPEX = 111145  # $
    else:
        if selected_price_signal_stem.upper().endswith("TOU_8"):
            baseline_OPEX = 124771
        elif selected_price_signal_stem.upper().endswith("RTP"):
            baseline_OPEX = 187020
        elif selected_price_signal_stem.upper().endswith("CPP"):
            baseline_OPEX = 123683
        else:
            baseline_OPEX = 120098  # $

    fs.calculate_replacement_costs(m)
    fs.calculate_flexibility_metrics(
        m,
        baseline_power=value(
            m.baseline_power
        ),  # kW, from the baseline with steady production and 12000 AF/yr water production
        baseline_OPEX=baseline_OPEX,
    )

    design_var_values = m.get_design_var_values()
    filtered_design_var_values = {
        k: v
        for k, v in design_var_values.items()
        if "flow_change" not in k and "flow_changed" not in k and "reduction" not in k
    }
    print(filtered_design_var_values)

    # Write optimal values of all operational variables to a csv file
    output_csv = script_dir / f"wrd_result_{output_suffix}.csv"
    m.get_operation_var_values().to_csv(output_csv)
    print(f"Saved operation variable results to: {output_csv}")

    plot_function(
        m,
        n_time_points=len(price_data),
        output_stem=script_dir / f"wrd_pricetaker_{output_suffix}",
        peak_hours=peak_hours,
    )

    # # Plot operational variables
    # fig, axs = m.plot_operation_profile(
    #     operation_vars=[
    #         "fixed_demand_rate",
    #         "variable_demand_rate",
    #         "posttreatment.product_flowrate",
    #         "num_skids_online",
    #     ],
    # )
    # fig.savefig(script_dir / f"wrd_operation_profile_{output_suffix}.png")

    return m


if __name__ == "__main__":
    seasons = ["summer"]
    flex_types = ["both"]
    num_flex_skids = [4]

    results_rows = []

    for season in seasons:
        for flex_type in flex_types:
            for num_skids in num_flex_skids:
                m = main(
                    season=season, flex_type=flex_type, num_flexible_trains=num_skids
                )
                results_rows.append(
                    {
                        "Season": season,
                        "Flexibility Type": flex_type,
                        "Num Flexible Trains": num_skids,
                        "Total Operational Cost": m.total_op_cost(),
                        "Total Water Production (m3)": m.total_water_production(),
                        "LCOW ($/m3)": m.LCOW(),
                        "Total Energy Cost": m.total_energy_cost(),
                        "Fixed Demand Cost": m.fixed_demand_cost(),
                        "Variable Demand Cost": m.variable_demand_cost(),
                        "Total Electricity Cost": m.total_energy_cost()
                        + m.total_demand_cost(),
                        "Total Feed Cost": m.total_feed_cost(),
                        "Total Brine Cost": m.total_brine_cost(),
                        "Total Chemical Cost": m.total_chemical_cost(),
                        "Total Replacement Cost": m.total_replacement_cost(),
                        "Total Demand Response Revenue": m.total_demand_response_revenue(),
                        "Total Cost": m.total_cost(),
                        "Maximum Power": m.maximum_power(),
                        "Discharge Energy Capacity": m.discharge_energy_capacity(),
                        "Discharge Power Capacity": m.discharge_power_capacity(),
                        "LVOF": m.LVOF(),
                        "Charge Energy Capacity": m.charge_energy_capacity(),
                        "Charge Power Capacity": m.charge_power_capacity(),
                    }
                )

    results_df = pd.DataFrame(results_rows)
    script_dir = Path(__file__).parent
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    results_csv = script_dir / f"wrd_pricetaker_summary_results_{timestamp}.csv"
    results_df.to_csv(results_csv, index=False)
    print(f"Saved summary results to: {results_csv}")
