"""Standalone plotting helpers for comparing model-fit energy to plant data."""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path


TOTAL_PLANT_PRODUCTION_CAPACITY = 53150 / 24  # m3/hr
MAX_TRAIN_FLOW = TOTAL_PLANT_PRODUCTION_CAPACITY / 4  # m3/hr
TRAIN_IDS = (1, 2, 3, 4)

DEFAULT_ACTUAL_ENERGY_CSV = Path(__file__).with_name("Aug_21_kW_hourly.csv")
DEFAULT_TRAIN_SCHEDULE_CSV = Path(__file__).with_name("real_operation_Aug_2021.csv")


def calculate_sim_energy_profile(train_schedule):
    """Calculate simulated energy from the fitted linear rules in the model."""

    if not isinstance(train_schedule, pd.DataFrame):
        schedule_path = Path(train_schedule)
        if not schedule_path.is_absolute():
            schedule_path = Path(__file__).resolve().parent / schedule_path
        train_schedule = pd.read_csv(schedule_path)

    required_cols = [
        *[f"train_{train_id}_on" for train_id in TRAIN_IDS],
        *[f"train_{train_id}_flow_pct" for train_id in TRAIN_IDS],
    ]
    missing_cols = [col for col in required_cols if col not in train_schedule.columns]
    if missing_cols:
        raise ValueError(f"Schedule is missing required columns: {missing_cols}")

    schedule = train_schedule.copy()
    on_data = schedule[[f"train_{train_id}_on" for train_id in TRAIN_IDS]].astype(float)
    pct_data = schedule[[f"train_{train_id}_flow_pct" for train_id in TRAIN_IDS]].astype(float)

    train_flows = MAX_TRAIN_FLOW * pct_data / 100.0
    print("Train flows (m3/hr):", train_flows.head())
    total_water_production = train_flows.sum(axis=1)
    uf_on = (on_data.sum(axis=1) > 0).astype(float)

    sim_energy_profile = (
        0.6343 * train_flows[f"train_1_flow_pct"].to_numpy()
        - 139.4 * on_data[f"train_1_on"].to_numpy()
        + 0.6343 * train_flows[f"train_2_flow_pct"].to_numpy()
        - 139.4 * on_data[f"train_2_on"].to_numpy()
        + 0.6343 * train_flows[f"train_3_flow_pct"].to_numpy()
        - 139.4 * on_data[f"train_3_on"].to_numpy()
        + 0.6343 * train_flows[f"train_4_flow_pct"].to_numpy()
        - 139.4 * on_data[f"train_4_on"].to_numpy()
        # So we don't have the flowrate to each individual UF pump, so this value is based on total flowrate to all UF pumps, using 4, 3, 2 RO trains. It's a bit different from what is encoded in the pricetaker model
        + 0.199 * total_water_production.to_numpy() 
        - 27.4 * uf_on.to_numpy()
        + 0.101 * total_water_production.to_numpy() # Post treatment
    )
    print((0.6343 * train_flows[f"train_1_flow_pct"].to_numpy() - 139.4 * on_data[f"train_1_on"].to_numpy())[:3])
    print((0.199 * total_water_production.to_numpy() - 27.4 * uf_on.to_numpy())[:3])
    print((0.101 * total_water_production.to_numpy())[:3])  # Post treatment

    return sim_energy_profile.tolist()


def validation_plot(
    train_schedule=DEFAULT_TRAIN_SCHEDULE_CSV,
    actual_energy_csv=DEFAULT_ACTUAL_ENERGY_CSV,
):
    """Plot hard-coded simulated energy against actual plant data."""
    csv_path = Path(actual_energy_csv)
    if not csv_path.is_absolute():
        csv_path = Path(__file__).resolve().parent / csv_path

    sim_energy_profile = calculate_sim_energy_profile(train_schedule)
    act_energy_profile = pd.read_csv(csv_path)["total_energy_kW"].to_list()

    if len(sim_energy_profile) != len(act_energy_profile):
        raise ValueError(
            "Simulated energy and actual energy profiles must have the same length."
        )

    n_time_points = len(sim_energy_profile)
    time = np.linspace(0, n_time_points - 1, n_time_points)

    fig, ax = plt.subplots(1, 1, figsize=(12, 6))

    sim_energy_line = ax.plot(
        time + 0.5,
        sim_energy_profile,
        label="Modeled Energy Consumption (kWh)",
        color="orange",
        marker="o",
    )

    act_energy_line = ax.plot(
        time + 0.5,
        act_energy_profile,
        label="Energy Consumption Data (kWh)",
        color="blue",
        marker="s",
    )

    ax.set_ylim(0, 2500)
    ax.set_ylabel("Energy Consumption (kWh)", fontsize=16)
    ax.set_xlabel("Hours", fontsize=16)
    ax.set_title("Energy Consumption - October 2021", fontsize=18, fontweight="bold")
    ax.grid(False)
    ax.xaxis.set_major_locator(plt.MaxNLocator(24))

    ax.legend(
        handles=[sim_energy_line[0], act_energy_line[0]],
        loc="lower left",
        framealpha=1.0,
        fontsize=11,
    )

    ax.set_xlim(0, n_time_points)
    ax.tick_params(axis="both", labelsize=11)
    fig.tight_layout()
    output_path = Path(__file__).resolve().parent / f"{csv_path.stem}_WRD_model_validation_Oct_2021.png"
    fig.savefig(output_path, dpi=600)
    plt.show()


def energy_data_plot(energy_csv, energy_col="total_energy_kW"):
    """Plot measured energy data from a CSV using the validation plot style."""

    csv_path = Path(energy_csv)
    if not csv_path.is_absolute():
        csv_path = Path(__file__).resolve().parent / csv_path

    energy_df = pd.read_csv(csv_path)
    if energy_col not in energy_df.columns:
        raise ValueError(
            f"Energy CSV is missing required column '{energy_col}'. "
            f"Available columns: {list(energy_df.columns)}"
        )

    energy_profile = energy_df[energy_col].to_list()
    n_time_points = len(energy_profile)
    time = np.linspace(0, n_time_points - 1, n_time_points)

    fig, ax = plt.subplots(1, 1, figsize=(12, 6))

    energy_line = ax.plot(
        time + 0.5,
        energy_profile,
        label="Energy Consumption Data (kWh)",
        color="blue",
        marker="s",
    )

    ax.set_ylim(0, 2500)
    ax.set_ylabel("Energy Consumption (kWh)", fontsize=16)
    ax.set_xlabel("Hours", fontsize=16)
    ax.set_title("Energy Consumption Data", fontsize=18, fontweight="bold")
    ax.grid(False)
    ax.xaxis.set_major_locator(plt.MaxNLocator(24))

    ax.legend(
        handles=[energy_line[0]],
        loc="lower left",
        framealpha=1.0,
        fontsize=11,
    )

    ax.set_xlim(0, n_time_points)
    ax.tick_params(axis="both", labelsize=11)
    fig.tight_layout()
    output_path = Path(__file__).resolve().parent / f"{csv_path.stem}_energy_data_plot.png"
    fig.savefig(output_path, dpi=600)
    plt.show()


def calc_energy_costs(energy_csv, season='summer', energy_col="total_energy_kW"):
    """Calculate total energy costs for a given CSV of energy data."""

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

        total_hours = int(n)
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
        # Delivery Pricing $/kWh
        on_peak_del = 0
        mid_peak_del = 0.01927
        off_peak_del = 0.01811
        super_off_peak_del = 0.01745

        # Generation Pricing $/kWh
        on_peak_gen = 0
        mid_peak_gen = 0.09639
        off_peak_gen = 0.09695
        super_off_peak_gen = 0.05329

        day_elec_price = np.ones(24)

        # off peak 12 AM - 8 AM
        day_elec_price[0:8] = off_peak_del + off_peak_gen
        # super off peak 8 AM - 4 PM
        day_elec_price[8:16] = super_off_peak_del + super_off_peak_gen
        # mid peak 4 PM - 9 PM
        day_elec_price[16:21] = mid_peak_del + mid_peak_gen
        # off peak 9 PM - 12 AM
        day_elec_price[21:24] = off_peak_del + off_peak_gen

        total_hours = int(n)
        if total_hours <= 0:
            return np.array([]), []

        elec_price = np.tile(day_elec_price, (total_hours + 23) // 24)[:total_hours]

        # Absolute hourly indices that fall in mid-peak window (16:00-21:00).
        peak_hours = []
        for h in range(total_hours):
            day_of_week = (h // 24) % 7
            hour_of_day = h % 24
            if 16 <= hour_of_day <= 20:
                peak_hours.append(h)

        return elec_price, peak_hours


    csv_path = Path(energy_csv)
    if not csv_path.is_absolute():
        csv_path = Path(__file__).resolve().parent / csv_path
    energy_df = pd.read_csv(csv_path)
    energy_profile = energy_df[energy_col].to_numpy()
    n = len(energy_profile)
    print(n)
    if season == 'summer':
        elec_price, peak_hours = build_elec_price_summer(n)
        fixed_demand_price = 19.94   # $/kW
        variable_demand_price = 36.78  # $/kW
    elif season == 'winter':
        elec_price, peak_hours = build_elec_price_winter(n)
        fixed_demand_price = 19.62   # $/kW
        variable_demand_price = 10.54  # $/kW

    month_to_week_factor = 7/30  # Average number of weeks in a month

    energy_cost = float((energy_profile * elec_price).sum())
    fixed_demand_cost = float(energy_profile.max() * fixed_demand_price) * month_to_week_factor
    variable_demand_cost = float(energy_profile[peak_hours].max() * variable_demand_price) * month_to_week_factor

    total_cost = energy_cost + fixed_demand_cost + variable_demand_cost
    total_energy_kWh = float(energy_profile.sum())

    print(f"Total Energy Cost ($): {total_cost:.1f}")
    print(f"Energy Cost ($): {energy_cost:.1f}")
    print(f"Fixed Demand Cost ($): {fixed_demand_cost:.1f}")
    print(f"Variable On-Peak Demand Cost ($): {variable_demand_cost:.1f}")
    print(f"Total Energy Consumption (kWh): {total_energy_kWh:.1f}")

    return {
        "total_cost": total_cost,
        "energy_cost": energy_cost,
        "fixed_demand_cost": fixed_demand_cost,
        "variable_demand_cost": variable_demand_cost,
        "total_energy_kWh": total_energy_kWh,
    }


if __name__ == "__main__":
    validation_plot(actual_energy_csv='Oct_21_kW_hourly_week.csv', train_schedule='real_operation_Oct_2021.csv')
    filename = "Oct_21_kW_hourly_week.csv"
    # energy_data_plot(filename)
    calc_energy_costs(filename,season='winter')

    