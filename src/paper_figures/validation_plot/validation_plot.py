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
    pct_data = schedule[
        [f"train_{train_id}_flow_pct" for train_id in TRAIN_IDS]
    ].astype(float)

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
        + 0.101 * total_water_production.to_numpy()  # Post treatment
    )
    print(
        (
            0.6343 * train_flows[f"train_1_flow_pct"].to_numpy()
            - 139.4 * on_data[f"train_1_on"].to_numpy()
        )[:3]
    )
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
    output_path = (
        Path(__file__).resolve().parent / f"{csv_path.stem}_WRD_model_validation.png"
    )
    fig.savefig(output_path, dpi=600)
    plt.show()


if __name__ == "__main__":
    validation_plot(
        actual_energy_csv="Oct_21_kW_hourly_week.csv",
        train_schedule="Oct_21_real_operation.csv",
    )
    filename = "Oct_21_kW_hourly_week.csv"
