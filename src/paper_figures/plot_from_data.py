## This script is meant to recreate the plotting results from using the csv files.
## The intention is for it to be edited to for customization of the plot formating.

from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def op_plot_from_data(filename, data_type="optimization_results"):
    script_dir = Path(__file__).resolve().parent
    file_path = Path(data_type + "/" + filename)
    if not file_path.is_absolute():
        file_path = script_dir / file_path

    df = pd.read_csv(file_path)

    # Column names depend on data_type
    if data_type == "optimization_results":
        column_names = {
            "total_energy": "net_power_consumption",
            "ro1_energy": "reverse_osmosis.ro_skid[1].power_consumption",
            "ro2_energy": "reverse_osmosis.ro_skid[2].power_consumption",
            "ro3_energy": "reverse_osmosis.ro_skid[3].power_consumption",
            "ro4_energy": "reverse_osmosis.ro_skid[4].power_consumption",
            "uf1_energy": "pretreatment.uf_pumps[1].power_consumption",
            "uf2_energy": "pretreatment.uf_pumps[2].power_consumption",
            "uf3_energy": "pretreatment.uf_pumps[3].power_consumption",
            "elec_price": "LMP",
            "prod": "posttreatment.product_flowrate",
            "train_1_flows": "reverse_osmosis.ro_skid[1].product_flowrate",
            "train_2_flows": "reverse_osmosis.ro_skid[2].product_flowrate",
            "train_3_flows": "reverse_osmosis.ro_skid[3].product_flowrate",
            "train_4_flows": "reverse_osmosis.ro_skid[4].product_flowrate",
            "peak_hours": "peak_hour",
            "demand_response_revenue": "demand_response_revenue",
            "working_hours": "working_hour",
        }
    elif data_type == "plant_data":
        # For plant data, flowrates are stored as percentages and power is pre-aggregated
        # max_perm_flowrate = 569.0  # m3/hr # I believe this value works for the August week data???
        max_perm_flowrate = 588.2  # m3/hr # For October week data.
        column_names = {
            "total_energy": None,  # Will be computed
            "ro1_energy": "RO_train_1_kW",
            "ro2_energy": "RO_train_2_kW",
            "ro3_energy": "RO_train_3_kW",
            "ro4_energy": "RO_train_4_kW",
            "uf1_energy": "UFFeedPumps_1_kW",
            "uf2_energy": "UFFeedPumps_2_kW",
            "uf3_energy": "UFFeedPumps_3_kW",
            "elec_price": "LMP",
            "prod": None,  # Will be computed from flowrate percentages
            "train_1_flows": ("train_1_flow_pct", max_perm_flowrate),
            "train_2_flows": ("train_2_flow_pct", max_perm_flowrate),
            "train_3_flows": ("train_3_flow_pct", max_perm_flowrate),
            "train_4_flows": ("train_4_flow_pct", max_perm_flowrate),
            "peak_hours": "peak_hour",
            "demand_response_revenue": "demand_response_revenue",
            "working_hours": "working_hour",
        }
        # Compute total energy as sum of all power columns
        power_cols = [
            "RO_train_1_kW",
            "RO_train_2_kW",
            "RO_train_3_kW",
            "RO_train_4_kW",
            "UFFeedPumps_1_kW",
            "UFFeedPumps_2_kW",
            "UFFeedPumps_3_kW",
            "UVAOP_kW",
        ]
        df["_computed_total_energy"] = df[power_cols].sum(axis=1)
        # Compute total flowrate from percentages
        df["_computed_prod"] = (
            df[
                [
                    "train_1_flow_pct",
                    "train_2_flow_pct",
                    "train_3_flow_pct",
                    "train_4_flow_pct",
                ]
            ].sum(axis=1)
            * max_perm_flowrate
            / 100.0
        )
    else:
        raise ValueError(
            f"Unsupported data_type '{data_type}'. Valid options are: "
            "'optimization_results' and 'plant_data'."
        )

    def get_series(series_name, required=True):
        spec = column_names[series_name]

        # Handle computed columns for plant_data
        if spec is None:
            if series_name == "total_energy" and data_type == "plant_data":
                return df["_computed_total_energy"].to_numpy()
            elif series_name == "prod" and data_type == "plant_data":
                return df["_computed_prod"].to_numpy()
            elif required:
                raise KeyError(
                    f"Column for '{series_name}' has not been configured or is not "
                    f"available for data_type='{data_type}'."
                )
            return None

        # Handle percentage-to-flowrate conversion for plant_data
        if isinstance(spec, tuple):
            pct_col, max_flow = spec
            if pct_col not in df.columns:
                raise KeyError(
                    f"Column '{pct_col}' for '{series_name}' was not found in {file_path}."
                )
            return df[pct_col].to_numpy() * max_flow / 100.0

        # Handle regular column lookups
        if spec not in df.columns:
            if required:
                raise KeyError(
                    f"Column '{spec}' for '{series_name}' was not found in {file_path}."
                )
            return None
        return df[spec].to_numpy()

    total_energy = get_series("total_energy")
    n_time_points = len(total_energy)
    output_name = file_path.stem
    output_stem = (
        output_name.split("wrd_result_", 1)[1]
        if "wrd_result_" in output_name
        else output_name
    )

    working_hours = get_series("working_hours", required=False)
    peak_hours_data = get_series("peak_hours", required=False)
    peak_hours = (
        None
        if peak_hours_data is None or np.all(peak_hours_data == 0)
        else peak_hours_data.astype(bool)
    )
    if peak_hours is None:
        print("Peak hours are missing or all zero; skipping peak-hour shading.")

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
                    zorder=-2,
                    hatch="///",
                    label=span_label,
                )
                ax_trains.axvspan(
                    i,
                    i + 1,
                    color="grey",
                    alpha=0.2,
                    linewidth=0,
                    zorder=-2,
                    hatch="///",
                    label=span_label,
                )
                peak_legend_added = True

    demand_response_event = get_series("demand_response_revenue", required=False) != 0
    if any(demand_response_event):
        DR_legend_added = False
        for i, is_DR in enumerate(demand_response_event):
            if is_DR:
                # Shade full hourly intervals where Demand Response Event occurs
                span_label = "DR Event" if not DR_legend_added else None
                ax_energy.axvspan(
                    i,
                    i + 1,
                    color="red",
                    alpha=0.2,
                    linewidth=0,
                    zorder=-1,
                    hatch="///",
                    label=span_label,
                )
                ax_trains.axvspan(
                    i,
                    i + 1,
                    color="red",
                    alpha=0.2,
                    linewidth=0,
                    zorder=-1,
                    hatch="///",
                    label=span_label,
                )
                DR_legend_added = True

    if working_hours is not None:
        working_legend_added = False
        for i, is_working in enumerate(working_hours):
            if is_working:
                # Shade full hourly intervals where the plant is operating.
                span_label = "Working Hours" if not working_legend_added else None
                ax_energy.axvspan(
                    i,
                    i + 1,
                    color="lightblue",
                    alpha=0.2,
                    linewidth=0,
                    zorder=-3,
                    hatch="\\\\\\",
                    label=span_label,
                )
                ax_trains.axvspan(
                    i,
                    i + 1,
                    color="lightblue",
                    alpha=0.2,
                    linewidth=0,
                    zorder=-3,
                    hatch="\\\\\\",
                    label=span_label,
                )
                working_legend_added = True

    # First subplot: Stacked energy consumption by major equipment
    ro1_energy = get_series("ro1_energy")
    ro2_energy = get_series("ro2_energy")
    ro3_energy = get_series("ro3_energy")
    ro4_energy = get_series("ro4_energy")

    uf1_energy = get_series("uf1_energy")
    uf2_energy = get_series("uf2_energy")
    uf3_energy = get_series("uf3_energy")

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
    elec_price = get_series("elec_price", required=False)
    if elec_price is not None and not np.all(np.isnan(elec_price)):
        ax_price.plot(
            time + 0.5,
            elec_price * 100,
            label="Elec. Price",
            color="orange",
            linewidth=2,
        )
        ax_price.set_ylabel("¢/kWh", fontsize=12)
        ax_price.set_ylim(0, np.nanmax(elec_price * 100) + 3)
    else:
        ax_price.set_ylabel("¢/kWh", fontsize=12)
        ax_price.set_ylim(0, 1)
        ax_price.set_visible(False)

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
    prod = get_series("prod")
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
    train_1_flows = get_series("train_1_flows")
    train_2_flows = get_series("train_2_flows")
    train_3_flows = get_series("train_3_flows")
    train_4_flows = get_series("train_4_flows")

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
    print(
        f"total water produced: {sum(sum([train_1_flows,train_2_flows,train_3_flows,train_4_flows]))} m3"
    )
    print(
        f"total water produced: {sum(sum([train_1_flows,train_2_flows,train_3_flows,train_4_flows]))/1233} AF"
    )

    fig.tight_layout()
    output_dir = script_dir / "paper_figs"
    output_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_dir / f"{output_name}.png", dpi=600)
    plt.show()


# Optimization Results
# filename = "wrd_result_summer_both_4_flexible_trains.csv"

# Plant Data
filename = "wrd_result_summer_full_flex_RTP.csv"

op_plot_from_data(filename, data_type="optimization_results")
