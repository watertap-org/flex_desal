import pandas as pd
import numpy as np
import os
from pathlib import Path

from pyomo.environ import (
    Var,
    value,
    ConcreteModel,
)

from idaes.core.surrogate.pysmo_surrogate import (
    PysmoRBFTrainer,
    PysmoPolyTrainer,
    PysmoSurrogate,
)

from idaes.core.util.scaling import calculate_variable_from_constraint

from idaes.core.surrogate.surrogate_block import SurrogateBlock
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from watertap.core.solvers import get_solver


# Helper Functions for pump curves that define plotting boundaries.
# Specific to scaled flow (1e-3) and head (1e-2) units.

__all__ = [
    "create_test_pairs",
    "filter_pump_test_points",
    "pump_param_sweep",
]


def head_limit(flow, pump_type):
    """Helper function to define minimum head for given flowrate."""
    # Coefficients from pump curve data
    # Flow must be units of scaled gpm. Head is in scaled ft.
    if pump_type == "RO_feed":
        b_0 = 3.127555
        b_1 = -0.09634
        b_2 = 0.051555
        b_3 = -0.026339
    elif pump_type == "RO_IS":
        b_0 = 1.02801
        b_1 = -0.21038
        b_2 = 0.284634
        b_3 = -0.253384
    elif pump_type == "UF":
        b_0 = 3.199441
        b_1 = -0.252764
        b_2 = 0.060008
        b_3 = -0.015999
    head = b_0 + b_1 * (flow) + b_2 * (flow) ** 2 + b_3 * (flow) ** 3
    return head


def MCSF(flow, pump_type):
    """Helper function to define max head for given flowrate along the minimum continuous stable flow."""
    # Coefficients from pump curve data
    # Flow must be in units of scaled gpm. Head is in scaled ft.
    if pump_type == "RO_feed":
        c_0 = 0.838
        c_1 = -1.869681
        c_2 = 2.477258
        c_3 = 0
    elif pump_type == "RO_IS":
        c_0 = 0.185777
        c_1 = -1.584001
        c_2 = 9.184983
        c_3 = 0
    elif pump_type == "UF":
        c_0 = -0.124779
        c_1 = 0.390064
        c_2 = 3.300655
        c_3 = 0
    head = c_0 + c_1 * (flow) + c_2 * (flow) ** 2 + c_3 * (flow) ** 3
    return head


def max_flows(flow, pump_type):
    """Helper function to define maximum flow for given flowrate along the maximum pump capacity across different speeds."""
    # Flow must be in units of scaled gpm. Head is in scaled ft.
    if pump_type == "RO_feed":
        d_0 = 0.19338
        d_1 = -0.21424
        d_2 = 0.23183
        d_3 = -0.00942
    elif pump_type == "RO_IS":
        # No multispeed head curve available
        # I think using affinity laws to determine the max flow for this pump
        d_0 = 0
        d_1 = 0
        d_2 = 0.387369
        d_3 = 0
    elif pump_type == "UF":
        # No multispeed head curve available
        # I think using affinity laws to determine the max flow for this pump
        d_0 = 0
        d_1 = 0
        d_2 = 0.04354
        d_3 = 0
    head = d_0 + d_1 * (flow) + d_2 * (flow) ** 2 + d_3 * (flow) ** 3
    return head


# User select Types


def create_surrogate(pump_type="UF", fittype="rbf"):
    # load data proudced by pump_param_sweep.py
    script_dir = Path(__file__).parent
    data_dir = script_dir / "surrogate_data"
    filename = f"{pump_type}_pump_aff_law_data_for_surrogate.csv"
    pump_data = pd.read_csv(data_dir / filename)
    input_labels = ["Flow (gpm)", "Head (ft)"]
    output_labels = ["total_efficiency"]
    input_data = pump_data[input_labels]
    output_data = pump_data[output_labels]

    # Scale Data
    Data_scaled = pump_data.copy()
    Data_scaled[output_labels[0]] = pump_data[output_labels[0]].mul(
        1e-2
    )  # Efficiency (%)
    Data_scaled[input_labels[0]] = pump_data[input_labels[0]].mul(1e-3)  # Flow (gpm)
    Data_scaled[input_labels[1]] = pump_data[input_labels[1]].mul(1e-2)  # Head (ft)
    min_flow = min(Data_scaled[input_labels[0]])
    max_flow = max(Data_scaled[input_labels[0]])
    min_head = min(Data_scaled[input_labels[1]])
    max_head = max(Data_scaled[input_labels[1]])
    input_bounds = {
        input_labels[0]: (min_flow, max_flow),
        input_labels[1]: (min_head, max_head),
    }

    # Create Surrogate Type and trainer
    # Create the trainer
    if fittype == "poly":
        trainer = PysmoPolyTrainer(
            input_labels=input_labels,
            output_labels=output_labels,
            training_dataframe=Data_scaled,  # NOT spliting data for training and validation because there's not that much data to begin with.
        )
        trainer.config.maximum_polynomial_order = 3

    if fittype == "rbf":
        trainer = PysmoRBFTrainer(
            input_labels=input_labels,
            output_labels=output_labels,
            training_dataframe=Data_scaled,
            basis_function="cubic",
        )

    # Train Data
    trained_surr = trainer.train_surrogate()

    # Create Surrogate Model
    Surrogate = PysmoSurrogate(
        trained_surr,
        input_labels,
        output_labels,
        input_bounds,
    )

    # Display Results of surrogate
    m = ConcreteModel()
    m.flowrate = Var()
    m.head = Var()
    m.eff = Var()

    m.surrogate_blk = SurrogateBlock(concrete=True)
    m.surrogate = Surrogate
    m.surrogate_blk.build_model(
        m.surrogate,
        input_vars=[m.flowrate, m.head],
        output_vars=[m.eff],
    )
    m.surrogate_blk.pysmo_constraint.display()  # display()

    # Get the absolute path of the current script
    script_dir = Path(__file__).parent
    data_dir = script_dir / "surrogate_data"
    surr_name = f"{pump_type}_total_efficiency_rbf_fit.json"
    Surrogate.save_to_file(os.path.join(data_dir, surr_name), overwrite=True)

    # Plotting Surrogate
    num_points = 50
    x_vals = np.linspace(min_flow, max_flow, num=num_points)
    y_vals = np.linspace(min_head, max_head, num=num_points)
    z_vals = np.zeros((num_points, num_points))
    for i in range(num_points):
        for j in range(num_points):
            # Filter out infeasible points
            if (
                y_vals[j] > head_limit(x_vals[i], pump_type)
                or y_vals[j] > MCSF(x_vals[i], pump_type)
                or y_vals[j] < max_flows(x_vals[i], pump_type)
            ):
                z_vals[i, j] = np.nan
                continue
            m.flowrate.fix(x_vals[i])
            m.head.fix(y_vals[j])
            calculate_variable_from_constraint(
                m.eff, m.surrogate_blk.pysmo_constraint["total_efficiency"]
            )
            z_vals[i, j] = value(m.eff)
    # Rescale the outputs
    flows = x_vals * 1e3  # gpm
    heads = y_vals * 1e2  # ft
    efficiencies = z_vals * 1e2  # %

    # Create Meshgrid to plot contour
    X, Y = np.meshgrid(flows, heads)
    Z = efficiencies.T  # Transpose Z to match X and Y dimensions

    fig1 = plt.contourf(
        X,
        Y,
        Z,
        levels=25,
        cmap="viridis",
    )

    plt.colorbar(label="Total Efficiency (%)")
    # Add countour lines for specific efficiencies in the orginial plot
    if pump_type == "RO_feed":
        plt.xlim(0, 4500)
        plt.ylim(0, 400)

    elif pump_type == "RO_IS":
        plt.xlim(0, 1600)
        plt.ylim(0, 150)

    elif pump_type == "UF":
        plt.xlim(0, 6000)
        plt.ylim(0, 400)

    plt.xlabel("Flow (gpm)")
    plt.ylabel("Head (ft)")
    plt.title(f"{pump_type.replace('_', ' ').title()} Total Efficiency Contour Plot")

    # Plotting the training data points to reveal if there are poorly sampled regions
    X_data = pump_data["Flow (gpm)"]
    Y_data = pump_data["Head (ft)"]
    # Z_data = pump_data["total_efficiency"]
    plt.scatter(X_data, Y_data, color="red", s=10, label="Training Data")
    plt.legend()
    plt.show()


if __name__ == "main":
    pump_type = "UF"
    fittype = "rbf"
    create_surrogate(pump_type, fittype)
