import pandas as pd
import numpy as np
import os

from pyomo.environ import (
    Param,
    Var,
    Suffix,
    NonNegativeReals,
    value,
    units as pyunits,
    ConcreteModel,
)

from idaes.core.surrogate.pysmo_surrogate import (
    PysmoRBFTrainer,
    PysmoPolyTrainer,
    PysmoSurrogate,
)
from idaes.core.surrogate.sampling.data_utils import (
    split_training_validation,
)  # Yes, it's a random split
from idaes.core.surrogate.plotting.sm_plotter import (
    surrogate_scatter2D,
    surrogate_parity,
    surrogate_residual,
)
from idaes.core.util.scaling import calculate_variable_from_constraint

from idaes.core.surrogate.surrogate_block import SurrogateBlock
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from watertap.core.solvers import get_solver


# Helper Function, stealing fit from pump head curve, but could also make a surrogate from json data
def head_limit(flow):
    """Helper function to define minimum head for given flowrate."""
    # Coefficients from pump curve data
    # Flow must be units of m3/s. Head is in ft.
    b_0 = 374.73
    b_1 = -1347.1
    b_2 = 8953.8
    b_3 = -26539
    head = b_0 + b_1 * flow + b_2 * flow**2 + b_3 * flow**3
    return head


def MCSF(flow):
    """Helper function to define max head for given flowrate along the minimum continuous stable flow."""
    # Coefficients from pump curve data
    # Flow must be in units of gpm. Head is in ft.
    c_0 = 0.838
    c_1 = -1.869681
    c_2 = 2.477258
    c_3 = 0
    head = 100 * (
        c_0 + c_1 * (flow / 1e3) + c_2 * (flow / 1e3) ** 2 + c_3 * (flow / 1e3) ** 3
    )
    return head


def min_head_limit(flow):
    """Helper function to define minimum head for given flowrate."""
    # Coefficients from pump curve data
    # Flow must be units of gpm. Head is in ft.
    a_0 = 1.8033
    a_1 = -0.202589
    a_2 = 0.148054
    a_3 = -0.060852
    head = 100 * (a_0 + a_1 * (flow / 1e3) + a_2 * (flow / 1e3) ** 2 + a_3 * (flow / 1e3) ** 3)
    return head


# load data
pump_data = pd.read_csv(
    os.path.join(os.path.dirname(__file__), "RO_feed_pump_eff_curve_data.csv")
)
input_labels = ["Flow (gpm)", "Head (ft)"]
output_labels = ["efficiency"]
input_data = pump_data[input_labels]
output_data = pump_data[output_labels]

# Scale Data
Data_scaled = pump_data
Data_scaled[output_labels[0]] = pump_data[output_labels[0]].mul(1e-2)  # Efficiency (%)
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
fittype = "rbf"
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
        basis_function="linear",
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

# minx, maxx = m.flowrate.bounds
# miny, maxy = m.power.bounds

num_points = 25
x_vals = np.linspace(min_flow, max_flow, num=num_points)
y_vals = np.linspace(min_head, max_head, num=num_points)
z_vals = np.zeros((num_points, num_points))
solver = get_solver()
for i in range(num_points):
    for j in range(num_points):
        # Filter out infeasible points
        if (
            y_vals[j] > head_limit(x_vals[i] / 264.2 / 60 * 1e3) / 1e2
            or y_vals[j] > MCSF(x_vals[i] * 1e3) / 1e2
            or y_vals[j] < min_head_limit(x_vals[i] * 1e3) / 1e2
        ): 
            z_vals[i, j] = np.nan
            continue
        m.flowrate.fix(x_vals[i])
        m.head.fix(y_vals[j])
        calculate_variable_from_constraint(
            m.eff, m.surrogate_blk.pysmo_constraint["efficiency"]
        )
        # print(m.flowrate.value, m.head.value)
        # print(value(m.eff))
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
plt.colorbar(label="Efficiency (%)")
# Add countour lines for specific efficiencies in the orginial plot
iso_eff = plt.contour(
    X,
    Y,
    Z,
    levels=[59,68,75,80,83],
    labels=["59%","68%","75%","80%","83%"],
    colors='red',
)
plt.clabel(iso_eff, inline=True, fontsize=8, fmt="%1.0f%%", colors='black')

plt.xlabel("Flow (gpm)")
plt.ylabel("Head (ft)")
plt.title("RO Feed Pump Efficiency Contour Plot")
plt.xlim(0, 4500)
plt.ylim(0, 400)
plt.show()


# # Create 3D scatter plot
# fig2 = plt.figure(figsize=(10, 8))
# ax2 = fig2.add_subplot(111, projection="3d")

# # Flatten the data for scatter plot
# X_flat = X.flatten()
# Y_flat = Y.flatten()
# Z_flat = Z.flatten()

# # Remove NaN values for cleaner plot
# mask = ~np.isnan(Z_flat)
# X_clean = X_flat[mask]
# Y_clean = Y_flat[mask]
# Z_clean = Z_flat[mask]

# # Create scatter plot with color mapping
# scatter = ax2.scatter(
#     X_clean, Y_clean, Z_clean, c=Z_clean, cmap="viridis", s=50, alpha=0.6
# )
# fig2.colorbar(scatter, ax=ax, label="Efficiency (%)")

# # Plot actual data points
# X = pump_data["Flow (gpm)"] * 1e3
# Y = pump_data["Head (ft)"] * 1e2
# Z = pump_data["efficiency"] * 1e2
# ax2.scatter(X, Y, Z, color="red", s=20, label="Data Points")

# ax2.set_xlabel("Flow (gpm)")
# ax2.set_ylabel("Head (ft)")
# ax2.set_zlabel("Efficiency (%)")
# ax2.set_title("RO Feed Pump Efficiency 3D Scatter Plot")
# ax2.set_xlim(0, 4500)
# ax2.set_ylim(0, 400)
# plt.show()


# # Get the absolute path of the current script
# current_script_path = os.path.abspath(__file__)
# # Get the directory containing the current script
# current_directory = os.path.dirname(current_script_path)
# surr_name = f"UV_power_poly_fit_order_{trainer.config.maximum_polynomial_order}.json"
# Surrogate.save_to_file(current_directory+surr_name, overwrite=True)
