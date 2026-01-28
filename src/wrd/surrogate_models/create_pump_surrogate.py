import pandas as pd
import numpy as np
import os

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


# Helper Function, stealing fit from pump head curve, but could also make a surrogate from json data
def min_head_limit(flow,pump_type):
    """Helper function to define minimum head for given flowrate."""
    # Coefficients from pump curve data
    # Flow must be units of gpm. Head is in ft.
    # This isn't valid anymore b/ it was for a smaller pump impeller diameter,not speed
    if pump_type == 'RO_feed':
        a_0 = 1.8033
        a_1 = -0.202589
        a_2 = 0.148054
        a_3 = -0.060852
    if pump_type == 'RO_IS':
        a_0 = 0.607309
        a_1 = -0.059714
        a_2 = -0.105637
        a_3 = -0.10102
    head = a_0 + a_1 * (flow) + a_2 * (flow) ** 2 + a_3 * (flow) ** 3
    return head


def head_limit(flow,pump_type):
    """Helper function to define minimum head for given flowrate."""
    # Coefficients from pump curve data
    # Flow must be units of scaled gpm. Head is in ft.
    if pump_type == 'RO_feed':
        b_0 = 3.127555
        b_1 = -0.09634
        b_2 = 0.051555
        b_3 = -0.026339
    if pump_type == 'RO_IS':
        b_0 = 1.02801
        b_1 = -0.21038
        b_2 = 0.284634
        b_3 = -0.253384
    head = (b_0 + b_1 * (flow) + b_2 * (flow)**2 + b_3 * (flow)**3)
    return head


def MCSF(flow,pump_type):
    """Helper function to define max head for given flowrate along the minimum continuous stable flow."""
    # Coefficients from pump curve data
    # Flow must be in units of scaled gpm. Head is in scaled ft.
    if pump_type == 'RO_feed':
        c_0 = 0.838
        c_1 = -1.869681
        c_2 = 2.477258
        c_3 = 0
    if pump_type == 'RO_IS':
        c_0 = 0.185777
        c_1 = -1.584001
        c_2 = 9.184983
        c_3 = 0
    head = (
        c_0 + c_1 * (flow) + c_2 * (flow) ** 2 + c_3 * (flow) ** 3
    )
    return head


def max_flows(flow,pump_type):
    """Helper function to define maximum flow for given flowrate along the maximum pump capacity across different speeds."""
    # Flow must be in units of scaled gpm. Head is in scaled ft.
    if pump_type == 'RO_feed':
        d_0 = 0.19338
        d_1 = -0.21424
        d_2 = 0.23183
        d_3 = -0.00942
    if pump_type == 'RO_IS':
        # No multispeed head curve available
        # I think using affinity laws to determine the max flow for this pump
        d_0 = 0
        d_1 = 0
        d_2 = 0.387369 
        d_3 = 0
    head = d_0 + d_1 * (flow) + d_2 * (flow) ** 2 + d_3 * (flow) ** 3
    return head


# Select Types
pump_type = 'RO_IS'
fittype = "rbf"

# load data
# filename = f"{pump_type}_pump_eff_curve_data.csv"
filename = "RO_IS_aff_laws_surr_2.csv"
pump_data = pd.read_csv(
    os.path.join(os.path.dirname(__file__), filename)
)
input_labels = ["Flow (gpm)", "Head (ft)"]
output_labels = ["total_efficiency"]
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

# minx, maxx = m.flowrate.bounds
# miny, maxy = m.power.bounds

num_points = 50
x_vals = np.linspace(min_flow, max_flow, num=num_points)
y_vals = np.linspace(min_head, max_head, num=num_points)
z_vals = np.zeros((num_points, num_points))
solver = get_solver()
for i in range(num_points):
    for j in range(num_points):
        # Filter out infeasible points
        if (
            y_vals[j] > head_limit(x_vals[i],pump_type) 
            or y_vals[j] > MCSF(x_vals[i],pump_type)
            # or y_vals[j] < min_head_limit(x_vals[i],pump_type)
            or y_vals[j] < max_flows(x_vals[i],pump_type)
        ): 
            z_vals[i, j] = np.nan
            continue
        m.flowrate.fix(x_vals[i])
        m.head.fix(y_vals[j])
        calculate_variable_from_constraint(
            m.eff, m.surrogate_blk.pysmo_constraint["total_efficiency"]
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
plt.colorbar(label="Total Efficiency (%)")
# Add countour lines for specific efficiencies in the orginial plot
if pump_type == 'RO_feed':
    plt.xlim(0, 4500)
    plt.ylim(0, 400)
    levels=[59,68,75,80,83]
    labels=["59%","68%","75%","80%","83%"]

elif pump_type == 'RO_IS':
    plt.xlim(0, 1600)
    plt.ylim(0, 150)
    levels=[52,62,71,77,81,82]
    labels=["52%","62%","71%","77%","81%","82%"]

plt.xlabel("Flow (gpm)")
plt.ylabel("Head (ft)")
plt.title(f"{pump_type.replace('_', ' ').title()} Pump Efficiency Contour Plot")

# Test value at 80% speed point
m.flowrate.fix(2778/1e3)
m.head.fix(150.8/1e2)
calculate_variable_from_constraint(
    m.eff, m.surrogate_blk.pysmo_constraint["total_efficiency"]
)
print(m.flowrate.value, m.head.value)
print(value(m.eff))

# BELOW CODE WON'T WORK ANYMORE. WE ONLY HAVE A FEW ACTUAL DATA POINTS (100% CURVE, PLUS ONE 80% SPEED POINT)
# # Data points vs. surrogate outputs
X_data = pump_data["Flow (gpm)"] * 1e3
Y_data = pump_data["Head (ft)"] * 1e2
Z_data = pump_data["total_efficiency"] * 1e2
plt.scatter(X_data, Y_data, color="red", s=10)

# for l in levels:
#     X_points = X_data[Z_data == l]
#     Y_points = Y_data[Z_data == l]    
#     plt.scatter(X_points, Y_points, color="red", s=10, label=f"{l}%")
#     plt.text(X_points.iloc[0]+20,Y_points.iloc[0]+3, f"{l}%",color="red")
plt.show()

# Data Validation Scatter plot

# Z_modeled = pd.DataFrame()

# for i in range(len(X_data)):
#     m.flowrate.fix(X_data[i]/1e3)
#     m.head.fix(Y_data[i]/1e2)
#     calculate_variable_from_constraint(m.eff, m.surrogate_blk.pysmo_constraint["efficiency"])
#     Z_modeled.at[i, "efficiency_modeled"] = value(m.eff)*1e2

# # # Create 3D scatter plot
# fig2 = plt.figure(figsize=(10, 8))
# ax2 = fig2.add_subplot(111, projection="3d")

# Flatten the data for scatter plot
# X_flat = X.flatten()
# Y_flat = Y.flatten()
# Z_flat = Z.flatten()

# Remove NaN values for cleaner plot
# mask = ~np.isnan(Z_flat)
# X_clean = X_flat[mask]
# Y_clean = Y_flat[mask]
# Z_clean = Z_flat[mask]

# # Create scatter plot with color mapping
# scatter = ax2.scatter(X_data, Y_data, Z_data, color='red', s=50, alpha=0.6)

# # Plot actual data points
# ax2.scatter(X_data,Y_data,Z_modeled["efficiency_modeled"], color="black", s=20, label="Modeled Points")

# ax2.set_xlabel("Flow (gpm)")
# ax2.set_ylabel("Head (ft)")
# ax2.set_zlabel("Efficiency (%)")
# ax2.set_title("RO Feed Pump Efficiency 3D Scatter Plot")
# ax2.set_xlim(0, 4500)
# ax2.set_ylim(0, 400)
# plt.show()


# Get the absolute path of the current script
current_script_path = os.path.abspath(__file__)
# Get the directory containing the current script
current_directory = os.path.dirname(current_script_path)
surr_name = f"{pump_type}_total_efficiency_rbf_fit.json"
Surrogate.save_to_file(os.path.join(current_directory, surr_name), overwrite=True)
