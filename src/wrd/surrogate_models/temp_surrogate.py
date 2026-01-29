import pandas as pd
import numpy as np
import os
from idaes.core.surrogate.pysmo_surrogate import (
    PysmoRBFTrainer,
    PysmoPolyTrainer,
    PysmoSurrogate,
)

# load data
pump_data = pd.read_csv(
    os.path.join(os.path.dirname(__file__), "RO_UF_pump_head_curves_data.csv")
)
input_labels = ["Flow (gpm)"]
output_labels = ["Head 100% Speed (ft)"]
pump_data.dropna(inplace=True, subset=output_labels)
pump_data.drop(columns=["MCSF (ft)"], inplace=True)
input_data = pump_data[input_labels]
output_data = pump_data[output_labels]

# Scale Data
Data_scaled = pump_data
Data_scaled[output_labels[0]] = pump_data[output_labels[0]].mul(1e-2)  # head (ft)
Data_scaled[input_labels[0]] = pump_data[input_labels[0]].mul(1e-3)  # Flow (gpm)
min_flow = min(Data_scaled[input_labels[0]])
max_flow = max(Data_scaled[input_labels[0]])
min_head = min(Data_scaled[output_labels[0]])
max_head = max(Data_scaled[output_labels[0]])
input_bounds = {
    input_labels[0]: (min_flow, max_flow),
}

# Create Surrogate Type and trainer
trainer = PysmoPolyTrainer(
    input_labels=input_labels,
    output_labels=output_labels,
    training_dataframe=Data_scaled,  # NOT spliting data for training and validation because there's not that much data to begin with.
)
trainer.config.maximum_polynomial_order = 3
# Train Data
trained_surr = trainer.train_surrogate()

# Create Surrogate Model
Surrogate = PysmoSurrogate(
    trained_surr,
    input_labels,
    output_labels,
    input_bounds,
)
