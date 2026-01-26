import pandas as pd
import json

# Load the JSON file as a dictionary first
with open(r'C:\Users\rchurchi\flex_desal\src\wrd\surrogate_models\wpd_project.json', 'r') as f:
    data = json.load(f)

datasets = data['datasetColl']

# Convert specific datasets to DataFrames
pump_data = pd.DataFrame()
for dataset in datasets:
    name = dataset['name']
    df = pd.DataFrame(dataset['data'])
    try:
        eff = int(name[:2])  # Extract first 2 characters and convert to int
        df['efficiency'] = eff  # add efficiency column with every value the same
        pump_data = pd.concat([pump_data, df], ignore_index=True)
    except:
        pass
print(f"\nColumns before cleanup: {pump_data.columns.tolist()}")

# Drop the 'value' column if it exists
if 'value' in pump_data.columns:
    pump_data.drop(columns=['value'], inplace=True)

# Rename columns
pump_data.rename(columns={'x': 'Flow (gpm)', 'y': 'Head (ft)'}, inplace=True)

print(f"\nCombined pump data shape: {pump_data.shape}")
print(f"Columns after cleanup: {pump_data.columns.tolist()}")
print(pump_data.head())
print(f"\nUnique efficiencies: {sorted(pump_data['efficiency'].unique())}")

