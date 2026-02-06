import pandas as pd
import json
from pathlib import Path

"""These functions are specificly designed to read output data 
   from WebPlotDigitizer for pump curves."""

# For RO pumps (no efficiency curve)
def load_ro_pump_data(file_name):
    # RO pump has multiple isoefficiency curves, so each must be extracted separately
    script_dir = Path(__file__).parent
    data_dir = script_dir / "pump_curve_data"
    file_path = data_dir / file_name
    with open(file_path, "r") as f:
        data = json.load(f)
    datasets = data["datasetColl"]

    # Convert specific datasets to DataFrames
    eff_data = pd.DataFrame()
    head_data = pd.DataFrame()
    for dataset in datasets:
        name = dataset["name"]
        df = pd.DataFrame(dataset["data"])
        try:
            eff = int(name[:2])  # Extract first 2 characters and convert to int
            df["efficiency"] = eff  # add efficiency column with every value the same
            df["Flow (gpm)"] = df["value"].apply(lambda x: x[0])
            df["Head (ft)"] = df["value"].apply(lambda x: x[1])
            df.drop(columns=["x", "y", "value"], inplace=True)
            eff_data = pd.concat([eff_data, df], ignore_index=True)
        except:
            df["Flow (gpm)"] = df["value"].apply(lambda x: x[0])
            df[name + " (ft)"] = df["value"].apply(lambda x: x[1])
            df.drop(columns=["x", "y", "value"], inplace=True)
            head_data = pd.concat([head_data, df], ignore_index=True)
        print(eff_data.head())
        eff_data.to_csv(data_dir / "RO_IS_pump_eff_curve_data.csv", index=False)
        print(head_data.head())
        head_data.to_csv(data_dir / "RO_IS_pump_head_curves_data.csv", index=False)

# For UF pump (with efficiency curve)
def load_uf_pump_data(file_name):
    # UF pump has one efficiency curve for 100% speed
    script_dir = Path(__file__).parent
    data_dir = script_dir / "pump_curve_data"
    file_path = data_dir / file_name
    with open(file_path, "r") as f:
        data = json.load(f)
    datasets = data["datasetColl"]
    # Convert specific datasets to DataFrames
    eff_data = pd.DataFrame()
    head_data = pd.DataFrame()
    for dataset in datasets:
        name = dataset["name"]
        print(name)
        df = pd.DataFrame(dataset["data"])
        if name == "Efficiency":
            df["Flow (gpm)"] = df["value"].apply(lambda x: x[0])
            df["Efficiency (%)"] = df["value"].apply(lambda x: x[1])
            df.drop(columns=["x", "y", "value"], inplace=True)
            eff_data = pd.concat([eff_data, df], ignore_index=True)
        else:
            df["Flow (gpm)"] = df["value"].apply(lambda x: x[0])
            df[name + " (ft)"] = df["value"].apply(lambda x: x[1])
            df.drop(columns=["x", "y", "value"], inplace=True)
            head_data = pd.concat([head_data, df], ignore_index=True)
    print(eff_data.head())
    eff_data.to_csv(
        data_dir / "RO_UF_pump_eff_curve_data.csv",
        index=False,
    )
    print(head_data.head())
    head_data.to_csv(
        data_dir / "RO_UF_pump_head_curves_data.csv",
        index=False,
    )


if __name__ == "__main__":
    filename = "UF_pump_curves.json"
    load_uf_pump_data(filename)
    # filename = "RO_IS_pump_curves.json"
    # load_ro_pump_data(filename)
