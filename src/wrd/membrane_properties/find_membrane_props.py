import pandas as pd
from wrd.membrane_properties.ro_memb_test import solve_ro_module
from pyomo.environ import (
    assert_optimal_termination,
    units as pyunits,
)
from watertap.core.solvers import get_solver

# Choose stage number to find properties
stage_num = 2

# Load Data from 2019 - 2024
if stage_num == 3:
    Data = pd.read_csv(
        "C:\\Users\\rchurchi\\flex_desal\\src\\wrd\\membrane_properties\\WRD_Data_TSRO1.csv"
    )
else:
    # Load Data source with inlet P, outlet P, feed salinity, flowrate, perm production, and perm salinity.
    Data = pd.read_csv(
        "C:\\Users\\rchurchi\\flex_desal\\src\\wrd\\membrane_properties\\WRD_Data_PRO1.csv"
    )

# Convert all columns except DateTime to numeric
for col in Data.columns:
    if col != "DateTime":
        Data[col] = pd.to_numeric(Data[col], errors="coerce")

# Exract only the month and year of interest
data_mask_date = (Data["DateTime"].str.startswith("8/")) & (
    Data["DateTime"].str.contains("/2021")
)
cleaned_data = Data[data_mask_date].reset_index(drop=True)

# Add variables
if stage_num == 3:
    # Remove low flow days (off or off spec)
    data_mask_perm_flow = cleaned_data["stage 3 permeate flowrate (gpm)"] >= 120
    cleaned_data = cleaned_data[data_mask_perm_flow].reset_index(drop=True)
    data_mask_pressure = cleaned_data["stage 3 feed pressure (psi)"] >= 50
    cleaned_data = cleaned_data[data_mask_pressure].reset_index(drop=True)
    cleaned_data["stage 3 concentrate salinity (g/L)"] = (
        cleaned_data["stage 3 concentrate conductivity (us/cm)"] * 0.0005
    )  # Conversion factor
    cleaned_data["stage 3 permeate salinity (g/L)"] = (
        cleaned_data["stage 3 permeate conductivity (us/cm)"] * 0.0005
    )
    cleaned_data["stage 3 feed flowrate (gpm)"] = (
        cleaned_data["stage 3 permeate flowrate (gpm)"]
        + cleaned_data["stage 3 concentrate flowrate (gpm)"]
    )
    cleaned_data["stage 3 feed salinity (g/L)"] = (
        cleaned_data["stage 3 permeate flowrate (gpm)"]
        * cleaned_data["stage 3 permeate salinity (g/L)"]
        + cleaned_data["stage 3 concentrate flowrate (gpm)"]
        * cleaned_data["stage 3 concentrate salinity (g/L)"]
    ) / cleaned_data["stage 3 feed flowrate (gpm)"]

else:
    # Remove low flow days (off or off spec)
    data_mask_perm_flow = cleaned_data["stage 1 permeate flowrate (gpm)"] >= 1000
    cleaned_data = cleaned_data[data_mask_perm_flow].reset_index(drop=True)
    # Add column for salinity from conductivity
    cleaned_data["stage 1 feed salinity (g/L)"] = (
        cleaned_data["stage 1 feed conductivity (us/cm)"] * 0.0005
    )  # Conversion factor
    cleaned_data["stage 1 permeate salinity (g/L)"] = (
        cleaned_data["stage 1 permeate conductivity (us/cm)"] * 0.0005
    )
    cleaned_data["stage 1 feed flowrate (gpm)"] = (
        cleaned_data["stage 1 permeate flowrate (gpm)"]
        + cleaned_data["stage 1 concentrate flowrate (gpm)"]
    )
    if stage_num == 2:
        cleaned_data["stage 2 feed salinity (g/L)"] = (
            cleaned_data["stage 1 feed flowrate (gpm)"]
            * cleaned_data["stage 1 feed salinity (g/L)"]
            - -cleaned_data["stage 1 permeate flowrate (gpm)"]
            * cleaned_data["stage 1 permeate salinity (g/L)"]
        ) / cleaned_data["stage 1 concentrate flowrate (gpm)"]
        cleaned_data["stage 2 feed flowrate (gpm)"] = (
            cleaned_data["stage 2 concentrate flowrate (gpm)"]
            + cleaned_data["stage 2 permeate flowrate (gpm)"]
        )
        cleaned_data["stage 2 permeate salinity (g/L)"] = (
            cleaned_data["stage 2 permeate conductivity (us/cm)"] * 0.0005
        )
print(cleaned_data)

if __name__ == "__main__":
    # Pass values to ro component model
    permability_values = {}
    for i in range(len(cleaned_data)):
        print()
        Qin = cleaned_data[f"stage {stage_num} feed flowrate (gpm)"][i]
        Qperm = cleaned_data[f"stage {stage_num} permeate flowrate (gpm)"][i]
        Cin = cleaned_data[f"stage {stage_num} feed salinity (g/L)"][i]
        Cperm = cleaned_data[f"stage {stage_num} permeate salinity (g/L)"][i]
        Pin = cleaned_data[f"stage {stage_num} feed pressure (psi)"][i]
        Pout = cleaned_data[f"stage {stage_num} concentrate pressure (psi)"][i]

        m = solve_ro_module(
            Qin=Qin,
            Cin=Cin,
            Tin=298,
            Pin=Pin,
            Pout=Pout,
            stage_num=stage_num,
        )

        # Intialize, then unfix A and fix recovery or perm flowrate. Also unfix B and fix permeate salinity.
        stage_RR = Qperm / Qin
        m.fs.ro.unit.A_comp.unfix()
        m.fs.ro.unit.B_comp.unfix()
        m.fs.ro.unit.recovery_vol_phase[0, "Liq"].fix(stage_RR)
        m.fs.ro.unit.mixed_permeate[0].conc_mass_phase_comp["Liq", "NaCl"].fix(Cperm)
        solver = get_solver()
        results = solver.solve(m)
        assert_optimal_termination(results)
        A_new = m.fs.ro.unit.A_comp[0, "H2O"].value
        B_new = m.fs.ro.unit.B_comp[0, "NaCl"].value

        # Tabulate all the A and B Values
        permability_values[cleaned_data["DateTime"][i]] = {"A": A_new, "B": B_new}
    
    print(permability_values)
    average_A = sum([v["A"] for v in permability_values.values()]) / len(
        permability_values
    )
    average_B = sum([v["B"] for v in permability_values.values()]) / len(
        permability_values
    )
    print(f"Average A: {average_A}, Average B: {average_B}")

    # Convert to DataFrame and export to CSV
    permeability_df = pd.DataFrame.from_dict(permability_values, orient="index")
    permeability_df.index.name = "DateTime"
    permeability_df.reset_index(inplace=True)
    if stage_num == 3:
        permeability_df.to_csv(
            f"C:\\Users\\rchurchi\\flex_desal\\src\\wrd\\membrane_properties\\memb_perm_values_TSRO1_aug.csv",
            index=False,
        )
    else:
        permeability_df.to_csv(
            f"C:\\Users\\rchurchi\\flex_desal\\src\\wrd\\membrane_properties\\memb_perm_values_PRO1_{stage_num}_aug.csv",
            index=False,
        )

