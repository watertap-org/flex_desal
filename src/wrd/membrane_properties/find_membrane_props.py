import pandas as pd
import datetime
from wrd.membrane_properties.ro_memb_test import solve_ro_module
from pyomo.environ import (
    assert_optimal_termination,
    units as pyunits,
)
from watertap.core.solvers import get_solver

# Load Data source with inlet P, outlet P, feed salinity, flowrate, perm production, and perm salinity.
Data = pd.read_csv("C:\\Users\\rchurchi\\flex_desal\\src\\wrd\\membrane_properties\\WRD_Data_PRO1.csv")

# Convert all columns except DateTime to numeric
for col in Data.columns:
    if col != "DateTime":
        Data[col] = pd.to_numeric(Data[col], errors='coerce')

# Do any required data cleaning
data_mask_date = (Data["DateTime"].str.startswith("3/")) & (Data["DateTime"].str.contains("/2021"))
cleaned_data = Data[data_mask_date].reset_index(drop=True)
data_mask_perm_flow = cleaned_data["stage 1 permeate flowrate (gpm)"] >= 1000
cleaned_data = cleaned_data[data_mask_perm_flow].reset_index(drop=True)

# Add column for salinity from conductivity
cleaned_data["stage 1 feed salinity (g/L)"] = cleaned_data["stage 1 feed conductivity (us/cm)"] * 0.0005 # Conversion factor
cleaned_data["stage 1 permeate salinity (g/L)"] = cleaned_data["stage 1 permeate conductivity (us/cm)"] * 0.0005
cleaned_data["stage 1 feed flowrate (gpm)"] = cleaned_data["stage 1 permeate flowrate (gpm)"] + cleaned_data["stage 1 concentrate flowrate (gpm)"]
print(cleaned_data)

if __name__ == '__main__':
    # Pass values to ro component model
    permability_values = {}
    for i in range(len(cleaned_data)):
        print()
        Qin=cleaned_data["stage 1 feed flowrate (gpm)"][i] #* pyunits.gallons/pyunits.minute
        Qperm = cleaned_data["stage 1 permeate flowrate (gpm)"][i] #* pyunits.gallons/pyunits.minute
        Cin = cleaned_data["stage 1 feed salinity (g/L)"][i] #* pyunits.g/pyunits.L
        Cperm = cleaned_data["stage 1 permeate salinity (g/L)"][i] #* pyunits.g/pyunits.L 
        Pin = cleaned_data["stage 1 feed pressure (psi)"][i]
        Pout = cleaned_data["stage 1 concentrate pressure (psi)"][i]

        m = solve_ro_module(
            Qin=Qin,
            Cin=Cin,
            Tin=298,
            Pin=Pin,
            Pout=Pout,
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
        A_new = m.fs.ro.unit.A_comp[0,"H2O"].value
        B_new = m.fs.ro.unit.B_comp[0,"NaCl"].value
        # print(f"DateTime: {cleaned_data['DateTime'][i]}")
        # print(f"Run {i+1}: A = {A_new}, B = {B_new}")
        # Tabulate all the A and B Values
        permability_values[cleaned_data['DateTime'][i]] = {"A": A_new, "B": B_new}
    # Export to csv?
    print(permability_values)
    
    # Convert to DataFrame and export to CSV
    permeability_df = pd.DataFrame.from_dict(permability_values, orient='index')
    permeability_df.index.name = 'DateTime'
    permeability_df.reset_index(inplace=True)
    permeability_df.to_csv("C:\\Users\\rchurchi\\flex_desal\\src\\wrd\\membrane_properties\\memb_perm_values_march.csv", index=False)
    average_A = sum([v["A"] for v in permability_values.values()]) / len(permability_values)
    average_B = sum([v["B"] for v in permability_values.values()]) / len(permability_values)
    print(f"Average A: {average_A}, Average B: {average_B}")
