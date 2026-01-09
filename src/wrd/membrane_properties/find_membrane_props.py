import pandas as pd
import datetime
from wrd.membrane_properties.ro_memb_test import *


# Load Data source with inlet P, outlet P, feed salinity, flowrate, perm production, and perm salinity.
Data = pd.read_csv("WRD_Data_PRO1.csv")

# Do any required data cleaning
data_mask = Data["DateTime"] == "2021-08-19 2:00:00"
cleaned_data = Data[data_mask].reset_index(drop=True)
# Add column for salinity from conductivity

# Pass values to ro component model

for i in range(len(cleaned_data)):
    m = main(
        Qin=cleaned_data["Feed Flowrate (gpm)"][i],
        Cin=cleaned_data["Feed Salinity (g/L)"][i],
        Tin=298,
        Pin=cleaned_data["Feed Pressure (psi)"][i],
        Pout=cleaned_data["Conc Pressure (psi)"][i],
    )
# Intialize, then unfix A and fix recovery or perm flowrate. Also unfix B and fix permeate salinity.

# Tabulate all the A and B Values

# Export to csv?