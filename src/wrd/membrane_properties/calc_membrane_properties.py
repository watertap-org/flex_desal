import pandas as pd
from wrd.membrane_properties.ro_for_membrane_permeabilities import solve_ro_module
from pyomo.environ import (
    assert_optimal_termination,
    units as pyunits,
)
from watertap.core.solvers import get_solver
import os
import calendar
from idaes.core.util.model_statistics import degrees_of_freedom

def read_and_clean_input_data(train=1, stage_num=1, month=8, year=2021):
    """
    This function reads and cleans the input data for the specified RO stage and month/year.
    """

    folderpath = os.path.dirname(os.path.abspath(__file__))
    print(folderpath)

    # Load Data from 2019 - 2024
    # Load Data source with inlet P, outlet P, feed salinity, flowrate, perm production, and perm salinity.
    if stage_num == 3:
        # Load TSRO data
        data = pd.read_csv(
            os.path.join(folderpath, f"input_data/WRD_Data_TSRO_train{train}.csv")
        )
    else:
        # Load PRO data
        data = pd.read_csv(
            os.path.join(folderpath, f"input_data/WRD_Data_PRO_train{train}.csv")
        )

    # Convert all columns except DateTime to numeric
    for col in data.columns:
        if col != "DateTime":
            data[col] = pd.to_numeric(data[col], errors="coerce")

    # Exract only the month and year of interest
    data_mask_date = (data["DateTime"].str.startswith(f"{month}/")) & (
        data["DateTime"].str.contains(f"/{year}")
    )

    # Final cleaned data
    cleaned_data = data[data_mask_date].reset_index(drop=True)

    return cleaned_data


def calculate_additional_columns(cleaned_data, stage_num):
    """
    This function filters the cleaned data and calculates additional columns based on the stage number.
    """

    # Filter rows and calculate additional columns based on stage number

    if stage_num == 3:
        # Remove low flow days (off or off spec)
        data_mask_perm_flow = cleaned_data["stage 3 permeate flowrate (gpm)"] >= 120
        cleaned_data = cleaned_data[data_mask_perm_flow].reset_index(drop=True)

        # Remove low pressure days (off or off spec)
        data_mask_pressure = cleaned_data["stage 3 feed pressure (psi)"] >= 50
        cleaned_data = cleaned_data[data_mask_pressure].reset_index(drop=True)

        # Add column for salinity from conductivity using conversion factor 0.0005
        cleaned_data["stage 3 concentrate salinity (g/L)"] = (
            cleaned_data["stage 3 concentrate conductivity (us/cm)"] * 0.0005
        )
        cleaned_data["stage 3 permeate salinity (g/L)"] = (
            cleaned_data["stage 3 permeate conductivity (us/cm)"] * 0.0005
        )

        # Calculate feed flowrate
        cleaned_data["stage 3 feed flowrate (gpm)"] = (
            cleaned_data["stage 3 permeate flowrate (gpm)"]
            + cleaned_data["stage 3 concentrate flowrate (gpm)"]
        )

        # Calculate stage 3 feed salinity
        cleaned_data["stage 3 feed salinity (g/L)"] = (
            cleaned_data["stage 3 permeate flowrate (gpm)"]
            * cleaned_data["stage 3 permeate salinity (g/L)"]
            + cleaned_data["stage 3 concentrate flowrate (gpm)"]
            * cleaned_data["stage 3 concentrate salinity (g/L)"]
        ) / cleaned_data["stage 3 feed flowrate (gpm)"]

    else:
        # if stage_num == 1:
        # Calculate feed flowrate
        cleaned_data["stage 1 feed flowrate (gpm)"] = (
            cleaned_data["stage 1 permeate flowrate (gpm)"]
            + cleaned_data["stage 1 concentrate flowrate (gpm)"]
        )

        # Remove low flow days (off or off spec)
        data_mask_perm_flow = cleaned_data["stage 1 permeate flowrate (gpm)"] >= 1000
        cleaned_data = cleaned_data[data_mask_perm_flow].reset_index(drop=True)

        # Add column for salinity from conductivity using conversion factor 0.0005
        cleaned_data["stage 1 feed salinity (g/L)"] = (
            cleaned_data["stage 1 feed conductivity (us/cm)"] * 0.0005
        )
        cleaned_data["stage 1 permeate salinity (g/L)"] = (
            cleaned_data["stage 1 permeate conductivity (us/cm)"] * 0.0005
        )

        if stage_num == 2:
            # Calculate feed flowrate
            cleaned_data["stage 2 feed flowrate (gpm)"] = (
                cleaned_data["stage 2 concentrate flowrate (gpm)"]
                + cleaned_data["stage 2 permeate flowrate (gpm)"]
            )

            # Calculate stage 2 feed salinity
            cleaned_data["stage 2 feed salinity (g/L)"] = (
                cleaned_data["stage 1 feed flowrate (gpm)"]
                * cleaned_data["stage 1 feed salinity (g/L)"]
                + cleaned_data["stage 1 permeate flowrate (gpm)"]
                * cleaned_data["stage 1 permeate salinity (g/L)"]
            ) / cleaned_data["stage 1 concentrate flowrate (gpm)"]

            # Add column for salinity from conductivity using conversion factor 0.0005
            cleaned_data["stage 2 permeate salinity (g/L)"] = (
                cleaned_data["stage 2 permeate conductivity (us/cm)"] * 0.0005
            )

    return cleaned_data


def calc_membrane_permeability(cleaned_data, train=1, stage_num=1, month=8, year=2021):
    """
    This function calculates membrane permeability values (A and B) for the specified RO stage.
    """
    # Pass values to ro component model
    permability_values = {}

    for i in range(len(cleaned_data)):
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
        assert degrees_of_freedom(m) == 0
        solver = get_solver()
        results = solver.solve(m)
        assert_optimal_termination(results)

        A_new = m.fs.ro.unit.A_comp[0, "H2O"].value
        B_new = m.fs.ro.unit.B_comp[0, "NaCl"].value

        # Tabulate all the A and B Values
        permability_values[cleaned_data["DateTime"][i]] = {"A": A_new, "B": B_new}

    # Convert to DataFrame and export to CSV
    permeability_df = pd.DataFrame.from_dict(permability_values, orient="index")
    permeability_df.index.name = "DateTime"
    permeability_df.reset_index(inplace=True)

    return permeability_df


def main(train=1, stage_num=1, month=3, year=2021):
    # Read and clean input data
    cleaned_data = read_and_clean_input_data(
        train=train, stage_num=stage_num, month=month, year=year
    )
    cleaned_data = calculate_additional_columns(cleaned_data, stage_num)

    # Calculate membrane permeability values
    permeability_df = calc_membrane_permeability(
        cleaned_data=cleaned_data,
        train=train,
        stage_num=stage_num,
        month=month,
        year=year,
    )

    return cleaned_data, permeability_df


if __name__ == "__main__":

    # Make Changes here for different trains, stages, and months/years
    train = 1  # Change to 1, 2, 3, or 4 for different trains
    stage_num = 3  # Change to 1, 2, or 3 for different RO stages
    month = 3
    year = 2021

    cleaned_data, permeability_df = main(
        train=train, stage_num=stage_num, month=month, year=year
    )

    # Save cleaned data to CSV (optional)
    folder_path = os.path.dirname(__file__)

    if stage_num == 3:
        cleaned_data.to_csv(
            folder_path
            + f"/cleaned_data/cleaned_data_TSRO_train{train}_{calendar.month_abbr[month].lower()}{year}.csv",
            index=False,
        )
    else:
        cleaned_data.to_csv(
            folder_path
            + f"/cleaned_data/cleaned_data_PRO_train{train}_stage{stage_num}_{calendar.month_abbr[month].lower()}{year}.csv",
            index=False,
        )

    # Save to permeate_df to CSV
    folderpath = os.path.dirname(os.path.abspath(__file__))

    if stage_num == 3:
        permeability_df.to_csv(
            os.path.join(
                folderpath,
                f"output_data/memb_perm_values_TSRO_train{train}_{calendar.month_abbr[month].lower()}{year}.csv",
            ),
            index=False,
        )
    else:
        permeability_df.to_csv(
            os.path.join(
                folderpath,
                f"output_data/memb_perm_values_PRO_train{train}_stage{stage_num}_{calendar.month_abbr[month].lower()}{year}.csv",
            ),
            index=False,
        )

    # Print average A and B values
    average_A = permeability_df["A"].mean()
    average_B = permeability_df["B"].mean()

    print(f"Average A: {average_A}, Average B: {average_B}")
