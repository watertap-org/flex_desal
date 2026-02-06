import pytest
import pandas as pd
import numpy as np
from pathlib import Path
from pyomo.environ import value, units as pyunits
from wrd.components.pump import main
from wrd.surrogate_models.create_pump_surrogate import head_limit, max_flows

"""pump_param_sweep function will perform a parameter sweep over head and flows using the pump affinity laws. 

# Inputs are:
    - number of points for each variable in the sweep (total points will be num_points^2, minus those filtered out)
    - additional_points: list of tuples of (flow, head) points to include in the sweep
    - Head and flow upper and lower bounds
    - Pump type: determines the surrogate curve for max head at 100% speed, and for efficiency at 100% speed 
        ** These curves are hard coded into the pump model at the moment. They were generated using temp_surrogate. 
        to make this more general, the data could be input and temp_surrogate called to gen curves.
        ** Specific scaling used for the flow and head values in the curve surrogates (1e-2 for head and 1e-3 for flow).
"""

__all__ = [
    "create_test_pairs",
    "filter_pump_test_points",
    "pump_param_sweep",
]


def create_test_pairs(
    flow_ub=2700,
    flow_lb=2600,
    head_ub=270,
    head_lb=250,
    num_points=2,
    additional_points=None,
):

    flow_vals = np.linspace(flow_lb, flow_ub, num=num_points)
    head_vals = np.linspace(head_lb, head_ub, num=num_points)
    flow_grid, head_grid = np.meshgrid(flow_vals, head_vals, indexing="ij")
    test_pairs = pd.DataFrame(
        {"flow": flow_grid.flatten(), "head": head_grid.flatten()}
    )
    if additional_points is not None:
        for point in additional_points:
            # Seems clunky
            new_row = pd.DataFrame({"flow": [point[0]], "head": [point[1]]})
            test_pairs = pd.concat([test_pairs, new_row], ignore_index=True)
    return test_pairs


def filter_pump_test_points(pump_data, pump_type="RO_feed"):
    # Filter out points above upper bound of pump curve
    pump_data = pump_data[
        pump_data["head"]
        <= head_limit(pump_data["flow"] / 1e3, pump_type=pump_type) * 1e2
    ]
    # Filter out points to the right of max flow at diff. speeds
    pump_data = pump_data[
        pump_data["head"]
        >= max_flows(pump_data["flow"] / 1e3, pump_type=pump_type) * 1e2
    ]
    return pump_data


def pump_param_sweep(test_pairs=None, pump_type="RO_feed", Pin=14.5):
    # RO Feed Pump
    # Flow and Head pairs
    if test_pairs is None:
        test_pairs = pd.DataFrame(
            [[2778.0, 271.1], [2778.0, 150.8]], columns=["flow", "head"]
        )
    outputs = pd.DataFrame(
        columns=["speed", "pump_efficiency", "total_efficiency", "power"]
    )
    if pump_type == "RO_feed":
        stage_num = 1
        uf = False
    elif pump_type == "RO_IS":
        stage_num = 2
        uf = False
    elif pump_type == "UF":
        stage_num = 1
        uf = True
    for i in range(len(test_pairs)):
        flow = test_pairs.iloc[i, 0]
        head = test_pairs.iloc[i, 1]
        m = main(
            Qin=flow,
            head=head,
            Cin=0.5,
            Tin=298.15,
            Pin=Pin,  # psi
            stage_num=stage_num,
            uf=uf,
            file="wrd_inputs_8_19_21.yaml",
            add_costing=True,
        )
        speed = value(m.fs.pump.unit.eff.speed)
        print(f"Pump speed for flow {flow} and head {head}: {speed}")
        rpm = value(m.fs.pump.unit.eff.speed * 1780)
        print(f"Pump RPM: {rpm}")
        fluid_efficiency = value(m.fs.pump.unit.eff.efficiency_fluid)
        total_efficiency = value(m.fs.pump.unit.efficiency_pump[0])
        power = value(
            pyunits.convert(m.fs.pump.unit.work_mechanical[0], to_units=pyunits.kW)
        )
        print(f"Fluid efficiency: {fluid_efficiency}")
        outputs.at[i, "pump_efficiency"] = fluid_efficiency * 100
        outputs.at[i, "total_efficiency"] = total_efficiency * 100
        outputs.at[i, "speed"] = speed * 100
        outputs.at[i, "power"] = power
    test_pairs.reset_index(drop=True, inplace=True)
    dataset = pd.concat([test_pairs, outputs], axis=1)
    return dataset


if __name__ == "__main__":
    pump_type = "RO_IS"
    num_points = 8

    if pump_type == "RO_feed":
        flow_ub = 3800
        flow_lb = 1000
        head_ub = 320
        head_lb = 100
        additional_points = [
            (1760, 290),
            (1840, 289),
            (1930, 288),
            (2110, 280),
            (2725, 263),
            (2000, 254),
            (3350, 225),
        ]  # RO Feed Pump
    elif pump_type == "RO_IS":
        flow_ub = 1400
        flow_lb = 400
        head_ub = 110
        head_lb = 40
        additional_points = [
            (778, 90),
            (840, 88),
            (590, 92),
            (715, 45),
            (830, 45),
            (965, 45),
        ]  # RO IS Pump
    elif pump_type == "UF":
        flow_ub = 5500
        flow_lb = 600
        head_ub = 280
        head_lb = 50
        additional_points = [
            (4690, 145),
            (4690, 110),
            (3585, 215),
            (3065, 235),
            (2000, 200),
            (2000, 250),
            (2000, 100),
            (2000, 75),
            (2600, 150),
            (2600, 75),
            (2600, 115),
        ]  # UF Pump

    test_pairs = create_test_pairs(
        flow_ub=flow_ub,
        flow_lb=flow_lb,
        head_ub=head_ub,
        head_lb=head_lb,
        num_points=num_points,
        additional_points=additional_points,
    )
    print("Test Pairs:")
    print(test_pairs)
    filtered_test_pairs = filter_pump_test_points(test_pairs, pump_type=pump_type)
    print("Filtered test pairs:")
    print(filtered_test_pairs)
    dataset = pump_param_sweep(
        test_pairs=filtered_test_pairs, pump_type=pump_type, Pin=150
    )  # Not sure Pin really matters here
    print("Dataset:")
    print(dataset)
    dataset.rename(columns={"flow": "Flow (gpm)"}, inplace=True)
    dataset.rename(columns={"head": "Head (ft)"}, inplace=True)

    data_dir = Path(__file__).parent / "surrogate_data"
    data_dir.mkdir(exist_ok=True)
    dataset.to_csv(
        data_dir / f"{pump_type}_pump_aff_law_data_for_surrogate.csv", index=False
    )

    # test_points = pd.DataFrame([[4690, 145], [2000, 75],[2600, 150]], columns=["flow", "head"])
    # dataset = pump_param_sweep(test_pairs=test_points, pump_type="UF")
