import pytest
import pandas as pd
import numpy as np
from pyomo.environ import value, units as pyunits
from wrd.components.pump import main

# TODO: Consider applying the param sweep tool that already exists in WaterTAP if using this to generate data for a surrogate model.


def create_test_pairs(
    flow_ub=2700,
    flow_lb=2600,
    head_ub=270,
    head_lb=250,
    num_points=2,
    additional_points=None,
):

    test_pairs = pd.DataFrame(columns=["flow", "head"])
    flow_vals = np.linspace(flow_lb, flow_ub, num=num_points)
    head_vals = np.linspace(head_lb, head_ub, num=num_points)
    test_pairs["flow"] = np.repeat(flow_vals, num_points)
    test_pairs["head"] = np.tile(head_vals, num_points)
    if additional_points is not None:
        for point in additional_points:
            # Seems clunky
            new_row = pd.DataFrame({"flow": [point[0]], "head": [point[1]]})
            test_pairs = pd.concat([test_pairs, new_row], ignore_index=True)
    return test_pairs


def test_pump_param_sweep(test_pairs=None, pump_type="RO_feed", Pin=14.5):
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


def head_limit(flow, pump_type):
    """Helper function to define minimum head for given flowrate."""
    # Coefficients from pump curve data
    # Flow must be units of scaled gpm. Head is in ft.
    if pump_type == "RO_feed":
        b_0 = 3.127555
        b_1 = -0.09634
        b_2 = 0.051555
        b_3 = -0.026339
    elif pump_type == "RO_IS":
        b_0 = 1.02801
        b_1 = -0.21038
        b_2 = 0.284634
        b_3 = -0.253384
    elif pump_type == "UF":
        b_0 = 3.199441
        b_1 = -0.252764
        b_2 = 0.060008
        b_3 = -0.015999
    head = b_0 + b_1 * (flow) + b_2 * (flow) ** 2 + b_3 * (flow) ** 3
    return head


def max_flows(flow, pump_type):
    """Helper function to define maximum flow for given flowrate along the maximum pump capacity across different speeds."""
    # Flow must be in units of scaled gpm. Head is in scaled ft.
    if pump_type == "RO_feed":
        # Using end points for multi-speed head curves
        d_0 = 0.19338
        d_1 = -0.21424
        d_2 = 0.23183
        d_3 = -0.00942
    elif pump_type == "RO_IS":
        # No multispeed head curve available
        # I think using affinity laws to determine the max flow for this pump
        d_0 = 0
        d_1 = 0
        d_2 = 0.387369
        d_3 = 0
    elif pump_type == "UF":
        # No multispeed head curve available
        # I think using affinity laws to determine the max flow for this pump
        d_0 = 0
        d_1 = 0
        d_2 = 0.04354
        d_3 = 0
    head = d_0 + d_1 * (flow) + d_2 * (flow) ** 2 + d_3 * (flow) ** 3
    return head


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
    # Filter out ????

    # Filter out data points with speed less than 50%
    # pump_data = pump_data[pump_data['speed'] >= 50]

    return pump_data


if __name__ == "__main__":
    pump_type = "RO_feed"
    if pump_type == "RO_feed":
        flow_ub=3800,
        flow_lb=1000,
        head_ub=320,
        head_lb=100,
        additional_points = [(3000,254),(3000,240),(3330,230),(2640,230),(2640,250),(1980,270),(2280,270)] # RO Feed Pump
    elif pump_type == "RO_IS":
        flow_ub=1350,
        flow_lb=400,
        head_ub=100,
        head_lb=40
        additional_points = [(1300,66),(1200,56),(1250,61),(850,55),(850,65)] # RO IS Pump
    elif pump_type == "UF":
        flow_ub=5500,
        flow_lb=600,
        head_ub=280,
        head_lb=50,
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
        num_points=1,
        additional_points=additional_points,
    )
    print("Test Pairs:")
    print(test_pairs)
    filtered_test_pairs = filter_pump_test_points(test_pairs, pump_type=pump_type)
    print("Filtered test pairs:")
    print(filtered_test_pairs)
    dataset = test_pump_param_sweep(
        test_pairs=filtered_test_pairs, pump_type=pump_type, Pin=150
    )  # Not sure Pin really matters here
    print("Dataset:")
    print(dataset)
    dataset.rename(columns={"flow": "Flow (gpm)"}, inplace=True)
    dataset.rename(columns={"head": "Head (ft)"}, inplace=True)
    dataset.to_csv(f"{pump_type}_aff_laws_surr_2.csv", index=False)