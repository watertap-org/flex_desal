from idaes.apps.grid_integration import PriceTakerModel
import pandas as pd
import pyomo.environ as pyo
from pathlib import Path
from pricetaker.flowsheets import wrd_ro_flowsheet as fs
from pricetaker.flowsheets import utils
from pricetaker.flowsheets.params import FlexDesalParams
from watertap.core.solvers import get_solver
from idaes.core.util.model_diagnostics import DiagnosticsToolbox


if __name__ == "__main__":
    # Get the directory where this script is located
    script_dir = Path(__file__).parent
    price_data = pd.read_csv(script_dir / "sbce_pricesignal_short.csv")
    # price_data["Energy Rate"] = (
    #     price_data["electric_energy_on_peak"]
    #     + price_data["electric_energy_mid_peak"]
    #     + price_data["electric_energy_off_peak"]
    #     + price_data["electric_energy_super_off_peak"]
    # )
    # price_data["Fixed Demand Rate"] = price_data["electric_demand_fixed_summer"]
    # price_data["Var Demand Rate"] = price_data["electric_demand_peak_adder_summer"]
    price_data["Energy Rate"] = (
        price_data["electric_energy_0_2022-07-05_2022-07-14_0"]
        + price_data["electric_energy_1_2022-07-05_2022-07-14_0"]
        + price_data["electric_energy_2_2022-07-05_2022-07-14_0"]
        + price_data["electric_energy_3_2022-07-05_2022-07-14_0"]
    )
    price_data["Fixed Demand Rate"] = price_data[
        "electric_demand_maximum_2022-07-05_2022-07-14_0"
    ]
    price_data["Var Demand Rate"] = price_data[
        "electric_demand_peak-summer_2022-07-05_2022-07-14_0"
    ]
    # units? + Don't really understand the Var demand rate... it's only applied at peak hours, but it's billed monthly? But so is the fixed demand rate, so is it in addition to that?
    price_data["Emissions Intensity"] = 0
    # price_data["Customer Cost"] = price_data["electric_customer_fixed_charge"]
    price_data["Customer Cost"] = price_data[
        "electric_customer_0_2022-07-05_2022-07-14_0"
    ]
    m = PriceTakerModel()

    # Instantiate an object containing the model parameters
    m.params = FlexDesalParams(
        start_date="2022-07-05 00:00:00",
        end_date="2022-07-05 02:15:00",
        annual_production_AF=1000,
        # fixed_monthly_cost = 10000,
        # customer_rate=price_data["Customer Cost"][1],  # acrft/yr
    )
    m.params.intake.nominal_flowrate = 1063.5  # m3/hr
    m.params.wrd_uf.update(
        {
            "surrogate_type": "constant_energy_intensity",
            "surrogate_a": 1.0,
            "surrogate_b": 0.0,
            "nominal_recovery": 1,
        }
    )
    m.params.wrd_ro.update(
        {
            "startup_delay": 8,  # hours
            "minimum_downtime": 4,  # hours
            "nominal_flowrate": 363.4,  # m3/hr #per skid
            "surrogate_type": "constant_energy_intensity",
            "surrogate_a": 1.0,
            "surrogate_b": 0.0,
            "nominal_recovery": 0.465,
            "num_ro_skids": 4,
        }
    )

    # Append LMP data to the model
    m.append_lmp_data(lmp_data=price_data["Energy Rate"])

    m.build_multiperiod_model(
        flowsheet_func=fs.build_desal_flowsheet,
        flowsheet_options={"params": m.params},
    )

    # Update the time-varying parameters other than the LMP, such as
    # demand costs and emissions intensity. LMP value is updated by default

    # First, discover what blocks exist in the model

    print("Discovering blocks in period[1,1]:")
    for component in m.period[1, 1].component_objects(pyo.Block, descend_into=False):
        print(f"  - {component.name}")

    m.update_operation_params(
        {
            "fixed_demand_rate": price_data["Fixed Demand Rate"],
            "variable_demand_rate": price_data["Var Demand Rate"],
            "emissions_intensity": price_data["Emissions Intensity"],
            "customer_cost": price_data["Customer Cost"],
        }
    )

    # Add demand cost and fixed cost calculation constraints
    fs.add_demand_and_fixed_costs(m)

    # Add the startup delay constraints
    fs.add_delayed_startup_constraints(m)

    m.total_water_production = pyo.Expression(
        expr=m.params.timestep_hours
        * sum(m.period[:, :].posttreatment.product_flowrate)
    )
    m.total_energy_cost = pyo.Expression(expr=sum(m.period[:, :].energy_cost))
    m.total_demand_cost = pyo.Expression(
        expr=m.fixed_demand_cost + m.variable_demand_cost
    )
    m.total_customer_cost = pyo.Expression(
        expr=sum(m.period[:, :].customer_cost) * m.params.num_months
    )
    m.total_electricity_cost = pyo.Expression(
        expr=m.total_energy_cost + m.total_demand_cost + m.total_customer_cost
    )

    # Feed flow to the intake does not vary with time
    m.fix_operation_var("intake.feed_flowrate", m.params.intake.nominal_flowrate)

    # Pretreatment is either active (1) or inactive (0) for the entire run
    # m.fix_operation_var("pretreatment.op_mode", 1)

    fs.constrain_water_production(m)

    # If water recovery is static, it must be fixed
    if not m.params.wrd_ro.allow_variable_recovery:
        utils.wrd_fix_recovery(m, recovery=m.params.wrd_ro.nominal_recovery)

    m.obj = pyo.Objective(
        expr=m.total_energy_cost + m.total_demand_cost,
        sense=pyo.minimize,
    )

    # Can't use gurobi because it requires a liciense for integer variables
    # So going to use ipopt, but may need to look into this further
    dt = DiagnosticsToolbox(m)
    solver = get_solver()
    results = solver.solve(m)

    pyo.assert_optimal_termination(results)

    # Write optimal values of all operational variables to a csv file
    m.get_operation_var_values().to_csv("dummy_wrd_result.csv")

    # Plot operational variables
    fig, axs = m.plot_operation_profile(
        operation_vars=[
            "fixed_demand_rate",
            "variable_demand_rate",
            "posttreatment.product_flowrate",
            "num_skids_online",
        ],
    )
    fig.savefig("operation_profile.png")
    # Return the values of all variables and expressions that do not vary with time
    print(m.get_design_var_values())

    # OK this runs and solves at least
