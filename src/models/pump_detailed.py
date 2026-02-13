from networkx import config
import pandas as pd
from numpy import polyfit
from pyomo.environ import (
    Suffix,
    check_optimal_termination,
)
from watertap.core.solvers import get_solver

from pyomo.common.config import ConfigValue, In
from pyomo.environ import Var, Param, Set, units as pyunits, Expr_if, value

from idaes.core.util.config import DefaultBool
from enum import Enum, auto

# Import IDAES cores
from idaes.models.unit_models.pressure_changer import PumpData
from idaes.core import declare_process_block_class
import idaes.core.util.scaling as iscale
from idaes.core.util.exceptions import InitializationError

import idaes.logger as idaeslog

from watertap.core import InitializationMixin
from watertap.costing.unit_models.pump import cost_pump
from watertap.costing.unit_models.energy_recovery_device import (
    cost_energy_recovery_device,
)

_log = idaeslog.getLogger(__name__)


class VariableEfficiency(Enum):
    Fixed = auto()  # default is constant efficiency
    Flow = auto()  # flow-only correlation


class PumpCurveDataType(Enum):
    DataSet = (
        auto(),
    ) 
    # user provides flow, head, and efficiency data to fit a curve and surrogate coefficients are calculated via polyfit
    SurrogateCoefficent = (
        auto()
    ) 
    # curve fit correlation based on flow and head is provided by the user as surrogate coefficients.


@declare_process_block_class("Pump")
class PumpIsothermalData(InitializationMixin, PumpData):
    """
    Detailed Isothermal Pump Unit Model Class

    * User can specify whether to use a variable efficiency based on flow and head, or to use a constant efficiency at the best efficiency point (BEP).
    * This function builds on the standard pump model by adding variables and constraints to represent pump performance curves.
    * The model includes options for defining the pump curve based on flow and head at the BEP, as well as scaling the efficiency curve accordingly.
    * This allows for more accurate representation of pump performance under varying operating conditions.
    """

    def _validate_curve_data(val):
        if val is None:
            return val
        if isinstance(val, str):
            try:
                df = pd.read_csv(val)
            except Exception as e:
                raise ValueError(f"Failed to read CSV file '{val}': {e}")

            required_cols = ["flow (m3/s)", "head (m)", "efficiency (-)"]
            if all(col in df.columns for col in required_cols):
                return df
            raise ValueError(
                f"CSV file must contain columns: {required_cols}, but found: {list(df.columns)}"
            )
        raise ValueError("Must be a string filepath to a CSV file or None")

    CONFIG = PumpData.CONFIG()

    CONFIG.declare(
        "variable_efficiency",
        ConfigValue(
            default=VariableEfficiency.Fixed,
            domain=In(VariableEfficiency),
            description="Variable pump efficiency flag",
            doc="""Indicates whether pump efficiency should be variable or fixed.
                ** Fixed: The default is constant efficiency at the best efficiency point (BEP)
                ** Flow: Efficiency is a function of flow only at reference speed based on pump performance curve data provided.
            """,
        ),
    )

    # Config defining if surrogate coefficients should be calculated based or hardcoded.
    CONFIG.declare(
        "pump_curve_data_type",
        ConfigValue(
            default=PumpCurveDataType.SurrogateCoefficent,
            domain=In(PumpCurveDataType),
            description="Type of pump curve data",
            doc="""Indicates whether the pump curve data is provided as a dataset or as surrogate coefficients.""",
        ),
    )

    CONFIG.declare(
        "head_surrogate_coeffs",
        ConfigValue(
            default={},
            domain=dict,
            description="Surrogate coefficients for the pump head curve based on flow only",
            doc="""Coefficients for the pump head curve surrogate based on flow only. 
            The head is calculated as a function of flow using a cubic polynomial with these coefficients.""",
        ),
    )

    CONFIG.declare(
        "eff_surrogate_coeffs",
        ConfigValue(
            default={},
            domain=dict,
            description="Surrogate coefficients for the pump efficiency curve based on flow only",
            doc="""Coefficients for the pump efficiency curve surrogate based on flow only. 
            The efficiency is calculated as a function of flow using a cubic polynomial with these coefficients.""",
        ),
    )

    CONFIG.declare(
        "pump_curves",
        ConfigValue(
            default=None,
            domain=_validate_curve_data,
            description="head curve at reference speed",
            doc="Data from digitization of head vs. flow and efficiency vs. flow curves at a reference speed from a pump datasheet",
        ),
    )

    def build(self):
        super().build()

        # Isothermal pump set up
        if hasattr(self.control_volume, "enthalpy_balances"):
            self.control_volume.del_component(self.control_volume.enthalpy_balances)

        @self.control_volume.Constraint(
            self.flowsheet().config.time, doc="Isothermal constraint"
        )
        def isothermal_balance(b, t):
            return b.properties_in[t].temperature == b.properties_out[t].temperature

        # Add efficiency variables for more and VFD efficiency
        self.vfd_efficiency = Param(
            initialize=0.97,
            doc="Variable frequency drive (VFD) efficiency",
            units=pyunits.dimensionless,
        )

        self.motor_efficiency = Param(
            initialize=0.95,
            doc="Motor efficiency",
            units=pyunits.dimensionless,
        )

        if self.config.variable_efficiency is not VariableEfficiency.Fixed:
            # Variable efficiency pump set-up
            #### Design point variables ####
            self.design_flow = Var(
                initialize=1.0,
                doc="""Design flowrate of the centrifugal pump. 
                This could be the flowrate at the best efficiency point (BEP) or another reference point.
                Used to build the system curve.""",
                units=pyunits.m**3 / pyunits.s,
            )

            self.design_head = Var(
                initialize=1.0,
                doc="""Design head of the centrifugal pump. 
                This could be the head at the best efficiency point (BEP) or another reference point.
                Used to build the system curve.""",
                units=pyunits.m,
            )

            self.design_efficiency = Var(
                initialize=0.8,
                doc="""Design efficiency of the centrifugal pump. 
                This could be the efficiency at the best efficiency point (BEP) or another reference point.""",
                units=pyunits.dimensionless,
            )

            self.design_speed_fraction = Var(
                initialize=0.8,
                doc="""Design speed fraction of the centrifugal pump.
                This could be the speed fraction at the best efficiency point (BEP) or another reference point.""",
                units=pyunits.dimensionless,
            )

            ### System curve variables ###
            # Head = system_curve_geometric_head +  system_curve_flow_constant * (flow)**2

            self.system_curve_geometric_head = Var(
                initialize=0.0,
                doc="""Geometric head constant for the pump, that represents the static head component used to define the system curve.""",
                units=pyunits.m,
            )

            self.system_curve_flow_constant = Var(
                initialize=1.0,
                doc="""Geometric flow constant for the pump, represents the major and minor losses in pump. Used to define the system curve.""",
                units=pyunits.m * (pyunits.m**3 / pyunits.s) ** (-2),
            )
            
            self.efficiency_adjustment_factor = Param(
                initialize = 0,
                doc = "Adjustment factor for design_efficiency_constraint, reducing pump efficiency at lower speeds.",
                units = pyunits.dimensionless
            )

            # Constraints connecting inlet and outlet conditions to the design point head and flow, used to solve for the system curve constants
            @self.Constraint(
                doc="Design head is the pressure difference across the pump at the design point"
            )
            def design_head_constraint(b):
                return b.design_head == b.control_volume.deltaP[0] / (
                    b.control_volume.properties_out[0].dens_mass_phase["Liq"]
                    * 9.81
                    * pyunits.m
                    / pyunits.s**2
                )

            @self.Constraint(
                doc="Design flow is the flow through the pump at the design point"
            )
            def design_flow_constraint(b):
                return b.design_flow == pyunits.convert(
                    b.control_volume.properties_in[0].flow_vol_phase["Liq"],
                    to_units=pyunits.m**3 / pyunits.s,
                )

            # Constraints to calculate system curve constants based on design point
            @self.Constraint(doc="System curve head constant calculation")
            def system_curve_head_constant_calculation(b):
                return (
                    b.design_head
                    == b.system_curve_geometric_head
                    + b.system_curve_flow_constant * (b.design_flow**2)
                )

            #### Pump curve variables ####
            self.ref_flow = Var(
                initialize=1.0,
                doc="Reference flowrate for the pump on the pump curve from specification sheet",
                units=pyunits.m**3 / pyunits.s,
            )

            self.ref_head = Var(
                initialize=10.0,
                doc="Reference head for the pump on the pump curve from specification sheet.",
                units=pyunits.m,
            )

            self.ref_efficiency = Var(
                initialize=0.8,
                doc="Reference efficiency for the pump on the pump curve from specification sheet.",
                units=pyunits.dimensionless,
            )

            self.ref_speed_fraction = Var(
                initialize=1,
                doc="Reference speed fraction for the pump on the pump curve from specification sheet.",
                units=pyunits.dimensionless,
            )

            if self.config.pump_curve_data_type == PumpCurveDataType.DataSet:
                # Read the dataset file path and create a constraint to fit the surrogate coefficients based on the dataset provided by the user
                # Check there is a dataset
                # if not self.config.pump_curves:

                self.surrogate_index = Set(
                    initialize=[0, 1, 2, 3], doc="Index for surrogate coefficients"
                )
                # pump_curves is converted from a filepath name to DataFrame from the validator
                curves_df = self.config.pump_curves

                p_head = polyfit(curves_df["flow (m3/s)"], curves_df["head (m)"], 3)
                p_eff = polyfit(curves_df["flow (m3/s)"], curves_df["efficiency (-)"], 3)

                head_surrogate_coeffs = {
                    i: float(p_head[3 - i]) for i in self.surrogate_index
                }               
                eff_surrogate_coeffs = {
                    i: float(p_eff[3 - i]) for i in self.surrogate_index
                }

            elif (
                self.config.pump_curve_data_type
                == PumpCurveDataType.SurrogateCoefficent
            ):
                # Surrogate coefficients based on flow only for head and efficiency calculation
                if (
                    not self.config.head_surrogate_coeffs
                    or not self.config.eff_surrogate_coeffs
                ):
                    raise ValueError(
                        "surrogate_coeffs must be provided for the pump head curve and efficiency curve when pump_curve_data_type is set to SurrogateCoefficent."
                    )

                self.surrogate_index = Set(
                    initialize=[0, 1, 2, 3], doc="Index for surrogate coefficients"
                )

                if list(self.config.head_surrogate_coeffs.keys()) != [0, 1, 2, 3]:
                    raise ValueError(
                        "head_surrogate_coeffs keys must match the surrogate_index set [0,1,2,3] where each key corresponds to the coefficient corespond to the order of each term. Ex: {0: 10, 1: 2, 2: 0.5, 3: 0.1} corresponds to the surrogate curve: = 10 + 2*flow + 0.5*flow^2 + 0.1*flow^3"
                    )
                if list(self.config.eff_surrogate_coeffs.keys()) != [0, 1, 2, 3]:
                    raise ValueError(
                        "eff_surrogate_coeffs keys must match the surrogate_index set [0,1,2,3] where each key corresponds to the coefficient corespond to the order of each term. Ex: {0: 10, 1: 2, 2: 0.5, 3: 0.1} corresponds to the surrogate curve: = 10 + 2*flow + 0.5*flow^2 + 0.1*flow^3"
                    )
                
                head_surrogate_coeffs = self.config.head_surrogate_coeffs
                eff_surrogate_coeffs = self.config.eff_surrogate_coeffs

            else:
                raise ValueError(
                    "Invalid pump curve data type specified. Must be either DataSet or SurrogateCoefficent."
                )

            # Also would be nice for fit to be in terms of ft and gpm, or m and m3/hr or m3/s
            # Note: Coefficients are dimensionless because they're applied to flow in m³/s (numeric value)
            # The polynomial produces head in meters
            self.head_surrogate_coefficients = Param(
                self.surrogate_index,
                initialize=head_surrogate_coeffs,
                doc="Coefficients for the head surrogate based on flow only (flow in m³/s, output in m)",
                units=pyunits.dimensionless,
            )

            self.efficiency_surrogate_coefficients = Param(
                self.surrogate_index,
                initialize=eff_surrogate_coeffs,
                doc="Coefficients for the efficiency surrogate based on flow only",
                units=pyunits.dimensionless,
            )

            # Constraints to calculate reference point variables by solving the system curve and pump curve simultaneously at the reference point.
            @self.Constraint(doc="Reference head equality")
            def ref_head_constraint(b):
                return (
                    b.ref_head
                    == b.system_curve_geometric_head
                    + b.system_curve_flow_constant * (b.ref_flow**2)
                )

            @self.Constraint(
                doc="Reference head calculation calculated using the pump curve surrogate"
            )
            def ref_head_surrogate(b):
                # Convert flow to dimensionless value in m³/s for polynomial
                # flow_val = pyunits.convert(b.ref_flow, to_units=pyunits.m**3/pyunits.s)
                return b.ref_head == (
                    b.head_surrogate_coefficients[0]
                    + b.head_surrogate_coefficients[1] * (b.ref_flow / (pyunits.m**3 /pyunits.s))
                    + b.head_surrogate_coefficients[2] * ((b.ref_flow / (pyunits.m**3 /pyunits.s))**2)
                    + b.head_surrogate_coefficients[3] * ((b.ref_flow / (pyunits.m**3 /pyunits.s))**3)
                ) * pyunits.m

            # Calculate the reference efficiency based on the surrogate coefficients and reference flow
            @self.Constraint(
                doc="Reference efficiency calculated calculated using the pump curve surrogate"
            )
            def ref_efficiency_constraint(b):
                # Convert flow to dimensionless value in m³/s for polynomial
                return b.ref_efficiency == (
                    b.efficiency_surrogate_coefficients[0]
                    + b.efficiency_surrogate_coefficients[1] * (b.ref_flow / (pyunits.m**3 / pyunits.s))
                    + b.efficiency_surrogate_coefficients[2] * ((b.ref_flow / ( pyunits.m**3 / pyunits.s))**2)
                    + b.efficiency_surrogate_coefficients[3] * ((b.ref_flow / ( pyunits.m**3 / pyunits.s))**3)                )

            @self.Constraint(
                doc="Design speed fraction calculation using affinity laws"
            )
            def design_speed_fraction_constraint(b):
                return b.design_flow == b.ref_flow * (
                    b.design_speed_fraction / b.ref_speed_fraction
                )


            # Expression to calculate the efficiency at the design point based as a function of speed fraction
            @self.Constraint(doc="Design efficiency calculation")
            def design_efficiency_constraint(b):
                return (
                    b.design_efficiency
                    == b.ref_efficiency
                    * (b.design_speed_fraction / b.ref_speed_fraction) ** b.efficiency_adjustment_factor
                )

            @self.Constraint(
                doc="Overall efficiency calculation including motor efficiency and VFD efficiency"
            )
            def overall_efficiency_constraint(b):
                return (
                    b.efficiency_pump[0]
                    == b.design_efficiency * b.motor_efficiency * b.vfd_efficiency
                )

        elif self.config.variable_efficiency is VariableEfficiency.Fixed:
            # Fixed efficiency pump set-up (efficiency is constant at the best efficiency point)
            @self.Constraint(
                doc="Overall efficiency calculation including motor efficiency and VFD efficiency"
            )
            def overall_efficiency_constraint(b):
                return (
                    b.efficiency_pump[0]
                    == b.design_efficiency * b.motor_efficiency * b.vfd_efficiency
                )

        else:
            raise ValueError(
                "Invalid variable efficiency option specified. Must be either VariableEfficiency.Fixed or VariableEfficiency.Flow."
            )

    def initialize_build(
        self,
        state_args=None,
        outlvl=idaeslog.NOTSET,
        solver=None,
        optarg=None,
    ):
        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(self.name, outlvl, tag="unit")
        opt = get_solver(solver, optarg)

        if state_args is None:
            state_args = {}
            state_dict = self.control_volume.properties_in[
                self.flowsheet().config.time.first()
            ].define_port_members()

            for k in state_dict.keys():
                if state_dict[k].is_indexed():
                    state_args[k] = {}
                    for m in state_dict[k].keys():
                        state_args[k][m] = state_dict[k][m].value
                else:
                    state_args[k] = state_dict[k].value

        flags = self.control_volume.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args,
        )

        init_log.info_high("Initialization Step 1 Complete.")

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(self, tee=slc.tee)

        init_log.info_high("Initialization Step 2 {}.".format(idaeslog.condition(res)))

        self.control_volume.release_state(flags, outlvl + 1)
        init_log.info("Initialization Complete: {}".format(idaeslog.condition(res)))

        if not check_optimal_termination(res):
            raise InitializationError(f"Unit model {self.name} failed to initialize")

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()
