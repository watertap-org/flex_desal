from pyomo.environ import (
    Suffix,
    check_optimal_termination,
)
from watertap.core.solvers import get_solver

from pyomo.common.config import ConfigValue, In
from pyomo.environ import Var, units as pyunits, Expr_if, value

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

# Unclear
# How to bring in head vs flow curve surrogate?
# How to bring in efficiency vs flow curve surrogate?

# Add VFD and motor efficiencies

class VariableEfficiency(Enum):
    Fixed = auto()  # default is constant efficiency
    Flow = auto()  # flow-only correlation
    FlowHead = auto()  # flow and head correlation


@declare_process_block_class("Pump")
class PumpIsothermalData(InitializationMixin, PumpData):
    """
    Detailed Isothermal Pump Unit Model Class

    * User can specify whether to use a variable efficiency based on flow and head, or to use a constant efficiency at the best efficiency point (BEP).
    * This function builds on the standard pump model by adding variables and constraints to represent pump performance curves. 
    * The model includes options for defining the pump curve based on flow and head at the BEP, as well as scaling the efficiency curve accordingly. 
    * This allows for more accurate representation of pump performance under varying operating conditions.
    
    """

    CONFIG = PumpData.CONFIG()

    CONFIG.declare(
        "variable_efficiency",
        ConfigValue(
            default = VariableEfficiency.Fixed,
            domain= In(VariableEfficiency),
            description="Variable pump efficiency flag",
            doc= '''Indicates whether pump efficiency should be variable or constant
                ** Fixed: The default is constant efficiency at the best efficiency point (BEP)
                ** Flow: Efficiency is a function of flow only, scaled to the BEP flow rate
                ** FlowHead: Efficiency is a function of both flow and head, scaled to the BEP flow rate and head
            ''',
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

        if self.config.variable_efficiency is not VariableEfficiency.Fixed:
            # Variable efficiency pump set-up
            #### Design point variables ####
            self.design_flow = Var(
                initialize=1.0,
                doc='''Design flowrate of the centrifugal pump. 
                This could be the flowrate at the best efficiency point (BEP) or another reference point.
                Used to build the system curve.''',
                units=pyunits.m**3 / pyunits.s,
            )

            self.design_head = Var(
                initialize=1.0,
                doc='''Design head of the centrifugal pump. 
                This could be the head at the best efficiency point (BEP) or another reference point.
                Used to build the system curve.''',
                units=pyunits.m,
            )

            # self.design_efficiency = Var(
            #     initialize=0.8,
            #     doc='''Design efficiency of the centrifugal pump. 
            #     This could be the efficiency at the best efficiency point (BEP) or another reference point.
            #     Used to scale the efficiency curve.''',
            #     units=pyunits.dimensionless,
            # )

            self.design_speed_fraction = Var(
                initialize=0.8,
                doc='''Design speed fraction of the centrifugal pump.
                This could be the speed fraction at the best efficiency point (BEP) or another reference point.
                Used to scale the efficiency curve.''',
                units=pyunits.dimensionless,
            )

            ### System curve variables ###
            # Head = system_curve_geometric_head +  system_curve_flow_constant * (flow)**2

            self.system_curve_geometric_head = Var(
                initialize=0.0,
                doc='''Geometric head constant for the pump, used to define the system curve.''',
                units=pyunits.m,
            )

            self.system_curve_flow_constant = Var(
                initialize=1.0,
                doc='''Geometric flow constant for the pump, used to define the system curve.''',
                units= pyunits.m * (pyunits.m**3 / pyunits.s)**(-2),
            )

            # Constraints connecting inlet and outlet conditions to the design point head and flow, used to solve for the system curve constants
            @self.Constraint(doc="Design head is the pressure difference across the pump at the design point")
            def design_head_constraint(b):
                return b.design_head == b.control_volume.deltaP[0] / (b.control_volume.properties_out[0].dens_mass_phase['Liq'] * 9.81 * pyunits.m / pyunits.s**2)
            
            @self.Constraint(doc="Design flow is the flow through the pump at the design point")
            def design_flow_constraint(b):
                return b.design_flow == pyunits.convert(b.control_volume.properties_out[0].flow_vol_phase['Liq'], to_units=pyunits.m**3 / pyunits.s)


            # Constraints to calculate system curve constants based on design point
            @self.Constraint(doc="System curve head constant calculation")
            def system_curve_head_constant_calculation(b):
                return b.design_head == b.system_curve_geometric_head +  b.system_curve_flow_constant * (b.design_flow**2)
            
            #### Pump curve variables ####
            self.ref_flow = Var(
                initialize=1.0,
                doc= "Reference flowrate for the pump, used to scale the efficiency curve.",
                units= pyunits.m**3 / pyunits.s,
            )

            self.ref_head = Var(
                initialize=10.0,
                doc="Reference head for the pump, used to scale the efficiency curve.",
                units=pyunits.m,
            )

            self.ref_efficiency = Var(
                initialize=0.8,
                doc="Reference efficiency for the pump, used to scale the efficiency curve.",
                units=pyunits.dimensionless,
            )

            self.ref_speed_fraction = Var(
                initialize=1,
                doc="Reference speed fraction for the pump, used to scale the efficiency curve.",
                units=pyunits.dimensionless,
            )

            # Constraints to calculate reference point variables based on design point
            @self.Constraint(doc="Reference head equality")
            def ref_head_constraint(b):
                return b.ref_head ==  b.system_curve_geometric_head + b.system_curve_flow_constant * (b.ref_flow**2)
            
            #TODO: Update to calculate the surrogate coefficients based on the reference point variables, rather than hardcoding the surrogate coefficients
            @self.Constraint(doc="Reference head calculation calculated using the pump curve surrogate")
            def ref_head_surrogate(b):
                return b.ref_head == -8089.1 * (b.ref_flow**3) + 2729.2 * b.ref_flow**2 + - 410.6 * b.ref_flow + 114.22
            
            @self.Constraint(doc="Reference efficiency calculated calculated using the pump curve surrogate")
            def ref_efficiency_constraint(b):
                return b.ref_efficiency == -138.82 * (b.ref_flow**3) + 41.373 * b.ref_flow**2 + -0.535* b.ref_flow +  0.389
            
            @self.Constraint(doc="Design speed fraction calculation using affinity laws")
            def design_speed_fraction_calculation(b):
                return b.design_speed_fraction == (b.design_flow / b.ref_flow) * b.ref_speed_fraction
            
        
            # Expression to calculate the efficiency at the design point based as a function of speed fraction
            @self.Constraint(doc="Design efficiency calculation")
            def design_efficiency_constraint(b):
                return b.efficiency_pump[0] == b.ref_efficiency * b.design_speed_fraction**0.1
            
    
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

        if iscale.get_scaling_factor(self.ref_flow) is None:
            iscale.set_scaling_factor(self.ref_flow, 1e2)
        
        





    