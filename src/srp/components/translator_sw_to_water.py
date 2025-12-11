from pyomo.environ import check_optimal_termination

from idaes.core import declare_process_block_class
from idaes.models.unit_models.translator import TranslatorData
from idaes.core.util.model_statistics import degrees_of_freedom
import idaes.logger as idaeslog
from idaes.core.util.exceptions import InitializationError

from watertap.core.solvers import get_solver


__author__ = "Kurban Sitterley"


@declare_process_block_class("TranslatorSWtoWater")
class TranslatorSWtoWaterData(TranslatorData):
    """
    Translator block for Seawater to Water property packages
    """

    CONFIG = TranslatorData.CONFIG()

    def build(self):

        super().build()

        @self.Constraint(doc="Isothermal")
        def eq_temperature(b):
            return b.properties_in[0].temperature == b.properties_out[0].temperature

        @self.Constraint(doc="Isobaric")
        def eq_pressure(b):
            return b.properties_in[0].pressure == b.properties_out[0].pressure

        @self.Constraint(
            doc="Equality mass flow water equation",
        )
        def eq_flow_mass_water(b):
            return (
                b.properties_in[0].flow_mass_phase_comp["Liq", "H2O"]
                == b.properties_out[0].flow_mass_phase_comp["Liq", "H2O"]
            )

    def initialize_build(
        self,
        state_args_in=None,
        state_args_out=None,
        outlvl=idaeslog.NOTSET,
        solver=None,
        optarg=None,
    ):
        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="unit")

        # Create solver
        opt = get_solver(solver, optarg)

        # ---------------------------------------------------------------------
        # Initialize state block
        flags = self.properties_in.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args_in,
            hold_state=True,
        )

        self.properties_out.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args_out,
        )

        if degrees_of_freedom(self) != 0:
            raise Exception(
                f"{self.name} degrees of freedom were not 0 at the beginning "
                f"of initialization. DoF = {degrees_of_freedom(self)}"
            )

        with idaeslog.solver_log(init_log, idaeslog.DEBUG) as slc:
            res = opt.solve(self, tee=slc.tee)

        self.properties_in.release_state(flags=flags, outlvl=outlvl)

        init_log.info(f"Initialization Complete: {idaeslog.condition(res)}")

        if not check_optimal_termination(res):
            raise InitializationError(
                f"{self.name} failed to initialize successfully. Please check "
                f"the output logs for more information."
            )
