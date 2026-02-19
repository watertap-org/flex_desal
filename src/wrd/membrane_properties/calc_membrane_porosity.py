import pandas as pd
from wrd.membrane_properties.ro_for_membrane_porosity import solve_ro_module
from pyomo.environ import (
    assert_optimal_termination,
    units as pyunits,
    value,
)
from watertap.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom

 
def calc_membrane_porosity(A,B,Qin,Qperm,Cin,Cperm,Pin,Pout):
    """
    This function calculates membrane permeability values (A and B) for the specified RO stage.
    """
    m = solve_ro_module(
        Qin=Qin,
        Cin=Cin,
        Tin=298,
        Pin=Pin,
        Pout=Pout,
        stage_num=1
    )

    # Intialize, then unfix membrane porosity and fix the recovery.
    stage_RR = Qperm / Qin
    m.fs.ro.unit.feed_side.spacer_porosity.unfix()
    m.fs.ro.unit.A_comp.fix(A)
    m.fs.ro.unit.B_comp.fix(B)

    m.fs.ro.unit.recovery_vol_phase[0, "Liq"].fix(stage_RR)
    m.fs.ro.unit.mixed_permeate[0].conc_mass_phase_comp["Liq", "NaCl"].fix(Cperm)

    # Then apply the new A and B values

    assert degrees_of_freedom(m) == 0
    solver = get_solver()
    results = solver.solve(m)
    assert_optimal_termination(results)

    return m


if __name__ == "__main__":
    # A and B vals calculated from August 2021 data with a porosity of 0.7211
    A = 7.32e-12 * pyunits.m / pyunits.s / pyunits.Pa
    B = 7.50e-8 * pyunits.m / pyunits.s
    # Test Conditions from the membrane specification sheet
    Qin = 56 # gpm      # 12100/.15 * pyunits.gal / pyunits.day
    Qperm = 8.4 # gpm   # 12100 * pyunits.gal / pyunits.day
    Cin = 2.0 # g/L
    Cperm = 0.006 # 99.7 % salt rejection
    Pin = 150 # psi
    Pout = 145 # psi # No idea what the exit pressure is... This is the meximum pressure drop across the membrane specified in the datasheet.

    m = calc_membrane_porosity(A,B,Qin=Qin,Qperm=Qperm,Cin=Cin,Cperm=Cperm,Pin=Pin,Pout=Pout)
 
    print(f"Calculated membrane porosity: {value(m.fs.ro.unit.feed_side.spacer_porosity)}")
    print(f"Calculated permeate salinity: {value( m.fs.ro.unit.mixed_permeate[0].conc_mass_phase_comp['Liq', 'NaCl'] )}")
