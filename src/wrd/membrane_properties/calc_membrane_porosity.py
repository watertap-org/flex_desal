import pandas as pd
from wrd.membrane_properties.ro_for_membrane_porosity import solve_ro_module
from pyomo.environ import (
    assert_optimal_termination,
    units as pyunits,
    value,
)
from watertap.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom

def report_membrane_properties(m,solve_num):
    print(f"--- SOLVE {solve_num} ---")
    print(f"Calculated membrane porosity: {value(m.fs.ro.unit.feed_side.spacer_porosity)}")
    print(f"Calculated water permeability (A): {value(m.fs.ro.unit.A_comp[0,'H2O'])} m/s/Pa")
    print(f"Calculated salt permeability (B): {value(m.fs.ro.unit.B_comp[0,'NaCl'])} m/s")
    print(f"Calculated pressure drop across membrane: {value(pyunits.convert(m.fs.ro.unit.deltaP[0], to_units=pyunits.psi))} psi")
 
def calc_membrane_porosity(Qin,Qperm,Cin,Cperm,Pin,Pout):
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

    # Intialize. Then re-calculate unfixing A 
    stage_RR = Qperm / Qin
    m.fs.ro.unit.feed_side.spacer_porosity.unfix()
    m.fs.ro.unit.A_comp.unfix()
    m.fs.ro.unit.B_comp.unfix()

    m.fs.ro.unit.recovery_vol_phase[0, "Liq"].fix(stage_RR)
    m.fs.ro.unit.mixed_permeate[0].conc_mass_phase_comp["Liq", "NaCl"].fix(Cperm)

    # assert degrees_of_freedom(m) == 0
    solver = get_solver()
    results = solver.solve(m)
    assert_optimal_termination(results)
    report_membrane_properties(m, solve_num=1)
    # Get new A and B values
    A_initial = m.fs.ro.unit.A_comp[0, "H2O"].value
    B_initial = m.fs.ro.unit.B_comp[0, "NaCl"].value
    porosity_initial = m.fs.ro.unit.feed_side.spacer_porosity.value
    
    # Now calculate porosity using the A and B values and the equations from the membrane properties document
    m.fs.ro.unit.A_comp.fix(A_initial)
    m.fs.ro.unit.B_comp.unfix()
    m.fs.ro.unit.feed_side.spacer_porosity.unfix()
    # m.fs.ro.unit.recovery_vol_phase[0, "Liq"].unfix()
    # m.fs.ro.unit.mixed_permeate[0].conc_mass_phase_comp["Liq", "NaCl"].unfix()
    # m.fs.ro.unit.deltaP.unfix() 

    # assert degrees_of_freedom(m) == 0
    results = solver.solve(m)
    assert_optimal_termination(results)
    report_membrane_properties(m, solve_num=2)

    # Now calculate porosity using the A and B values and the equations from the membrane properties document
    m.fs.ro.unit.A_comp.unfix()
    m.fs.ro.unit.B_comp.fix(B_initial)
    m.fs.ro.unit.feed_side.spacer_porosity.unfix()
    # m.fs.ro.unit.deltaP.unfix() 

    # assert degrees_of_freedom(m) == 0
    results = solver.solve(m)
    assert_optimal_termination(results)
    report_membrane_properties(m, solve_num=3)


    # Now solve with unfixed A 
    m.fs.ro.unit.A_comp.unfix()
    m.fs.ro.unit.B_comp.unfix()
    m.fs.ro.unit.feed_side.spacer_porosity.fix(porosity_initial)

    # assert degrees_of_freedom(m) == 0
    results = solver.solve(m)
    assert_optimal_termination(results)
    report_membrane_properties(m, solve_num=4)

    return m


if __name__ == "__main__":
    # A and B vals calculated from August 2021 data with a porosity of 0.7211
    # A = 7.32e-12 * pyunits.m / pyunits.s / pyunits.Pa
    # B = 7.50e-8 * pyunits.m / pyunits.s
    # Test Conditions from the membrane specification sheet
    Qin = 56 # gpm      # 12100/.15 * pyunits.gal / pyunits.day
    Qperm = 8.4 # gpm   # 12100 * pyunits.gal / pyunits.day
    Cin = 2.0 # g/L
    Cperm = 0.006 # 99.7 % salt rejection
    Pin = 150 # psi
    Pout = 135 # psi # No idea what the exit pressure is... max pressure drop across the membrane is 15 psi specified in the datasheet.

    m = calc_membrane_porosity(Qin=Qin,Qperm=Qperm,Cin=Cin,Cperm=Cperm,Pin=Pin,Pout=Pout)
 
