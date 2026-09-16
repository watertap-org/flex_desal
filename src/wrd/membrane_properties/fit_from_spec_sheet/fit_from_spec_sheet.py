import matplotlib.pyplot as plt
from pyomo.environ import (
    ConcreteModel,
    check_optimal_termination,
    value,
    units as pyunits,
)

from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core import FlowsheetBlock
import idaes.core.util.scaling as iscale

from watertap.unit_models.reverse_osmosis_1D import (
    ReverseOsmosis1D as RO1D,
    PressureChangeType,
    MassTransferCoefficient,
    ConcentrationPolarizationType,
)
from watertap.property_models import NaCl_T_dep_prop_pack as props
from watertap.core.solvers import get_solver

__author__ = "Kurban Sitterley"

solver = get_solver()

rho = 997.0 * pyunits.kg / pyunits.m**3

bw30_4040 = {
    "feed_conc": 2 * pyunits.g / pyunits.liter,
    "recovery": 0.15,
    "pressure": 225 * pyunits.psi,
    "perm_vol_flow": 9.1 * pyunits.m**3 / pyunits.day,
    "salt_rej": 0.995,
    "temperature": 25,
    "mem_length": (1.016 - 2 * 0.0267) * pyunits.m,
    "mem_area": 7.2 * pyunits.m**2,
    "pressure_loss": -15 * pyunits.psi,
    "spacer_thickness": 34 * 1e-3 * pyunits.inch,  # mil
    "spec_sheet_link": "https://www.lenntech.com/Data-sheets/DuPont-FilmTec-BW30-4040-L.pdf",
    "membrane_name": "BW30-4040",
}

tmg20d_400 = {
    "feed_conc": 2 * pyunits.g / pyunits.liter,
    "recovery": 0.15,
    "pressure": 150 * pyunits.psi,
    "perm_vol_flow": 45.8 * pyunits.m**3 / pyunits.day,
    "salt_rej": 0.997,
    "temperature": 25,
    "mem_length": 1.016 * pyunits.m,
    "mem_area": 37 * pyunits.m**2,
    "pressure_loss": -15 * pyunits.psi,
    "spacer_thickness": 34 * 1e-3 * pyunits.inch,  # mil
    "spec_sheet_link": "https://www.streamlinefiltration.com/wp-content/uploads/2015/05/TMGD.pdf",
    "membrane_name": "TMG20D-400",
}

ag8040f_400 = {
    "feed_conc": 2 * pyunits.g / pyunits.liter,
    "recovery": 0.15,
    "pressure": 225 * pyunits.psi,
    "perm_vol_flow": 41.6 * pyunits.m**3 / pyunits.day,
    "salt_rej": 0.995,
    "temperature": 25,
    "mem_length": 1.016 * pyunits.m,
    "mem_area": 37.2 * pyunits.m**2,
    "pressure_loss": -15 * pyunits.psi,
    "spacer_thickness": 34 * 1e-3 * pyunits.inch,  # mil, assumed, not on spec sheet
    "spec_sheet_link": "https://www.lenntech.com/Data-sheets/Veolia-AG-Series-L.pdf",
    "membrane_name": "AG8040F-400",
}


def estimate_params(
    water_perm=4.2e-12, salt_perm=3.5e-8, porosity=0.95, spec_sheet_params=dict()
):
    """
    Run one iteration of the parameter estimation
    for water permeability, salt permeability, and spacer porosity
    using the provided specification sheet parameters.
    """

    perm_vol_flow = spec_sheet_params.get("perm_vol_flow", None)
    feed_conc = spec_sheet_params.get("feed_conc", None)
    recovery = spec_sheet_params.get("recovery", None)
    pressure = spec_sheet_params.get("pressure", None)
    salt_rej = spec_sheet_params.get("salt_rej", None)
    temperature = spec_sheet_params.get("temperature", None)
    mem_length = spec_sheet_params.get("mem_length", None)
    mem_area = spec_sheet_params.get("mem_area", None)
    pressure_loss = spec_sheet_params.get("pressure_loss", None)
    spacer_thickness = spec_sheet_params.get("spacer_thickness", None)

    if any(
        x is None
        for x in [
            perm_vol_flow,
            feed_conc,
            recovery,
            pressure,
            salt_rej,
            temperature,
            mem_length,
            mem_area,
            pressure_loss,
            spacer_thickness,
        ]
    ):
        raise ValueError(
            "One or more required specification sheet parameters are missing."
        )

    perm_mass_flow = pyunits.convert(
        rho * perm_vol_flow, to_units=pyunits.kg / pyunits.s
    )

    feed_vol_flow = pyunits.convert(
        perm_vol_flow / recovery, to_units=pyunits.m**3 / pyunits.s
    )
    feed_mass_flow_water = value(perm_mass_flow / recovery)
    feed_mass_flow_salt = value(
        pyunits.convert(feed_vol_flow * feed_conc, to_units=pyunits.kg / pyunits.s)
    )

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = props.NaClParameterBlock()

    m.fs.RO = RO1D(
        property_package=m.fs.properties,
        has_pressure_change=True,
        pressure_change_type=PressureChangeType.calculated,
        mass_transfer_coefficient=MassTransferCoefficient.calculated,
        concentration_polarization_type=ConcentrationPolarizationType.calculated,
        transformation_scheme="BACKWARD",
        transformation_method="dae.finite_difference",
        module_type="spiral_wound",
        finite_elements=10,
        has_full_reporting=True,
    )

    print("RO DOF:", degrees_of_freedom(m.fs.RO))

    m.fs.RO.inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].fix(feed_mass_flow_salt)
    m.fs.RO.inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(feed_mass_flow_water)
    m.fs.RO.inlet.pressure[0].fix(pressure)
    m.fs.RO.inlet.temperature[0].fix(temperature + 273.15)

    m.fs.RO.permeate.pressure[0].fix(101325)
    # m.fs.RO.feed_side.channel_height.fix(1e-3)
    m.fs.RO.feed_side.channel_height.fix(spacer_thickness)
    m.fs.RO.length.fix(mem_length)
    m.fs.RO.area.fix(mem_area)

    m.fs.RO.A_comp.fix(water_perm)
    m.fs.RO.B_comp.fix(salt_perm)
    m.fs.RO.feed_side.spacer_porosity.fix(porosity)
    m.fs.RO.feed_side.spacer_porosity.setlb(0.65)

    print("DOF = ", degrees_of_freedom(m))
    print("RO DOF = ", degrees_of_freedom(m.fs.RO))

    iscale.set_scaling_factor(m.fs.RO.area, 1e-2)
    iscale.set_scaling_factor(m.fs.RO.feed_side.area, 1e-2)
    iscale.set_scaling_factor(m.fs.RO.width, 1e-2)

    m.fs.properties.set_default_scaling("flow_mass_phase_comp", 1, index=("Liq", "H2O"))
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e2, index=("Liq", "NaCl")
    )

    iscale.calculate_scaling_factors(m)

    m.fs.RO.initialize()

    results = solver.solve(m)
    if check_optimal_termination(results):
        display_results(m)
    else:
        raise RuntimeError("Solver did not terminate optimally.")

    # Unfix A variable
    print("unfix A, fix perm flow...")
    m.fs.RO.A_comp.unfix()

    # Fix the permeate flow
    m.fs.RO.mixed_permeate[0.0].flow_mass_phase_comp["Liq", "H2O"].fix(perm_mass_flow)

    print("DOF = ", degrees_of_freedom(m))

    results = solver.solve(m)
    if check_optimal_termination(results):
        display_results(m)
    else:
        raise RuntimeError("Solver did not terminate optimally.")

    # Unfix B variable
    print("unfix B, fix rejection...")
    m.fs.RO.B_comp.unfix()

    # Fix the salt rejection
    m.fs.RO.rejection_phase_comp[0, "Liq", "NaCl"].fix(salt_rej)

    print("DOF = ", degrees_of_freedom(m))

    results = solver.solve(m)
    if check_optimal_termination(results):
        display_results(m)
    else:
        raise RuntimeError("Solver did not terminate optimally.")

    # Unfix the permeate flow rate and concentration
    m.fs.RO.mixed_permeate[0.0].flow_mass_phase_comp["Liq", "H2O"].unfix()
    m.fs.RO.rejection_phase_comp[0, "Liq", "NaCl"].unfix()

    # Fix the new A & B values. This will fix them at the solved values from the previous steps
    m.fs.RO.A_comp.fix()
    m.fs.RO.B_comp.fix()

    # Fix the new flow rate
    m.fs.RO.inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(feed_mass_flow_water)
    # Fix the pressure drop
    m.fs.RO.deltaP.fix(pressure_loss)

    # Unfix the spacer porosity
    print("Unfix spacer porosity...")
    m.fs.RO.feed_side.spacer_porosity.unfix()

    print("DOF = ", degrees_of_freedom(m))

    results = solver.solve(m)
    if check_optimal_termination(results):
        display_results(m)
    else:
        raise RuntimeError("Solver did not terminate optimally.")

    water_perm = m.fs.RO.A_comp[0, "H2O"].value
    salt_perm = m.fs.RO.B_comp[0, "NaCl"].value
    porosity = m.fs.RO.feed_side.spacer_porosity.value

    return water_perm, salt_perm, porosity


def iterate_estimation(spec_sheet_params=dict(), num_iter=3, close_figs=False):
    """
    Iterate the estimation of membrane parameters (water permeability, salt permeability, and spacer porosity)
    for a specified number of iterations using the provided specification sheet parameters.
    """

    wp = list()
    sp = list()
    por = list()
    num = list()

    # Initial values for estimation routine
    water_perm = 4.2e-12
    salt_perm = 3.5e-8
    porosity = 0.95

    for n in range(1, num_iter + 1):
        water_perm, salt_perm, porosity = estimate_params(
            water_perm=water_perm,
            salt_perm=salt_perm,
            porosity=porosity,
            spec_sheet_params=spec_sheet_params,
        )
        wp.append(water_perm)
        sp.append(salt_perm)
        por.append(porosity)
        num.append(n)

    f1, ax1 = plt.subplots(figsize=(3, 3))
    ax1.scatter(num, wp, marker="o", color="r")
    ax1.set_title(f"{spec_sheet_params['membrane_name']}\nWater Permeability")

    f2, ax2 = plt.subplots(figsize=(3, 3))
    ax2.scatter(num, sp, marker="o", color="g")
    ax2.set_title(f"{spec_sheet_params['membrane_name']}\nSalt Permeability")

    f3, ax3 = plt.subplots(figsize=(3, 3))
    ax3.scatter(num, por, marker="o", color="k")
    ax3.set_title(f"{spec_sheet_params['membrane_name']}\nPorosity")

    for f, a in [(f1, ax1), (f2, ax2), (f3, ax3)]:
        a.label_outer()
        a.set_axisbelow(True)
        a.set_xlabel("Iteration")
        a.grid(visible=True)
        f.tight_layout()

    if close_figs:
        plt.close("all")

    # return final values
    return water_perm, salt_perm, porosity


def display_results(m):

    print("\n--------- OPTIMAL SOLVE!!! ---------\n")
    print(f"water_perm = {m.fs.RO.A_comp[0,'H2O']():.3e} m/s/Pa")
    print(f"salt_perm = {m.fs.RO.B_comp[0,'NaCl']():.3e} m/s")
    print(f"porosity = {m.fs.RO.feed_side.spacer_porosity():.3f}")
    print("\n")


if __name__ == "__main__":

    water_perm, salt_perm, porosity = iterate_estimation(
        spec_sheet_params=bw30_4040, close_figs=True
    )
