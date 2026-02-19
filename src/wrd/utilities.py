import yaml
import os
from pyomo.environ import value, units as pyunits
from idaes.models.unit_models import Mixer, Separator


__all__ = [
    "report_head_loss",
    "report_sj",
    "report_separator",
    "report_mixer",
    "get_config_value",
    "load_config",
    "get_config_file",
]


def report_head_loss(hl, w=25):

    title = hl.name.replace("fs.", "").replace("_", " ").upper()
    cv = hl.control_volume

    side = int(((3 * w) - len(title)) / 2) - 1
    header = "=" * side + f" {title} " + "=" * side
    print(f"\n{header}\n")
    dp = value(
        pyunits.convert(
            cv.deltaP[0],
            to_units=pyunits.psi,
        )
    )
    print(f'{"Pressure Drop":<{w}s}{f"{dp:<{w}.1f}"}{"psi":<{w}s}')
    flow_in = value(
        pyunits.convert(
            cv.properties_in[0].flow_vol_phase["Liq"],
            to_units=pyunits.gallons / pyunits.minute,
        )
    )
    conc_in = value(
        pyunits.convert(
            cv.properties_in[0].conc_mass_phase_comp["Liq", "NaCl"],
            to_units=pyunits.mg / pyunits.L,
        )
    )
    p_in = value(
        pyunits.convert(
            cv.properties_in[0].pressure,
            to_units=pyunits.psi,
        )
    )
    print(f'{"INLET Flow":<{w}s}{f"{flow_in:<{w},.1f}"}{"gpm":<{w}s}')
    print(f'{"INLET NaCl":<{w}s}{f"{conc_in:<{w},.1f}"}{"mg/L":<{w}s}')
    print(f'{"INLET Pressure":<{w}s}{f"{p_in:<{w},.1f}"}{"psi":<{w}s}')

    flow_out = value(
        pyunits.convert(
            cv.properties_out[0].flow_vol_phase["Liq"],
            to_units=pyunits.gallons / pyunits.minute,
        )
    )
    conc_out = value(
        pyunits.convert(
            cv.properties_out[0].conc_mass_phase_comp["Liq", "NaCl"],
            to_units=pyunits.mg / pyunits.L,
        )
    )
    p_out = value(
        pyunits.convert(
            cv.properties_out[0].pressure,
            to_units=pyunits.psi,
        )
    )
    print(f'{"OUTLET Flow":<{w}s}{f"{flow_out:<{w},.1f}"}{"gpm":<{w}s}')
    print(f'{"OUTLET NaCl":<{w}s}{f"{conc_out:<{w},.1f}"}{"mg/L":<{w}s}')
    print(f'{"OUTLET Pressure":<{w}s}{f"{p_out:<{w},.1f}"}{"psi":<{w}s}')


def report_sj(sj, w=25):

    title = sj.name.replace("fs.", "").replace("_", " ").upper()

    side = int(((3 * w) - len(title)) / 2) - 1
    header = "=" * side + f" {title} " + "=" * side
    print(f"\n{header}\n")
    flow_in = value(
        pyunits.convert(
            sj.properties[0].flow_vol_phase["Liq"],
            to_units=pyunits.gallons / pyunits.minute,
        )
    )
    conc_in = value(
        pyunits.convert(
            sj.properties[0].conc_mass_phase_comp["Liq", "NaCl"],
            to_units=pyunits.mg / pyunits.L,
        )
    )
    print(f'{"INLET/OUTLET Flow":<{w}s}{f"{flow_in:<{w},.1f}"}{"gpm":<{w}s}')
    print(f'{"INLET/OUTLET NaCl":<{w}s}{f"{conc_in:<{w},.1f}"}{"mg/L":<{w}s}')


def report_mixer(mixer, w=25):

    ms = mixer.find_component("mixed_state")
    title = mixer.name.replace("fs.", "").replace("_", " ").upper()
    side = int(((3 * w) - len(title)) / 2) - 1
    header = "=" * side + f" {title} " + "=" * side
    print(f"\n{header}\n")
    tot_flow_in = sum(
        value(
            pyunits.convert(
                mixer.find_component(f"{x}_state")[0].flow_vol_phase["Liq"],
                to_units=pyunits.gallons / pyunits.minute,
            )
        )
        for x in mixer.config.inlet_list
    )
    print(f'{"TOTAL INLET FLOW":<{w}s}{f"{tot_flow_in:<{w},.1f}"}{"gpm":<{w}s}')
    for x in mixer.config.inlet_list:
        sb = mixer.find_component(f"{x}_state")
        flow_in = value(
            pyunits.convert(
                sb[0].flow_vol_phase["Liq"],
                to_units=pyunits.gallons / pyunits.minute,
            )
        )
        conc_in = value(
            pyunits.convert(
                sb[0].conc_mass_phase_comp["Liq", "NaCl"],
                to_units=pyunits.mg / pyunits.L,
            )
        )
        p_in = value(
            pyunits.convert(
                sb[0].pressure,
                to_units=pyunits.psi,
            )
        )
        print(
            f'{"   Flow " + x.replace("_", " ").title():<{w}s}{f"{flow_in:<{w},.1f}"}{"gpm":<{w}s}'
        )
        print(
            f'{"   NaCl " + x.replace("_", " ").title():<{w}s}{f"{conc_in:<{w},.1f}"}{"mg/L":<{w}s}'
        )
        print(
            f'{"   Pressure " + x.replace("_", " ").title():<{w}s}{f"{p_in:<{w},.1f}"}{"psi":<{w}s}'
        )
    flow_out = value(
        pyunits.convert(
            ms[0].flow_vol_phase["Liq"],
            to_units=pyunits.gallons / pyunits.minute,
        )
    )
    conc_out = value(
        pyunits.convert(
            ms[0].conc_mass_phase_comp["Liq", "NaCl"],
            to_units=pyunits.mg / pyunits.L,
        )
    )
    p_out = value(
        pyunits.convert(
            ms[0].pressure,
            to_units=pyunits.psi,
        )
    )
    print(f'{"Outlet Flow":<{w}s}{f"{flow_out:<{w},.1f}"}{"gpm":<{w}s}')
    print(f'{"Outlet NaCl":<{w}s}{f"{conc_out:<{w},.1f}"}{"mg/L":<{w}s}')
    print(f'{"Outlet Pressure":<{w}s}{f"{p_out:<{w},.1f}"}{"psi":<{w}s}')


def report_separator(sep, w=25):

    title = sep.name.replace("fs.", "").replace("_", " ").upper()

    side = int(((3 * w) - len(title)) / 2) - 1
    header = "=" * side + f" {title} " + "=" * side
    print(f"\n{header}\n")
    ms = sep.find_component("mixed_state")
    flow_in = value(
        pyunits.convert(
            ms[0].flow_vol_phase["Liq"],
            to_units=pyunits.gallons / pyunits.minute,
        )
    )
    conc_in = value(
        pyunits.convert(
            ms[0].conc_mass_phase_comp["Liq", "NaCl"],
            to_units=pyunits.mg / pyunits.L,
        )
    )
    print(f'{"INLET Flow":<{w}s}{f"{flow_in:<{w},.1f}"}{"gpm":<{w}s}')
    print(f'{"INLET NaCl":<{w}s}{f"{conc_in:<{w},.1f}"}{"mg/L":<{w}s}')
    tot_flow_out = sum(
        value(
            value(
                pyunits.convert(
                    sep.find_component(f"{x}_state")[0].flow_vol_phase["Liq"],
                    to_units=pyunits.gallons / pyunits.minute,
                )
            )
        )
        for x in sep.config.outlet_list
    )
    print(f'{"TOTAL OUTLET FLOW":<{w}s}{f"{tot_flow_out:<{w},.1f}"}{"gpm":<{w}s}')
    for x in sep.config.outlet_list:
        sb = sep.find_component(f"{x}_state")
        flow_out = value(
            pyunits.convert(
                sb[0].flow_vol_phase["Liq"],
                to_units=pyunits.gallons / pyunits.minute,
            )
        )
        conc_out = value(
            pyunits.convert(
                sb[0].conc_mass_phase_comp["Liq", "NaCl"],
                to_units=pyunits.mg / pyunits.L,
            )
        )
        print(
            f'{"   Flow " + x.replace("_", " ").title():<{w}s}{f"{flow_out:<{w},.1f}"}{"gpm":<{w}s}'
        )
        print(
            f'{"   NaCl " + x.replace("_", " ").title():<{w}s}{f"{conc_out:<{w},.1f}"}{"mg/L":<{w}s}'
        )


def load_config(config):
    with open(config, "r") as file:
        return yaml.safe_load(file)


def get_config_value(
    config,
    key,
    section,
    subsection=None,
):
    """
    Get a value from the configuration file.
    """

    if section in config:
        if subsection:
            if subsection in config[section]:
                if key in config[section][subsection]:
                    if (
                        isinstance(config[section][subsection][key], dict)
                        and "value" in config[section][subsection][key]
                        and "units" in config[section][subsection][key]
                    ):
                        return config[section][subsection][key]["value"] * getattr(
                            pyunits, config[section][subsection][key]["units"]
                        )
                    return config[section][subsection][key]
                else:
                    raise KeyError(
                        f"Key '{key}' not found in subsection '{subsection}' of section '{section}' of the configuration."
                    )
            else:
                raise KeyError(
                    f"Section '{section}' or subsection '{subsection}' not found in the configuration."
                )
        else:
            if key in config[section]:
                if (
                    isinstance(config[section][key], dict)
                    and "value" in config[section][key]
                    and "units" in config[section][key]
                ):
                    return config[section][key]["value"] * getattr(
                        pyunits, config[section][key]["units"]
                    )
                return config[section][key]
            else:
                raise KeyError(
                    f"Key '{key}' not found in section '{section}' of the configuration."
                )
    else:
        raise KeyError(f"Section '{section}' not found in the configuration.")


def get_config_file(yaml_name):
    # WRD RO configurations input file. References to all values included in yml file
    # Get the absolute path of the current script
    current_script_path = os.path.abspath(__file__)
    # Get the directory containing the current script
    current_directory = os.path.dirname(current_script_path)
    # Get the parent directory of the current directory (one folder prior)
    # parent_directory = os.path.dirname(current_directory)
    config_file_name = os.path.join(current_directory, "meta_data", yaml_name)
    return config_file_name


def get_chem_list(yaml_name, section):
    # Section must be "pre_treatment" or "post_treatment"
    chem_list = []
    config_file_name = get_config_file(yaml_name)
    config = load_config(config_file_name)
    if section in config:
        for subsection in config[section]:
            chem_list.append(subsection)
        # chem_list.remove("default")
    return chem_list


def touch_flow_and_conc(b):
    """
    Touch flow and conc variables for construction
    """
    props = b.find_component("properties")
    if props is not None:
        props[0].flow_vol_phase
        props[0].conc_mass_phase_comp
        return

    cv = b.find_component("control_volume")
    if cv is not None:
        cv.properties_in[0].flow_vol_phase
        cv.properties_in[0].conc_mass_phase_comp
        cv.properties_out[0].flow_vol_phase
        cv.properties_out[0].conc_mass_phase_comp
        return

    b.mixed_state[0].flow_vol_phase
    b.mixed_state[0].conc_mass_phase_comp

    if isinstance(b, Mixer):
        for x in b.config.inlet_list:
            sb = b.find_component(f"{x}_state")
            sb[0].flow_vol_phase
            sb[0].conc_mass_phase_comp
    if isinstance(b, Separator):
        for x in b.config.outlet_list:
            sb = b.find_component(f"{x}_state")
            sb[0].flow_vol_phase
            sb[0].conc_mass_phase_comp


def ft_head_to_psi(head):
    """
    Convert head (in feet) to pressure (in psi)
    Default fluid density is for water at room temperature (998.2 kg/m^3)
    Default gravity is 9.81 m/s^2
    """
    # head should have units of ft already
    fluid_density = 62.4 * pyunits.lb / pyunits.ft**3
    gravity = 32.174 * pyunits.ft / pyunits.s**2
    pressure = fluid_density * gravity * head
    return pyunits.convert(pressure, to_units=pyunits.psi)


def psi_to_ft_head(pressure):
    """
    Convert pressure (in psi) to head (in feet)
    Default fluid density is for water at room temperature (998.2 kg/m^3)
    Default gravity is 9.81 m/s^2
    """
    fluid_density = 62.4 * pyunits.lb / pyunits.ft**3
    gravity = 32.174 * pyunits.ft / pyunits.s**2
    head = pyunits.convert(pressure / (fluid_density * gravity), to_units=pyunits.ft)
    return head
