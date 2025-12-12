import yaml
import os
from pyomo.environ import units as pyunits
from idaes.models.unit_models import Mixer, Separator

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