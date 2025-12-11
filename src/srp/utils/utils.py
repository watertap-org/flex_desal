from idaes.models.unit_models import Mixer, Separator


__all__ = [
    "touch_flow_and_conc",
]


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
