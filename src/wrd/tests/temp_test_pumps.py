import pytest
from pyomo.environ import value, units as pyunits
from wrd.components.pump import main


def test_pump_param_sweep():
    # RO Feed Pump
    # Flow and Head pairs
    test_pairs = [(2675,276.4)]
    for flow, head in test_pairs:
        m = main(Qin= flow,
            head = head,
            Cin=0.5,
            Tin=302,
            Pin=14.5, #psi
            stage_num=1,
            uf=False,
            file="wrd_inputs_8_19_21.yaml",
            add_costing=True,
            )
        speed = value(m.fs.pump.unit.eff.speed)
        print(f"Pump speed for flow {flow} and head {head}: {speed}")
        rpm = value(m.fs.pump.unit.eff.speed * 1780)
        print(f"Pump RPM: {rpm}")
        fluid_efficiency = value(m.fs.pump.unit.eff.efficiency_fluid)
        print(f"Fluid efficiency: {fluid_efficiency}")



if __name__ == "__main__":
    test_pump_param_sweep()