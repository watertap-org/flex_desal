import pytest
from pyomo.environ import value, units as pyunits
from wrd.components.pump import main


@pytest.mark.component
def test_pump_main():

    # August 19, 2021 Data
    # Stage 1 is default
    m = main()
    assert pytest.approx(value(m.fs.costing.SEC), rel=1e-3) == 0.33924  # kWh/m3
    # Stage 2
    m = main(Qin=1029, Pin=131.2 * pyunits.psi, stage_num=2)
    assert pytest.approx(value(m.fs.costing.SEC), rel=1e-3) == 0.07753
    # Stage 3
    m = main(Qin=384, Pin=(112.6 - 41.9) * pyunits.psi, stage_num=3)
    assert pytest.approx(value(m.fs.costing.SEC), rel=1e-3) == 0.40512

# TODO: Add tests for different pump speeds and compare to both power use and head. 

@pytest.mark.component
def test_pump_PRO_stage2_2_20_21():
    m = main(
    Qin=983.2,
    Cin=1.2,
    Tin=302,
    Pin=143.3 * pyunits.psi,
    speed=.756,
    stage_num=2,
    file="wrd_inputs_2_20_21.yaml",
)
    assert pytest.approx(value(pyunits.convert(m.fs.pump.unit.work_mechanical[0], to_units=pyunits.kW)), rel=.15) == 9.7