import pytest
from pyomo.environ import value, units as pyunits
import pandas as pd
from idaes.core.util.exceptions import InitializationError
from wrd.components.pump import main
from wrd.surrogate_models.pump_param_sweep import *


# Add pump tests to test robustness of soving
@pytest.mark.component
def test_pump_UF():
    # This test should cover a range of flows and heads that could cause an issue with solving
    test_points = pd.DataFrame(
        [[4690, 145], [2000, 75], [2600, 150]], columns=["flow", "head"]
    )
    dataset = pump_param_sweep(test_pairs=test_points, pump_type="UF")
    expected_effs = [66.921, 75.834, 74.9799]
    for i in range(len(test_points)):
        assert dataset["total_efficiency"][i] == pytest.approx(
            expected_effs[i], rel=1e-3
        )


@pytest.mark.component
def test_pump_speed_too_high():
    test_point = (3000, 280)
    # Check that InitializationError is raised with expected message pattern
    with pytest.raises(InitializationError) as exc_info:
        m = main(
            Qin=test_point[0],
            head=test_point[1],
            Cin=0,
            Tin=298,
            stage_num=1,
        )
    assert (
        exc_info.value.args[0]
        == "Pump speed ratio too high during initialization: 1.02. Check head and flow inputs."
    )

@pytest.mark.component
def test_speed_and_head_inputs():
    m = main(head=250, speed =.95, stage_num=1, Pin=35.4)
    expected_flow = 2485
    assert value(pyunits.convert(m.fs.pump.feed.properties[0].flow_vol_phase["Liq"],to_units=pyunits.gallon / pyunits.minute)) == pytest.approx(expected_flow, rel=1e-3)