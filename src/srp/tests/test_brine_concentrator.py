import pytest

from pyomo.environ import value, units as pyunits

import srp.components.brine_concentrator as bc


@pytest.mark.component
def test_brine_concentrator_95_recov():
    m = bc.main(recovery_vol=0.95)

    assert pytest.approx(value(m.fs.costing.LCOW), rel=1e-3) == 4.7559
    assert pytest.approx(value(m.fs.costing.SEC), rel=1e-3) == 23.8141
    assert pytest.approx(value(m.fs.bc.compressor.pressure_ratio), rel=1e-3) == 1.67
    assert pytest.approx(value(m.fs.bc.evaporator.area), rel=1e-3) == 641.9
    assert pytest.approx(value(m.fs.bc.evaporator.lmtd), rel=1e-3) == 25.53
    assert pytest.approx(value(m.fs.bc.hx_brine.area), rel=1e-3) == 6.162
    assert pytest.approx(value(m.fs.bc.hx_distillate.area), rel=1e-3) == 80.92


@pytest.mark.component
def test_brine_concentrator_50_recov():
    m = bc.main(recovery_vol=0.5)

    assert pytest.approx(value(m.fs.costing.LCOW), rel=1e-3) == 4.7846
    assert pytest.approx(value(m.fs.costing.SEC), rel=1e-3) == 22.298
    assert pytest.approx(value(m.fs.bc.compressor.pressure_ratio), rel=1e-3) == 1.612
    assert pytest.approx(value(m.fs.bc.evaporator.area), rel=1e-3) == 322.51
    assert pytest.approx(value(m.fs.bc.evaporator.lmtd), rel=1e-3) == 26.66
    assert pytest.approx(value(m.fs.bc.hx_brine.area), rel=1e-3) == 98.54
    assert pytest.approx(value(m.fs.bc.hx_distillate.area), rel=1e-3) == 107.67
