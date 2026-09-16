import pytest
from wrd.membrane_properties.fit_from_spec_sheet import fit_from_spec_sheet as fit


@pytest.mark.unit
def test_no_spec_sheet_data():
    with pytest.raises(
        ValueError,
        match="One or more required specification sheet parameters are missing.",
    ):
        fit.estimate_params(spec_sheet_params={})


@pytest.mark.component
def test_estimate_params():
    for d in [fit.bw30_4040, fit.ag8040f_400, fit.tmg20d_400]:

        water_perm, salt_perm, porosity = fit.estimate_params(
            water_perm=4.2e-12,
            salt_perm=3.5e-8,
            porosity=0.95,
            spec_sheet_params=d,
        )

        if d["membrane_name"] == "BW30-4040":
            assert pytest.approx(water_perm, rel=1e-3) == 1.217e-11
            assert pytest.approx(salt_perm, rel=1e-3) == 4.726e-08
            assert pytest.approx(porosity, rel=1e-3) == 0.708
        if d["membrane_name"] == "AG8040F-400":
            assert pytest.approx(water_perm, rel=1e-3) == 1.070e-11
            assert pytest.approx(salt_perm, rel=1e-3) == 4.321e-08
            assert pytest.approx(porosity, rel=1e-3) == 0.704
        if d["membrane_name"] == "TMG20D-400":
            assert pytest.approx(water_perm, rel=1e-3) == 2.094e-11
            assert pytest.approx(salt_perm, rel=1e-3) == 2.802e-08
            assert pytest.approx(porosity, rel=1e-3) == 0.721


@pytest.mark.component
def test_iterate_estimation():
    for d in [fit.bw30_4040, fit.ag8040f_400, fit.tmg20d_400]:

        water_perm, salt_perm, porosity = fit.iterate_estimation(
            spec_sheet_params=d, num_iter=3, close_figs=True
        )

        if d["membrane_name"] == "BW30-4040":
            assert pytest.approx(water_perm, rel=1e-3) == 1.229e-11
            assert pytest.approx(salt_perm, rel=1e-3) == 5.614e-08
            assert pytest.approx(porosity, rel=1e-3) == 0.708
        if d["membrane_name"] == "AG8040F-400":
            assert pytest.approx(water_perm, rel=1e-3) == 1.084e-11
            assert pytest.approx(salt_perm, rel=1e-3) == 5.061e-08
            assert pytest.approx(porosity, rel=1e-3) == 0.704
        if d["membrane_name"] == "TMG20D-400":
            assert pytest.approx(water_perm, rel=1e-3) == 2.134e-11
            assert pytest.approx(salt_perm, rel=1e-3) == 3.289e-08
            assert pytest.approx(porosity, rel=1e-3) == 0.721
