from pathlib import Path
import pandas as pd


def _load_results_dataframe(filename, data_type="optimization_results"):
    script_dir = Path(__file__).resolve().parent
    file_path = Path(data_type + "/" + filename)
    if not file_path.is_absolute():
        file_path = script_dir / file_path

    return pd.read_csv(file_path), file_path


def _calc_flex_degree(df):
    num_days = (
        167 / 24
    )  # THIS IS ONLY TRUE FOR RESULTS OVER A WEEK. IDK why this number is off by an hour, but just roll with it.
    column_names = [
        "reverse_osmosis.ro_skid[1].shutdown",
        "reverse_osmosis.ro_skid[2].shutdown",
        "reverse_osmosis.ro_skid[3].shutdown",
        "reverse_osmosis.ro_skid[4].shutdown",
    ]
    num_shutdowns = df[column_names].sum().sum()
    return min(1, num_shutdowns / (4 * num_days)), num_days


def calc_new_replacement_cost(
    filename, max_degradation=0.1, data_type="optimization_results"
):
    # Example calculation for new replacement cost based on given parameters
    membrane_total_cost = 500 * 4 * (72 + 30 + 15)
    membrane_lifetime = 5  # years
    motor_total_cost = (
        125000 * 4
    )  # This is only considering the RO motors... and not the UF motors.
    motor_lifetime = 17.5  # years

    df, _ = _load_results_dataframe(filename, data_type=data_type)
    flex_degree, num_days = _calc_flex_degree(df)
    num_months = num_days / 31

    memb_rep_cost = (
        membrane_total_cost
        * (num_months / 12)
        / (membrane_lifetime * (1 - max_degradation * flex_degree))
    )
    motor_rep_cost = (
        motor_total_cost
        * (num_months / 12)
        / (motor_lifetime * (1 - max_degradation * flex_degree))
    )

    total_replacement_cost = memb_rep_cost + motor_rep_cost

    return total_replacement_cost


def replacement_cost_residual(max_degradation, filename, target_cost=4000.0):
    return (
        calc_new_replacement_cost(filename, max_degradation=max_degradation)
        - target_cost
    )


def solve_max_degradation_for_target(
    filename,
    data_type="optimization_results",
    target_cost=4000.0,
    lower=0.0,
    upper=0.999999,
    tol=1e-8,
    max_iter=200,
):
    df, _ = _load_results_dataframe(filename, data_type=data_type)
    flex_degree, _ = _calc_flex_degree(df)

    if flex_degree <= 0:
        flat_cost = calc_new_replacement_cost(filename, max_degradation=0.0)
        if abs(flat_cost - target_cost) < tol:
            return lower
        raise ValueError(
            "No shutdown flexibility detected (flex_degree=0), so replacement cost "
            "does not vary with max_degradation. "
            f"Current cost is ${flat_cost:,.2f}, target is ${target_cost:,.2f}."
        )

    max_feasible_deg = (1.0 / flex_degree) * (1 - 1e-9)
    effective_upper = min(upper, max_feasible_deg)

    f_low = replacement_cost_residual(lower, filename, target_cost)
    f_high = replacement_cost_residual(effective_upper, filename, target_cost)

    # If needed, expand bracketing up to the asymptote where cost diverges.
    if f_low * f_high > 0 and effective_upper < max_feasible_deg:
        effective_upper = max_feasible_deg
        f_high = replacement_cost_residual(effective_upper, filename, target_cost)

    if f_low == 0:
        return lower
    if f_high == 0:
        return effective_upper

    if f_low * f_high > 0:
        raise ValueError(
            "Could not bracket a root for max_degradation in "
            f"[{lower}, {effective_upper}]. Residuals are f(lower)={f_low:.4f}, "
            f"f(upper)={f_high:.4f}."
        )

    lo = lower
    hi = effective_upper
    for _ in range(max_iter):
        mid = 0.5 * (lo + hi)
        f_mid = replacement_cost_residual(mid, filename, target_cost)
        if abs(f_mid) < tol or abs(hi - lo) < tol:
            return mid
        if f_low * f_mid < 0:
            hi = mid
            f_high = f_mid
        else:
            lo = mid
            f_low = f_mid

    raise RuntimeError(
        f"Bisection did not converge in {max_iter} iterations for target_cost={target_cost}."
    )


# Plant Data
filename = "wrd_result_summer_DR_1_hr_startup.csv"

rep_cost = calc_new_replacement_cost(
    filename, data_type="optimization_results", max_degradation=0.4
)
print(f"Calculated replacement cost based on shutdowns: ${rep_cost:,.2f}")

target_cost = 3718 + 1410  # 1410 is the baseline replacement cost (full lifetime)
max_deg_solution = solve_max_degradation_for_target(filename, target_cost=target_cost)
check_cost = calc_new_replacement_cost(filename, max_degradation=max_deg_solution)
print(
    f"max_degradation for target replacement cost ${target_cost:,.2f}: "
    f"{max_deg_solution:.6f}"
)
print(f"Check replacement cost at this value: ${check_cost:,.2f}")
