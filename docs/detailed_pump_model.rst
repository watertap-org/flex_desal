Detailed Pump Model
===================

Overview
--------
The detailed pump model (`PumpIsothermalData`) provides an advanced isothermal pump unit for process modeling. It extends the standard pump model by incorporating pump performance curves, variable efficiency options, and system curve calculations for more accurate simulation of pump behavior under varying conditions.

Key Features
------------
- Supports both fixed and variable pump efficiency (based on flow and head).
- Allows pump curve definition via dataset (CSV) or surrogate coefficients.
- Calculates system curve constants and reference points for pump operation.
- Integrates motor and variable frequency drive (VFD) efficiency.
- Provides constraints for isothermal operation and pump performance.

Configuration Options
---------------------
The model is configured via the following options:

- ``variable_efficiency``: Selects between fixed efficiency (default) or variable efficiency based on flow.
- ``pump_curve_data_type``: Specifies whether pump curve data is provided as a dataset or as surrogate coefficients.
- ``head_surrogate_coeffs``: Surrogate coefficients for the pump head curve (cubic polynomial).
- ``efficiency_surrogate_coeffs``: Surrogate coefficients for the pump efficiency curve (cubic polynomial).
- ``pump_curves``: DataFrame or CSV file containing flow, head, and efficiency data for curve fitting.

Pump Curve Data
---------------
Pump curves can be provided in two ways:

1. **Dataset (CSV):**
	- CSV file must contain columns: ``flow (m3/s)``, ``head (m)``, ``efficiency (-)``.
	- The model fits cubic polynomials to the data for head and efficiency.
2. **Surrogate Coefficients:**
	- User provides cubic polynomial coefficients for head and efficiency as dictionaries.
	- Example: ``{0: a, 1: b, 2: c, 3: d}`` for ``head = a + b*flow + c*flow^2 + d*flow^3``.

Variables and Parameters
-----------------------
- ``design_flow``: Design flowrate (m³/s).
- ``design_head``: Design head (m).
- ``design_efficiency``: Design efficiency (dimensionless).
- ``design_speed_fraction``: Design speed fraction (dimensionless).
- ``system_curve_geometric_head``: Static head constant (m).
- ``system_curve_flow_constant``: Flow constant for losses (m·(m³/s)⁻²).
- ``ref_flow``: Reference flow from pump datasheet (m³/s).
- ``ref_head``: Reference head from pump datasheet (m).
- ``ref_speed_fraction``: Reference speed fraction from pump datasheet (dimensionless).
- ``vfd_efficiency``: VFD efficiency (default 0.97).
- ``motor_efficiency``: Motor efficiency (default 0.95).


System Curve Equation
--------------------
The system curve describes the relationship between the required pump head and the flow rate in the system. It accounts for both the static (geometric) head and the dynamic head losses due to friction and fittings. The general form of the system curve equation is:

.. math::
	H_{system}(Q) = H_{static} + K Q^2

where:

	- :math:`H_{system}(Q)`: Total head required by the system at flow rate :math:`Q` (m)
	- :math:`H_{static}`: Static (geometric) head, representing elevation difference or constant pressure requirement (m)
	- :math:`K`: System loss coefficient, representing friction and minor losses (m·(m³/s)⁻²)
	- :math:`Q`: Volumetric flow rate (m³/s)

In the model, these are represented by the variables ``system_curve_geometric_head`` (for :math:`H_{static}`) and ``system_curve_flow_constant`` (for :math:`K`). The system curve is used to determine the required pump head for a given flow rate and is essential for matching the pump performance to the process requirements.

Constraints
-----------
- Isothermal balance: Ensures inlet and outlet temperatures are equal.
- Design head and flow constraints: Link pump operation to system curve.
- System curve calculation: Defines geometric and flow constants.
- Surrogate constraints: Calculate head and efficiency using cubic polynomials.
- Reference point constraints: Solve system and pump curve simultaneously.
- Overall efficiency: Combines pump, motor, and VFD efficiency.


References
----------
- WaterTAP core and costing modules
- IDAES core and unit models


Variables to Fix for 0 Degrees of Freedom
-----------------------------------------
To achieve zero degrees of freedom (DOF) and enable the model to solve, the following variables must be fixed, as demonstrated in the test file:

**For variable efficiency (Efficiency.Flow):**
- Inlet stream variables:
	- ``inlet.flow_mass_phase_comp[0, "Liq", "TDS"]`` (when flow is specified)
	- ``inlet.flow_mass_phase_comp[0, "Liq", "H2O"]`` (when flow is specified)
	- ``inlet.pressure[0]``
	- ``inlet.temperature[0]``
- Outlet pressure:
	- ``outlet.pressure[0]`` (when head is specified)
- System curve and pump parameters:
	- ``system_curve_geometric_head``
	- ``ref_speed_fraction``
	- Optionally, ``design_speed_fraction`` (when speed is specified)

**For fixed efficiency (Efficiency.Fixed):**
- Inlet stream variables:
	- ``inlet.flow_mass_phase_comp[0, "Liq", "TDS"]``
	- ``inlet.flow_mass_phase_comp[0, "Liq", "H2O"]``
	- ``inlet.pressure[0]``
	- ``inlet.temperature[0]``
- Pump parameters:
	- ``efficiency_pump``
	- ``deltaP``

These settings ensure the model is fully specified and ready for initialization and solving. See the test file for example code.
