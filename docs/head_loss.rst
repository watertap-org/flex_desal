Head Loss Model
===============

Overview
--------
The head loss model (`HeadLossData`) represents pressure drops in system components such as pipes, fittings, or other unit operations. It is a simple unit model that enforces a specified pressure drop between inlet and outlet streams, while maintaining material and (optionally) energy balances.

Key Features
------------
- Models pressure drop (head loss) across a unit.
- Supports flexible configuration of material, energy, and momentum balances.
- Assumes isothermal operation by default.
- Integrates with WaterTAP and IDAES property packages.

Configuration Options
---------------------
The model is configured via the following options:

- ``dynamic``: Must be False (steady-state only).
- ``has_holdup``: Must be False (no holdup volume).
- ``property_package``: Property package for state variables.
- ``property_package_args``: Arguments for property package construction.
- ``material_balance_type``: Type of material balance (default: use property package default).
- ``is_isothermal``: If True, enforces isothermal conditions (default: True).
- ``energy_balance_type``: Type of energy balance (default: none).
- ``momentum_balance_type``: Type of momentum balance (default: pressureTotal).

Variables and Parameters
-----------------------
- ``deltaP``: Pressure drop across the unit (Pa).
- Inlet and outlet state variables (flow, pressure, temperature, composition) as defined by the property package.

Ports
-----
- ``inlet``: Inlet port for the unit.
- ``outlet``: Outlet port for the unit.

Constraints
-----------
- Material balance (configurable type).
- Energy balance (optional, default: none).
- Momentum balance (default: pressure drop).
- Isothermal constraint (if enabled).

Example Usage
-------------
.. code-block:: python

	from flex_desal.src.models.head_loss import HeadLoss
	# Configure and add to flowsheet
	head_loss = HeadLoss(property_package=..., ...)
	# Fix deltaP to specify the pressure drop
	head_loss.deltaP.fix(1e5)  # Pressure drop of 100 kPa

See WaterTAP documentation for further details on unit model integration and flowsheet setup.
