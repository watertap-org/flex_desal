Source Model
============

Overview
--------
The source model (`SourceData`) represents a feed or source stream entering the process. It is a simple unit model that provides specified inlet conditions (flow, composition, pressure, temperature) to the flowsheet. The model is based on the standard IDAES Feed unit and is extended for WaterTAP integration and costing.

Key Features
------------
- Provides a fixed inlet stream to the process.
- Supports specification of all inlet state variables.
- Integrates with WaterTAP and IDAES property packages.
- Includes a default costing method for source water.

Configuration Options
---------------------
The model is configured via the same options as the IDAES Feed unit:

- ``property_package``: Property package for state variables.
- ``property_package_args``: Arguments for property package construction.

Variables and Parameters
-----------------------
- Inlet state variables (flow, pressure, temperature, composition) as defined by the property package.

Ports
-----
- ``outlet``: Outlet port for the unit.

Example Usage
-------------
.. code-block:: python

	from flex_desal.src.models.source import Source
	# Configure and add to flowsheet
	source = Source(property_package=..., ...)
	# Fix inlet state variables as needed
	source.outlet.flow_vol.fix(1.0)
	source.outlet.pressure.fix(101325)
	source.outlet.temperature.fix(298.15)


Costing
-------

The model includes a default costing method (`cost_source`) that assigns a variable operating cost per unit volume of source water. The key features of the costing implementation are:

- **Variable Operating Cost:**
	- The variable operating cost is calculated as:
		``variable_operating_cost = unit_cost * base_period_flow``
		where ``unit_cost`` is the cost per cubic meter of source water (default: 0.15 in base currency), and ``base_period_flow`` is the volumetric flow of water over the base period.
- **Capital Cost:**
	- Capital cost is set to zero by default, as the source is assumed to have no capital investment in this model.
- **Costing Variables and Constraints:**
	- ``unit_cost``: Source cost per cubic meter (can be adjusted in the costing package).
	- ``variable_operating_cost``: Total variable operating cost for the base period.
	- ``capital_cost``: Set to zero.
	- Constraints are included to enforce these relationships.
- **Integration:**
	- The costing method is compatible with WaterTAP and IDAES costing frameworks, and can be extended or customized as needed.

To use costing, ensure the unit is attached to a costing package and the relevant variables are referenced in reports or objective functions as needed.

See WaterTAP documentation for further details on unit model integration and flowsheet setup.
