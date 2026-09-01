FILE DESCRIPTIONS:

calc_membrane_properties --> script that reads input pressure, flowrate, and salinity data from each ro stage, then creates a 1-D ro unit model and initializes using sample A and B values. Then the permeate flowrate and salinity are fixed with A and B unfixed. The resulting values are recorded. The A and B values are tabulated for each operational data point and saved to output_data.

ro_for_membrane_properties --> Adapted version of ro.py in components that is compatable with the calc_membrane_properties script. Changes to how data is passed, not with the unit model itself.

NOTES:
- The tests folder is currently empty, but I was considered writing a pytest file. That would not be in the documentation PR though.
