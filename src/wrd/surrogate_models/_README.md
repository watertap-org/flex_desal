FILE DESCRIPTIONS:

energy_surrogate_ro_train.py -> script that outputs the polynomial fit .json for the ro train energy consumption using PYSMO.
                                Figures for the surrogate model validation are included in this script

param_sweep_fxs_ro_train.py -> Functions required for the param sweep. Building the ro_train flowsheet, the required outputs, the initialization and optimization routines. 

param_sweep_fxs_uf_pump.py -> Functions required for the param sweep. Building the uf_pump flowsheet, the required outputs, the initialization and optimization routines. 

param_sweep_ro_train.py -> Implementation of the parameter sweep fuction (https://watertap.readthedocs.io/en/latest/how_to_guides/how_to_use_parameter_sweep.html#module-documentation). 
                         Determines the specific energy of an RO train based on recovery and flowrate. 
                         This includes three pumps, one per RO stage.

param_sweep_uf_pump.py -> Same parameter sweep function for one UF pump. 
                                 This is only a function of flowrate.
                                 It is assumed the pressure required by the RO is constant and that the pressure drop will be constant across the UF membranes.

ro_SEC_poly_fit_order_X.json -> PYSMO output file for one train with specific energy given as a function of the recovery and flowrate with polynomial order X.

NOTES:
- There is no energy_surrogate_uf_pump script. A second order polynomial fit between SEC and was performed in excel using data in sweep_results/PT_uf_pump_sweep.csv
- In addition to total energy, the ro_train sweep tracked the pressures, first elemement product flow, pump speed, and permeate concentration at each stage. These are recorded in sweep_results/PT_ro_train_sweep.csv

