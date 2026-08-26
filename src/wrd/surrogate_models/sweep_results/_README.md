FILE DESCRIPTIONS:

RESULTS: Data used to create the surrogates used in Pricetaker (PT)
PT_ro_train_sweep -> The PYSMO model is created in energy_surrogate_ro_train.py. Will be output by param_sweep script in its current state.

PT_uf_pump_sweep -> There is no PYSMO surrogate required because enegry is treated as a function of one variable (flow).

TESTS: Other sweeps used to validate some assumptions

S3_70_eff_sweep -> third stage pump efficiency of 70% was used. Unclear why the outputs of this file are different from PT_ro_train_sweep.

S3_40_RR_lim -> changed the recovery limit for the third stage to 40%.

S1_S2_equal_RR -> Additional constraint that stage one and two must have the same recovery. This is close to how the plant typically operates, based on data.

S1_S2_not_equal_RR -> Result with the above constraint removed

S1_S2_equal_S3_40_RR_lim -> combination of including a recovery limit for the third stage and recovery equality between the first two stages

NOTE:
- In some of the TEST files above, the energy consumption of the third stage is negative because a constraint on work > 0 had not yet been included.