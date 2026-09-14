FILE DESCRIPTIONS:
- plot_from_data => script that recreates the plot from the pricetaker formulation. This allows tweaking formatting of the figures without re-running the optimization. Also has option to take properly formatted plant data. 
- replacement_costs.py => script that calculates the replacement costs based on the degree of flexibility of an operation scheme.

- paper_fig => will ultimately include all the plots that area directly used in the paper. plot_from data saves plots to this folder
- optimization_results => csv files with hourly operational states. Output from pricetaker and the optimization. 
- validation_plot => see README in that folder.
- plant_data => redundant with files in the "Comparison to Real Week"

NOTES:
- Currently, an additional column has to be added to the csvs to represent the peak hours. This can be easily copied over from one of the existing csvs.