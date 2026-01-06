import pyomo.environ as pyo
from watertap.costing.util import (
    register_costing_parameter_block,
    make_capital_cost_var,
    make_fixed_operating_cost_var,
)


def build_default_chem_cost_param_block(blk):

    blk.cost = pyo.Var(
        initialize=1,
        doc="Default chem cost",
        units=pyo.units.USD_2020 / pyo.units.kg,
    )
    blk.purity = pyo.Var(
        initialize=1,
        doc="Default chem purity",
        units=pyo.units.dimensionless,
    )
    blk.capital_A_parameter = pyo.Var(
        initialize=15408,
        doc="Default chem addition capital cost A parameter",
        units=pyo.units.USD_2007,
    )
    blk.capital_b_parameter = pyo.Var(
        initialize=0.5479,
        doc="Default chem addition capital cost b parameter",
        units=pyo.units.dimensionless,
    )

    costing = blk.parent_block()
    costing.register_flow_type("default", blk.cost / blk.purity)


def build_alum_cost_param_block(blk):
    # CatCost v 1.1.1
    # Aluminum sulphate, 5-lb. bgs., c.l., works, frt. equald., 17% Al203, W. Coast
    blk.cost = pyo.Var(
        initialize=0.54,
        doc="Alum cost",
        units=pyo.units.USD_2020 / pyo.units.kg,
    )
    blk.purity = pyo.Var(
        initialize=1,
        doc="Alum purity",
        units=pyo.units.dimensionless,
    )
    blk.capital_A_parameter = pyo.Var(
        initialize=15408,
        doc="Alum addition capital cost A parameter",
        units=pyo.units.USD_2007,
    )
    blk.capital_b_parameter = pyo.Var(
        initialize=0.5479,
        doc="Alum addition capital cost b parameter",
        units=pyo.units.dimensionless,
    )

    costing = blk.parent_block()
    costing.register_flow_type("alum", blk.cost / blk.purity)


def build_ammonium_sulfate_cost_param_block(blk):
    blk.cost = pyo.Var(
        initialize=1.02,
        doc="Ammonium sulfate cost",
        units=pyo.units.USD_2021 / pyo.units.gallon,
    )
    blk.purity = pyo.Var(
        initialize=1,
        doc="Ammonium sulfate purity",
        units=pyo.units.dimensionless,
    )
    blk.capital_A_parameter = pyo.Var(
        initialize=6699.1,
        doc="Ammonium sulfate addition capital cost A parameter",
        units=pyo.units.USD_2007,
    )
    blk.capital_b_parameter = pyo.Var(
        initialize=0.4219,
        doc="Ammonium sulfate addition capital cost b parameter",
        units=pyo.units.dimensionless,
    )

    costing = blk.parent_block()
    costing.register_flow_type("ammonium_sulfate", blk.cost / blk.purity)


def build_ammonia_cost_param_block(blk):
    # CatCost v 1.1.1
    # Ammonia, US Gulf, spot c.f.r. Tampa
    blk.cost = pyo.Var(
        initialize=0.76,
        doc="Ammonia cost",
        units=pyo.units.USD_2020 / pyo.units.kg,
    )
    blk.purity = pyo.Var(
        initialize=1,
        doc="Ammonia purity",
        units=pyo.units.dimensionless,
    )
    blk.capital_A_parameter = pyo.Var(
        initialize=6699.1,
        doc="Ammonia addition capital cost A parameter",
        units=pyo.units.USD_2007,
    )
    blk.capital_b_parameter = pyo.Var(
        initialize=0.4219,
        doc="Ammonia addition capital cost b parameter",
        units=pyo.units.dimensionless,
    )

    costing = blk.parent_block()
    costing.register_flow_type("ammonia", blk.cost / blk.purity)


def build_calcium_hydroxide_param_block(blk):

    blk.cost = pyo.Var(
        initialize=2.3,
        doc="Calcium hydroxide cost",
        units=pyo.units.USD_2021 / pyo.units.kg,
    )
    blk.purity = pyo.Var(
        initialize=1,  # assumed
        doc="Calcium hydroxide purity",
        units=pyo.units.dimensionless,
    )
    blk.capital_A_parameter = pyo.Var(
        initialize=2262.8,
        doc="Calcium hydroxide addition capital cost A parameter",
        units=pyo.units.USD_2007,
    )
    blk.capital_b_parameter = pyo.Var(
        initialize=0.6195,
        doc="Calcium hydroxide addition capital cost b parameter",
        units=pyo.units.dimensionless,
    )

    costing = blk.parent_block()
    costing.register_flow_type("calcium_hydroxide", blk.cost / blk.purity)


def build_sodium_hydroxide_cost_param_block(blk):
    # CatCost v 1.1.1
    # Caustic soda (sodium hydroxide), liq., dst contract f.o.b.
    blk.cost = pyo.Var(
        initialize=2.37,
        doc="Sodium hydroxide cost",
        units=pyo.units.USD_2021 / pyo.units.kg,
    )
    blk.purity = pyo.Var(
        initialize=1,  # assumed
        doc="Sodium hydroxide purity",
        units=pyo.units.dimensionless,
    )
    blk.capital_A_parameter = pyo.Var(
        initialize=2262.8,
        doc="Sodium hydroxide addition capital cost A parameter",
        units=pyo.units.USD_2007,
    )
    blk.capital_b_parameter = pyo.Var(
        initialize=0.6195,
        doc="Sodium hydroxide addition capital cost b parameter",
        units=pyo.units.dimensionless,
    )

    costing = blk.parent_block()
    costing.register_flow_type("sodium_hydroxide", blk.cost / blk.purity)


def build_ferric_chloride_cost_param_block(blk):
    # CatCost v 1.1.1
    # Ferric chloride, technical grade, 100% basis, tanks, f.o.b. works
    blk.cost = pyo.Var(
        initialize=1053.68,
        doc="Ferric chloride cost",
        units=pyo.units.USD_2020 / pyo.units.kg,
    )
    blk.purity = pyo.Var(
        initialize=1,
        doc="Ferric chloride purity",
        units=pyo.units.dimensionless,
    )
    blk.capital_A_parameter = pyo.Var(
        initialize=34153,
        doc="Ferric chloride addition capital cost A parameter",
        units=pyo.units.USD_2007,
    )
    blk.capital_b_parameter = pyo.Var(
        initialize=0.319,
        doc="Ferric chloride addition capital cost b parameter",
        units=pyo.units.dimensionless,
    )

    costing = blk.parent_block()
    costing.register_flow_type("ferric_chloride", blk.cost / blk.purity)


def build_hydrochloric_acid_cost_param_block(blk):
    # CatCost v 1.1.1
    # Hydrochloric acid, 22 deg. Be, US Gulf dom. ex-works US NE
    blk.cost = pyo.Var(
        initialize=0.12,
        doc="HCl cost",
        units=pyo.units.USD_2020 / pyo.units.kg,
    )
    blk.purity = pyo.Var(
        initialize=0.36,  # 22 deg. Be
        doc="HCl purity",
        units=pyo.units.dimensionless,
    )
    blk.capital_A_parameter = pyo.Var(
        initialize=900.97,
        doc="Hydrochloric acid addition capital cost A parameter",
        units=pyo.units.USD_2007,
    )
    blk.capital_b_parameter = pyo.Var(
        initialize=0.6179,
        doc="Hydrochloric acid addition capital cost b parameter",
        units=pyo.units.dimensionless,
    )

    costing = blk.parent_block()
    costing.register_flow_type("hydrochloric_acid", blk.cost / blk.purity)


def build_lime_cost_param_block(blk):
    # CatCost v 1.1.1
    # Lime, hydrated, bulk, t.l., f.o.b. works
    blk.cost = pyo.Var(
        initialize=0.11,
        doc="Lime cost",
        units=pyo.units.USD_2020 / pyo.units.kg,
    )
    blk.purity = pyo.Var(
        initialize=1,
        doc="Lime purity",
        units=pyo.units.dimensionless,
    )
    blk.capital_A_parameter = pyo.Var(
        initialize=12985,
        doc="Lime addition capital cost A parameter",
        units=pyo.units.USD_2007,
    )
    blk.capital_b_parameter = pyo.Var(
        initialize=0.5901,
        doc="Lime addition capital cost b parameter",
        units=pyo.units.dimensionless,
    )

    costing = blk.parent_block()
    costing.register_flow_type("lime", blk.cost / blk.purity)


# def build_polymer_cost_param_block(blk):

#     blk.cost = pyo.Var(
#         initialize=0.17,
#         doc="Polymer cost",
#         units=pyo.units.USD_2020 / pyo.units.kg,
#     )
#     blk.purity = pyo.Var(
#         initialize=1,
#         doc="Polymer purity",
#         units=pyo.units.dimensionless,
#     )
#     blk.capital_A_parameter = pyo.Var(
#         initialize=15408,
#         doc="Polymer addition capital cost A parameter",
#         units=pyo.units.USD_2007,
#     )
#     blk.capital_b_parameter = pyo.Var(
#         initialize=0.5479,
#         doc="Polymer addition capital cost b parameter",
#         units=pyo.units.dimensionless,
#     )

#     costing = blk.parent_block()
#     costing.register_flow_type("polymer", blk.cost / blk.purity)


def build_soda_ash_cost_param_block(blk):
    # CatCost v 1.1.1
    # Soda ash (sodium carbonate), dense, US Gulf, f.o.b. bulk
    blk.cost = pyo.Var(
        initialize=0.23,
        doc="Soda ash cost",
        units=pyo.units.USD_2020 / pyo.units.kg,
    )
    blk.purity = pyo.Var(
        initialize=1,
        doc="Soda ash purity",
        units=pyo.units.dimensionless,
    )
    # TODO: Check these values; adopted from ferric chloride
    blk.capital_A_parameter = pyo.Var(
        initialize=34153,
        doc="Soda ash addition capital cost A parameter",
        units=pyo.units.USD_2007,
    )
    blk.capital_b_parameter = pyo.Var(
        initialize=0.319,
        doc="Soda ash addition capital cost b parameter",
        units=pyo.units.dimensionless,
    )

    costing = blk.parent_block()
    costing.register_flow_type("soda_ash", blk.cost / blk.purity)


def build_sodium_bisulfite_cost_param_block(blk):

    blk.cost = pyo.Var(
        initialize=1.405,
        doc="Sodium bisulfite cost",
        units=pyo.units.USD_2021 / pyo.units.gallon,
    )
    blk.purity = pyo.Var(
        initialize=1,
        doc="Sodium bisulfite purity",
        units=pyo.units.dimensionless,
    )
    blk.capital_A_parameter = pyo.Var(
        initialize=900.97,
        doc="Sodium bisulfite addition capital cost A parameter",
        units=pyo.units.USD_2007,
    )
    blk.capital_b_parameter = pyo.Var(
        initialize=0.6179,
        doc="Sodium bisulfite addition capital cost b parameter",
        units=pyo.units.dimensionless,
    )

    costing = blk.parent_block()
    costing.register_flow_type("sodium_bisulfite", blk.cost / blk.purity)


def build_sodium_hypochlorite_cost_param_block(blk):

    blk.cost = pyo.Var(
        initialize=1.29,
        doc="Sodium hypochlorite cost",
        units=pyo.units.USD_2021 / pyo.units.gallon,
    )
    blk.purity = pyo.Var(
        initialize=1,
        doc="Sodium hypochlorite purity",
        units=pyo.units.dimensionless,
    )
    blk.capital_A_parameter = pyo.Var(
        initialize=900.97,
        doc="Hypochlorite addition capital cost A parameter",
        units=pyo.units.USD_2007,
    )
    blk.capital_b_parameter = pyo.Var(
        initialize=0.6179,
        doc="Hypochlorite addition capital cost b parameter",
        units=pyo.units.dimensionless,
    )

    costing = blk.parent_block()
    costing.register_flow_type("sodium_hypochlorite", blk.cost / blk.purity)


def build_sulfuric_acid_cost_param_block(blk):

    blk.cost = pyo.Var(
        initialize=0.17,
        doc="Sulfuric acid cost",
        units=pyo.units.USD_2020 / pyo.units.kg,
    )
    blk.purity = pyo.Var(
        initialize=1,
        doc="Sulfuric acid purity",
        units=pyo.units.dimensionless,
    )
    blk.capital_A_parameter = pyo.Var(
        initialize=900.97,
        doc="Sulfuric acid addition capital cost A parameter",
        units=pyo.units.USD_2007,
    )
    blk.capital_b_parameter = pyo.Var(
        initialize=0.6179,
        doc="Sulfuric acid addition capital cost b parameter",
        units=pyo.units.dimensionless,
    )

    costing = blk.parent_block()
    costing.register_flow_type("sulfuric_acid", blk.cost / blk.purity)


def build_scale_inhibitor_cost_param_block(blk):

    blk.cost = pyo.Var(
        initialize=2,
        doc="Scale inhibitor cost",
        units=pyo.units.USD_2020 / pyo.units.gallon,
    )
    blk.purity = pyo.Var(
        initialize=1,
        doc="Scale inhibitor purity",
        units=pyo.units.dimensionless,
    )
    blk.capital_A_parameter = pyo.Var(
        initialize=900.97,
        doc="Scale inhibitor addition capital cost A parameter",
        units=pyo.units.USD_2007,
    )
    blk.capital_b_parameter = pyo.Var(
        initialize=0.6179,
        doc="Scale inhibitor addition capital cost b parameter",
        units=pyo.units.dimensionless,
    )

    costing = blk.parent_block()
    costing.register_flow_type("scale_inhibitor", blk.cost / blk.purity)


def cost_chemical_addition(blk, cost_capital=False):

    chem_build_rule_dict = {
        "default": build_default_chem_cost_param_block,
        "ammonia": build_ammonia_cost_param_block,
        "ammonium_sulfate": build_ammonium_sulfate_cost_param_block,
        "lime": build_lime_cost_param_block,
        "ferric_chloride": build_ferric_chloride_cost_param_block,
        "soda_ash": build_soda_ash_cost_param_block,
        "alum": build_alum_cost_param_block,
        # "polymer": build_polymer_cost_param_block,
        "calcium_hydroxide": build_calcium_hydroxide_param_block,
        "sodium_hydroxide": build_sodium_hydroxide_cost_param_block,
        "sulfuric_acid": build_sulfuric_acid_cost_param_block,
        "scale_inhibitor": build_scale_inhibitor_cost_param_block,
        "sodium_hypochlorite": build_sodium_hypochlorite_cost_param_block,
        "sodium_bisulfite": build_sodium_bisulfite_cost_param_block,
        "hydrochloric_acid": build_hydrochloric_acid_cost_param_block,
    }

    chemical = blk.unit_model.config.chemical
    chem_build_rule = chem_build_rule_dict.get(chemical, None)
    if chem_build_rule is None:
        raise ValueError(f"Unrecognized chemical type {chemical} in ChemAddition")
    

    @register_costing_parameter_block(
        build_rule=chem_build_rule, parameter_block_name=chemical
    )
    def cost_chem_addition(blk):
        chem_addition_param_blk = blk.costing_package.find_component(f"{chemical}")
        if cost_capital:
            make_capital_cost_var(blk)
            blk.costing_package.add_cost_factor(blk, "TPEC")
            chem_flow_mass_dim = pyo.units.convert(
                blk.unit_model.chemical_flow_mass / (pyo.units.lb / pyo.units.day),
                to_units=pyo.units.dimensionless,
            )

            blk.capital_cost_constraint = pyo.Constraint(
                expr=blk.capital_cost
                == blk.cost_factor
                * pyo.units.convert(
                    chem_addition_param_blk.capital_A_parameter
                    * chem_flow_mass_dim**chem_addition_param_blk.capital_b_parameter,
                    to_units=blk.costing_package.base_currency,
                )
            )
            
        cost_units = pyo.units.get_units(chem_addition_param_blk.cost)

        if any(_x in cost_units.to_string() for _x in ["gal", "gallon", "m**3"]):
            blk.costing_package.cost_flow(
                blk.unit_model.chemical_soln_flow_vol, chemical
            )
        else:
            blk.costing_package.cost_flow(
                blk.unit_model.chemical_flow_mass, chemical
            )

    cost_chem_addition(blk)
