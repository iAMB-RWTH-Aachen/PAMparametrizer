import matplotlib.pyplot as plt
from PAModelpy.utils import set_up_pam
from Scripts.i2_parametrization.pam_parametrizer_iML1515 import set_up_pamparametrizer
from Scripts.i3_analysis.PAMparametrizer_progress_cleaned_figure import plot_valid_data, plot_simulation

import pandas as pd


def set_up_ecoli_pam_parametrizer_and_get_substrate_uptake_rates():
    parametrizer = set_up_pamparametrizer(-12, -0.1, kcat_increase_factor=3)
    parametrizer._init_results_objects()
    substrate_rates = parametrizer._init_validation_df([parametrizer.min_substrate_uptake_rate,
                                                        parametrizer.max_substrate_uptake_rate])['EX_glc__D_e']
    substrate_rates = sorted(substrate_rates)
    return parametrizer, substrate_rates

sut, substrate_rates = set_up_ecoli_pam_parametrizer_and_get_substrate_uptake_rates()

substrate_rates = sorted(substrate_rates)
fig, axs = plot_valid_data(sut, fontsize=16)
fluxes, _ = sut.run_simulations_to_plot(substrate_uptake_id='EX_glc__D_e',
                                                 substrate_rates=substrate_rates,
                                                 sensitivity=False)
fig, axs = plot_simulation(fig, axs, fluxes, [abs(rate) for rate in substrate_rates]+[abs(substrate_rates[-1])],
                           sut.validation_data.get_by_id('EX_glc__D_e')._reactions_to_plot,
                           iteration=1, max_iteration=2, label='original')

data = {
    "id": [
        'P0A6E6_P0AB98_P0ABA0_P0ABA4_P0ABA6_P0ABB0_P0ABB4_P0ABC0_P68699',
        "P0A6P9", "P0A6P9", "P0A825", "P0A9B2", "P0A9B2",
        "P0A9P0_P0AFG3_P0AFG6", "P0A9P0_P0AFG3_P0AFG6",
        "P0AB80", "P0AB89", "P0AB89",
        "P0ABI8_P0ABJ1_P0ABJ3_P0ABJ6",
        "P36683", "P36683",
        "P77399", "P77399"
    ],
    "direction": ["f", "f", "b", "f", "f", "b", "f", "b", "f", "f", "b", "f", "f", "b", "f", "b"],
    "rxn_id": [
        "CE_ATPS4rpp_P0A6E6_P0AB98_P0ABA0_P0ABA4_P0ABA6_P0ABB0",
        "CE_ENO_P0A6P9", "CE_ENO_P0A6P9", "CE_THFAT_P0A825",
        "CE_GAPD_P0A9B2", "CE_GAPD_P0A9B2",
        "CE_AKGDH_P0A9P0_P0AFG3_P0AFG6", "CE_AKGDH_P0A9P0_P0AFG3_P0AFG6",
        "TYRTA", "CE_ADSL1r_P0AB89", "CE_ADSL1r_P0AB89",
        "CE_CYTBO3_4pp_P0ABI8_P0ABJ1_P0ABJ3_P0ABJ6",
        "CE_ACONTb_P36683", "CE_ACONTb_P36683",
        "CE_ECOAH1_P77399", "CE_ECOAH1_P77399"
    ],
    "type": ["mutation"] * 16,
    "value": [
        0.118800, 0.025059, 0.197559, 0.453631,
        0.012494, 0.012358, 0.002590, 0.111740,
        0.004976, 0.001467, 0.138531, 0.083330,
        0.025669, 0.005615, 0.001484, 0.009297
    ]
}

# Create the DataFrame
df = pd.DataFrame(data)

for i, row in df.iterrows():
    # current value is coefficient calculated as follows:
    # coeff = 1 / (kcat * 3600 * 1e-6) #3600 to convert to /h to /s *1e-6 to make calculations more accurate
    # need to convert the coefficient to the actual kcat
    new_kcat = row["value"]/ 3600 * 1e6
    direction = row["direction"]

    for pamodel in [sut._pamodel, sut.pamodel_no_sensitivity]:
        print(1/row['value'], new_kcat)
        try:
            print(pamodel.constraints[f"EC_{row['id']}_{direction}"])
            print(pamodel.enzymes.get_by_id(row['id']).rxn2kcat[row['rxn_id']])
        except:
            print('couldnt find enzyme', row)
            print(pamodel.enzymes.query(row['id'].split('_')[0]))
        pamodel.change_kcat_value(enzyme_id=row['id'], kcats={row["rxn_id"]:
                                                                           {direction: new_kcat}})

fluxes, _ = sut.run_simulations_to_plot(substrate_uptake_id='EX_glc__D_e',
                                                 substrate_rates=substrate_rates,
                                                 sensitivity=False)
fig, axs = plot_simulation(fig, axs, fluxes, [abs(rate) for rate in substrate_rates]+[abs(substrate_rates[-1])],
                           sut.validation_data.get_by_id('EX_glc__D_e')._reactions_to_plot,
                           iteration=2, max_iteration=2, label='modified')

plt.savefig('Results/test_mutation2.png')
