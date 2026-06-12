import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from typing import List
import os
import matplotlib.ticker as ticker
from matplotlib.lines import Line2D
import seaborn as sns


from cobra.io import read_sbml_model
from PAModelpy.utils import set_up_pam

from Modules.PAMparametrizer.utils.pam_generation import create_pamodel_from_diagnostics_file
from Modules.PAMparametrizer.utils.pamparametrizer_visualization import RXN_NAME_MAPPER
from Modules.PAMparametrizer.utils.pamparametrizer_analysis import calculate_error_for_reactions

PARAM_FILE_OLD = os.path.join('Results', '1_preprocessing','proteinAllocationModel_iML1515_EnzymaticData_250912.xlsx')
PARAM_FILE_PREPROC = os.path.join('Results','2_parametrization','proteinAllocationModel_iML1515_EnzymaticData_multi.xlsx')
PARAM_FILE_PREPROC_6 = os.path.join('Results','2_parametrization','proteinAllocationModel_iML1515_EnzymaticData_multi_6.xlsx')
PARAM_FILE_PREPROC_8 = os.path.join('Results','2_parametrization','proteinAllocationModel_iML1515_EnzymaticData_multi_8.xlsx')
RANDOM_PARAMETERS = [os.path.join(
    'Results','3_analysis', 'parameter_files',f'proteinAllocationModel_EnzymaticData_iML1515_random_{i}.xlsx')
    for i in range(1, 4)
]
EXCHANGE_FLUX_FILE_PATH = os.path.join('Data', 'Ecoli_phenotypes', 'Ecoli_phenotypes_py_rev.xls')
MFA_DATA_FILE_PATH = os.path.join('Data', 'Ecoli_phenotypes', 'fluxomics_datasets_gerosa.xlsx')

MODEL_FILE = os.path.join('Models', 'iML1515.xml')
SUBSTRATE_ID = 'EX_glc__D_e'
GEROSA_GLC_UPTAKE = -9.654

RESULT_FLUX_PATH = os.path.join('Results', '3_analysis', 'iML1515_alternative_models_predictions_diff_growth_rates.xlsx')

def save_flux_predictions(model_files: List[str],
                         substrate_uptake_rates: List[float],
                         diagnostic_files: List[str],
                         label_names:List[str],
                         simulated_fluxes_result_file: str = RESULT_FLUX_PATH,
                         ):
    pam = set_up_pam(pam_info_file = PARAM_FILE_PREPROC, model = MODEL_FILE, sensitivity = False)
    models = [set_up_pam(pam_info_file = file, model = MODEL_FILE, sensitivity = False)
              if 'xml' not in file else read_sbml_model(file)
              for file in model_files
    ]
    alt_pams = [create_pamodel_from_diagnostics_file(diag_file, pam.copy_with_pickle()) for diag_file in diagnostic_files]

    for label, model in zip(label_names, models+alt_pams):
        print('------------------------------------------------------------------------')
        print('Analyzing flux distribution of model', label, '\n')
        all_fluxes = pd.DataFrame()
        for substrate_uptake_rate in substrate_uptake_rates:
            model.reactions.get_by_id(SUBSTRATE_ID).lower_bound = substrate_uptake_rate
            model.optimize()
            if model.solver.status != 'optimal': continue

            all_reactions = [rxn for rxn in model.reactions if 'CE_' not in rxn.id]#CE to ignore the catalytic events
            fluxes = (pd.DataFrame.from_dict({rxn.id: rxn.flux for rxn in all_reactions},
                                  columns =['flux'], orient='index',)
                      .assign(substrate_uptake_rate = substrate_uptake_rate)
                      .reset_index(names=['rxn_id'])
                      )
            all_fluxes = pd.concat([all_fluxes, fluxes])

        kwargs = {'mode': 'a', 'if_sheet_exists': 'replace'} \
            if os.path.exists(simulated_fluxes_result_file) else {'mode': 'w'}
        with pd.ExcelWriter(simulated_fluxes_result_file,
                            engine='openpyxl', **kwargs) as writer:
            all_fluxes.to_excel(writer, sheet_name=label, index= False)

def perform_and_save_simulations(num_alternative_models: int,
                                 substrate_uptake_rates: List[float],
                                 simulated_fluxes_result_file: str = RESULT_FLUX_PATH,
                                 ):
    diagnostic_files = list()
    param_files = list()
    all_model_labels = ['iML1515', 'GotEnzymes', 'Scaled 6x','Scaled 7x','Scaled 8x',
                        'Randomized 1','Randomized 2','Randomized 3'
                        ] \
                       + [f'alternative {i}' for i in range(1, num_alternative_models + 1)]
    for file_nmbr in range(1, num_alternative_models + 1):
        diagnostic_files += [os.path.join('Results', '2_parametrization', 'diagnostics',
                                          f'pam_parametrizer_diagnostics_{file_nmbr}.xlsx')]
        param_files += [os.path.join('Results', '3_analysis', 'parameter_files',
                                     f'proteinAllocationModel_EnzymaticData_iML1515_{file_nmbr}.xlsx')]

    # save_flux_predictions(model_files=[MODEL_FILE, PARAM_FILE_OLD, PARAM_FILE_PREPROC_6,
    #                                   PARAM_FILE_PREPROC, PARAM_FILE_PREPROC_8
    #                                    ]+RANDOM_PARAMETERS,
    #                       substrate_uptake_rates=substrate_uptake_rates,
    #                       diagnostic_files=diagnostic_files,
    #                       label_names=all_model_labels,
    #                       simulated_fluxes_result_file=simulated_fluxes_result_file,
    #                       )
    save_flux_predictions(model_files=[PARAM_FILE_PREPROC_6,
                                       PARAM_FILE_PREPROC, PARAM_FILE_PREPROC_8
                                       ],
                          substrate_uptake_rates=substrate_uptake_rates,
                          diagnostic_files=[],
                          label_names=['Scaled 6x','Scaled 7x','Scaled 8x'],
                          simulated_fluxes_result_file=simulated_fluxes_result_file,
                          )


def plot_exchange_fluxes_vs_experiments(all_fluxes: pd.DataFrame,
                                        exchange_fluxes: pd.DataFrame,
                                        axs: List[plt.Axes] = None,
                                        reactions_to_plot: List[str]= ['BIOMASS_Ec_iML1515_core_75p37M',
                                                                       'EX_ac_e', 'EX_co2_e', 'EX_o2_e'],
                                        fontsize = 12,
                                        ):
    exchange_fluxes = exchange_fluxes.sort_values(by='EX_glc__D_e', ascending=False)
    all_fluxes = all_fluxes.sort_values(by='substrate_uptake_rate', ascending=False)
    all_alternative_ids = [f'alternative {i}' for i in range(1,11)]
    cmap = {
            **{'iML1515': 'black', 'GotEnzymes': 'black', 'Randomized': 'midnightblue'},
        **{f'Scaled {i}x': color for i, color in zip(list(range(6,9)), sns.color_palette('plasma', n_colors=3))},
        **{alt: 'lightblue' for alt in all_alternative_ids},
            **{f'Randomized {i+1}':'lightcoral' for i in range(3)},}
    linestyles = {**{alt_id:'-' for alt_id in all_alternative_ids},
                  **{f'Randomized {i}':":" for i in range(1,4)},
                  **{'iML1515': '--', 'GotEnzymes': '-.', 'Randomized': ':', 'After preprocessing':'-'},
                  **{f'Scaled {i}x':'--' for i in range(6,9)},}
    model_order = (['iML1515', 'GotEnzymes', 'Randomized', 'After preprocessing']
                   +[f'Scaled {i}x' for i in range(6,9)]
                   + all_alternative_ids
                   + [f'Randomized {i+1}' for i in range(3)])
    if axs is None:
        fig, axs = plt.subplots(2, 2, dpi=100, figsize = (21/2.54, 15/2.54))
        axs = axs.flatten()

        for r, ax in zip(reactions_to_plot, axs):
            if r not in RXN_NAME_MAPPER.keys(): continue
            # plot data
            x = [abs(glc) for glc in exchange_fluxes['EX_glc__D_e']]
            y = [abs(data) for data in exchange_fluxes[r]]
            ax.set_ylabel(RXN_NAME_MAPPER[r].replace('[mmol', '\n[mmol'), fontsize=fontsize, labelpad =5)
            ax.scatter(x, y,
                       color='black', marker='o', s=20,
                       facecolors=None, zorder=0,
                       label='Literature')

            reaction_fluxes = all_fluxes[all_fluxes['rxn_id'] == r]
            for model in model_order:
                flux_df = reaction_fluxes[reaction_fluxes['model'] == model]
                if model == 'Randomized' or model == 'After preprocessing': continue
                ax.plot(flux_df['substrate_uptake_rate'].abs(), flux_df['flux'].abs(),
                        color = cmap[model],
                        linestyle = linestyles[model],
                        alpha = 1 if not 'alternative' in model else 0.8,
                        label = model
                )

            # ax.set_xlabel(RXN_NAME_MAPPER['EX_glc__D_e'], fontsize=fontsize)
            # ax.set_ylim([-0.05,reaction_fluxes.flux.max()*1.2])
            ax.tick_params(axis='both', labelsize=fontsize)
            ax.yaxis.set_major_formatter(ticker.ScalarFormatter(useOffset=False, useMathText=False))
            ax.yaxis.get_major_formatter().set_scientific(False)  # Disable scientific notation
            ax.yaxis.get_major_formatter().set_useLocale(True)  # Ensure proper decimal formatting

        fig.text(0.5, 0.15, RXN_NAME_MAPPER['EX_glc__D_e'], ha='center', fontsize=fontsize)
        h, l = ax.get_legend_handles_labels()
        handles, labels = [],[]
        for handle, label in zip(h, l):
            if 'alternative' in label or 'Randomized' in label: continue
            handles.append(handle)
            labels.append(label)
        handles+= [Line2D([0], [0], color='lightblue', alpha = 0.8), Line2D([0], [0], color = 'indianred')]
        labels += ['Alternative PAMs', 'Randomized PAMs']

        fig.legend(handles = handles, labels=labels,
                   loc ='lower center', bbox_to_anchor=(0.5, 0), fontsize=11,ncols=3)

        fig.tight_layout()
        fig.subplots_adjust(bottom = 0.25, wspace=0.3)


        return fig, axs

def determine_coeff_of_corr_sim_vs_measurements(all_exp_data: pd.DataFrame,
                                                all_fluxes: pd.DataFrame,):
    errors_for_models_rows = []
    rxns_to_validate_intracell = [rxn for rxn in list(all_exp_data.columns)
                            if not any([substr in rxn for substr in ['Yield', '_ub', '+', 'EX_lac_e', 'EX_']])
                        ]
    rxns_to_validate_extracell = [rxn for rxn in list(all_exp_data.columns)
                                  if not any([substr in rxn for substr in ['Yield', '_ub', '+', 'EX_lac_e']])
                                  and 'EX_' in rxn
                                  ]
    rxns_to_validate_extracell = rxns_to_validate_intracell + rxns_to_validate_extracell

    for model, flux_df in all_fluxes.groupby('model'):
        flux_df_wide =(flux_df.drop_duplicates(['rxn_id', 'substrate_uptake_rate'])
                       .pivot(columns='rxn_id', values='flux', index = 'substrate_uptake_rate')
                       .reset_index()
                       .rename(columns={'substrate_uptake_rate':'substrate'})
                       )
        flux_df_wide['substrate'] = flux_df_wide['EX_glc__D_e']

        r_squared_intracell  = calculate_error_for_reactions(validation_df=all_exp_data,
                                                   flux_df=flux_df_wide.assign(substrate_id = 'EX_glc__D_e'),
                                                   rxns_to_validate=[rxn for rxn in rxns_to_validate_intracell if rxn in flux_df_wide.columns],)
        r_squared_extracell  = calculate_error_for_reactions(validation_df=all_exp_data,
                                                   flux_df=flux_df_wide.assign(substrate_id = 'EX_glc__D_e'),
                                                   rxns_to_validate=[rxn for rxn in rxns_to_validate_extracell if rxn in flux_df_wide.columns],)
        r_squared = r_squared_intracell + r_squared_extracell

        for i, rxn in enumerate([rxn for rxn in rxns_to_validate if rxn in flux_df_wide.columns]):
            errors_for_models_rows.append({'model': model,
                                           'rxn_id':rxn,
                                           'r_squared': r_squared[i],
                                           'mean_rsquared_intra': np.mean(r_squared_intracell),
                                           'mean_rsquared_extra': np.mean(r_squared_extracell),
                                           })
    return pd.DataFrame(errors_for_models_rows)

if __name__ == '__main__':
    mfa_data = pd.read_excel(MFA_DATA_FILE_PATH)
    mfa_data_glc = (mfa_data[mfa_data.condition == 'Glucose'][['reaction', 'measured', 'condition']]
                    .pivot(columns=['reaction'], values='measured', index='condition')
                    .reset_index(drop=True))

    exchange_fluxes = pd.read_excel(EXCHANGE_FLUX_FILE_PATH, sheet_name='Yields')
    exchange_fluxes = exchange_fluxes[exchange_fluxes.Strain == 'MG1655'][[col for col in exchange_fluxes.columns
                                                                           if not any([substr in col for substr in ['Yield', 'Strain', 'Reference']])]]

    all_exp_data = pd.concat([mfa_data_glc, exchange_fluxes])

    rxns_to_validate = [rxn for rxn in list(all_exp_data.columns)
                            if not any([substr in rxn for substr in ['Yield', '_ub', '+', 'lac']])]

    substrate_uptake_rates = list(exchange_fluxes['EX_glc__D_e']) + [mfa_data_glc['EX_glc__D_e'].iloc[0]]

    # perform_and_save_simulations(num_alternative_models=10,
    #                              substrate_uptake_rates=substrate_uptake_rates,)

    simulated_fluxes = pd.read_excel(RESULT_FLUX_PATH, sheet_name=None)

    all_fluxes = pd.concat([fluxes.assign(model=model) for model, fluxes in simulated_fluxes.items()])
    all_fluxes = all_fluxes[all_fluxes.rxn_id.isin(rxns_to_validate)]

    corr_coeff_df = determine_coeff_of_corr_sim_vs_measurements(all_exp_data, all_fluxes=all_fluxes)
    corr_coeff_df['diff_rsquared_intra'] = corr_coeff_df.mean_rsquared_intra - corr_coeff_df[corr_coeff_df.model == 'GotEnzymes'].mean_rsquared_intra.iloc[0]
    corr_coeff_df['diff_rsquared_extra'] =  corr_coeff_df.mean_rsquared_extra - corr_coeff_df[corr_coeff_df.model == 'GotEnzymes'].mean_rsquared_extra.iloc[0]

    print(corr_coeff_df.drop_duplicates(['model'])[['model', 'mean_rsquared_intra', 'mean_rsquared_extra']].to_latex())


    # fig, axs = plot_exchange_fluxes_vs_experiments(all_fluxes=all_fluxes,
    #                                                exchange_fluxes=exchange_fluxes,
    #                                                )
    # plt.show()
    # plt.savefig('Figures/SuppFig_ablation_analysis.png')

    # all_alternative_ids = [f'alternative {i}' for i in range(1,11)]
    #
    # cmap = {
    #     **{'iML1515': 'black', 'GotEnzymes': 'grey', 'Randomized': 'midnightblue'},
    #     **{f'Scaled {i}x': color for i, color in zip(list(range(6, 9)), sns.color_palette('plasma', n_colors=3))},
    #     **{alt: 'lightblue' for alt in all_alternative_ids},
    #     **{f'Randomized {i + 1}': 'lightcoral' for i in range(3)}, }
    #
    # fig, ax = plt.subplots()
    # for model, corr_df in corr_coeff_df.groupby('model'):
    #     if  not any ([substr in model for substr in ['alternative', 'Scaled']]): continue
    #     ax.scatter(corr_df['mean_rsquared_intra'], corr_df['mean_rsquared_extra'],
    #                color = cmap[model])
    #
    # ax.set_xlabel(r'R$^2$ intracellular fluxes')
    # ax.set_ylabel('R$^2$ extracellular fluxes')
    # plt.show()
    #


