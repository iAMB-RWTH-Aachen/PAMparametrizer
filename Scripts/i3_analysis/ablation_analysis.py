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
from PAModelpy import PAModel

from Modules.PAMparametrizer.utils.sector_config_functions import get_protein_sector_config
from Modules.PAMparametrizer.utils.pam_generation import create_pamodel_from_diagnostics_file
from Modules.PAMparametrizer.utils.pamparametrizer_visualization import RXN_NAME_MAPPER
from Modules.PAMparametrizer.utils.error_calculation import calculate_r_squared_for_reaction,calulate_pearson_correlation_simulation_vs_experiment

from Figures.Scripts.Figure4_sensitivity_error import *
from Scripts.pam_generation import setup_ecoli_pam as set_up_ecoli_pam_curated

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

RESULT_FLUX_PATH = os.path.join('Results', '3_analysis', 'iML1515_alternative_models_predictions_csources.xlsx')

def save_flux_predictions(model_files: List[str],
                         substrate_uptake_rates: List[float],
                         diagnostic_files: List[str],
                         label_names:List[str],
                          rxns_to_save: List[str],
                          substrate_uptake_id: str='EX_glc__D_e',
                         simulated_fluxes_result_file: str = RESULT_FLUX_PATH,
                         ):
    pam = set_up_pam(pam_info_file = PARAM_FILE_PREPROC, model = MODEL_FILE, sensitivity = False)
    models = [set_up_pam(pam_info_file = file, model = MODEL_FILE, sensitivity = False)
              if 'xml' not in file else read_sbml_model(file)
              for file in model_files
    ]
    alt_pams = [create_pamodel_from_diagnostics_file(diag_file, pam.copy_with_pickle()) for diag_file in diagnostic_files]
    curated_pam =  set_up_ecoli_pam_curated(
        pam_data_file_path=os.path.join('Data', 'proteinAllocationModel_iML1515_EnzymaticData_py.xls'),
        sensitivity=False)

    if substrate_uptake_id != 'EX_glc__D_e':
        trans_sector_config = get_protein_sector_config(pam,
                                                        substrate_id=substrate_uptake_id,
                                                        substrate_range=np.arange(-4, -1, 1)
                                                        )

        unused_sector_config = get_protein_sector_config(pam,
                                                        substrate_id=substrate_uptake_id,
                                                        substrate_range=np.arange(-4, -1, 1),
                                                         protein_sector='UnusedEnzymeSector'
                                                         )
    flux_df = None
    for label, model in zip(label_names+['Curated'], models+alt_pams+[curated_pam]):
        print('------------------------------------------------------------------------')
        print(f'Analyzing flux distribution with {substrate_uptake_id} of model', label, '\n')
        model.reactions.EX_glc__D_e.bounds = (0, 1e3)
        if substrate_uptake_id!='EX_glc__D_e' and isinstance(model, PAModel):
            if label == "Curated":
                usc = {'UnusedProteinSector':get_protein_sector_config(model,
                                                 substrate_id=substrate_uptake_id,
                                                 substrate_range=np.arange(-4, -1, 1),
                                                 protein_sector='UnusedProteinSector'
                                                 )}
            else: usc = {'UnusedEnzymeSector':unused_sector_config}
            for sectorid, config in {'TranslationalProteinSector': trans_sector_config, **usc}.items():
                change_sector_parameters_with_config_dict(pamodel=model,
                                                          sector_config=config,
                                                          substrate_uptake_id=substrate_uptake_id,
                                                          sector_id=sectorid)

        for substrate_uptake_rate in substrate_uptake_rates:
            model.reactions.get_by_id(substrate_uptake_id).lower_bound = substrate_uptake_rate
            model.optimize()
            if model.solver.status != 'optimal': continue

            fluxes_df = pd.Series({**{'model': label, 'substrate':substrate_uptake_rate},
                                   **{rxnid:model.reactions.get_by_id(rxnid).flux for rxnid in rxns_to_save if rxnid in model.reactions}}).to_frame().T
            kwargs = {'substrate_ids': [substrate_uptake_id],
                          'substrate_rates': [rate for rate in substrate_uptake_rates if rate == rate], #remove NaN, NaN does not equal itself
                          'fluxes_to_save': rxns_to_save
                          }
            if substrate_uptake_id != 'EX_glc__D_e':
                kwargs = {**kwargs,
                          'sectors_config':{
                              'TranslationalProteinSector': trans_sector_config,
                              'UnusedEnzymeSector': unused_sector_config
                          }}

            # fluxes_df = get_results_from_simulations(pam,**kwargs)['fluxes'].assign(model = label)
            if flux_df is not None:
                flux_df = pd.concat([flux_df, fluxes_df], ignore_index=True)
            else:
                flux_df = fluxes_df

    kwargs = {'mode': 'a', 'if_sheet_exists': 'replace'} \
        if os.path.exists(simulated_fluxes_result_file) else {'mode': 'w'}
    with pd.ExcelWriter(simulated_fluxes_result_file,
                        engine='openpyxl', **kwargs) as writer:
        flux_df.to_excel(writer, sheet_name=substrate_uptake_id, index= False)

def perform_and_save_simulations(num_alternative_models: int,
                                 substrate_uptake_rates: List[float],
                                 rxns_to_save: List[str],
                                 substrate_id: str='EX_glc__D_e',
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

    save_flux_predictions(model_files=[MODEL_FILE, PARAM_FILE_OLD, PARAM_FILE_PREPROC_6,
                                      PARAM_FILE_PREPROC, PARAM_FILE_PREPROC_8
                                       ]+RANDOM_PARAMETERS,
                          substrate_uptake_rates=substrate_uptake_rates,
                          substrate_uptake_id=substrate_id,
                          diagnostic_files=diagnostic_files,
                          label_names=all_model_labels,
                          rxns_to_save=rxns_to_save,
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
            **{'iML1515': 'black', 'GotEnzymes': 'black', 'Randomized': 'midnightblue', 'Curated': 'darkgrey'},
        **{f'Scaled {i}x': color for i, color in zip(list(range(6,9)), sns.color_palette('plasma', n_colors=3))},
        **{alt: 'lightblue' for alt in all_alternative_ids},
            **{f'Randomized {i+1}':'lightcoral' for i in range(3)},}
    linestyles = {**{alt_id:'-' for alt_id in all_alternative_ids},
                  **{f'Randomized {i}':":" for i in range(1,4)},
                  **{'iML1515': '--', 'GotEnzymes': '-.', 'Randomized': ':', 'After preprocessing':'-', 'Curated':':'},
                  **{f'Scaled {i}x':'--' for i in range(6,9)},}
    model_order = (['iML1515', 'GotEnzymes', 'Randomized', 'After preprocessing', 'Curated']
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
                            if not any([substr in rxn for substr in ['Yield', '_ub', '+', 'EX_lac_e', 'EX_', 'substrate','EX_glc__D']])
                        ]
    rxns_to_validate_extracell = [rxn for rxn in list(all_exp_data.columns)
                                  if not any([substr in rxn for substr in ['Yield', '_ub', '+', 'EX_lac_e', 'substrate', 'EX_glc__D']])
                                  and 'EX_' in rxn
                                  ]

    for (substrate_id, model), flux_df in all_fluxes.groupby(['substrate_id', 'model']):
        flux_df_wide =(flux_df.drop_duplicates(['rxn_id', 'substrate'])
                       .pivot(columns='rxn_id', values='flux', index = 'substrate')
                       .reset_index()
                       # .rename(columns={'substrate_uptake_rate':'substrate'})
                       .dropna(axis =1, how='all')
                       .assign(substrate_id = substrate_id)
                       .sort_values(['substrate'])
                       )
        exp_data_csource = all_exp_data[all_exp_data['substrate_id'] == substrate_id].sort_values(substrate_id)
        r_squared_intracell, rxns_intracell  = calculate_error_for_reactions(validation_df=exp_data_csource,
                                                   flux_df=flux_df_wide,
                                                             substrate_id=substrate_id,
                                                   rxns_to_validate=[rxn for rxn in rxns_to_validate_intracell if rxn in flux_df_wide.columns],)
        r_squared_extracell, rxns_extracell  = calculate_error_for_reactions(validation_df=exp_data_csource,
                                                   flux_df=flux_df_wide,
                                                             substrate_id=substrate_id,
                                                   rxns_to_validate=[rxn for rxn in rxns_to_validate_extracell if rxn in flux_df_wide.columns],)
        r_squared = r_squared_intracell + r_squared_extracell
        pearson_intra,pval_intra = calulate_pearson_correlation_simulation_vs_experiment(validation_df=exp_data_csource,
                                                                              flux_df=flux_df_wide,
                                                                              rxns_to_validate=[rxn for rxn in rxns_to_validate_intracell if rxn in flux_df_wide.columns],
                                                                              substr_rxn=substrate_id)
        pearson_extra,pval_extra = calulate_pearson_correlation_simulation_vs_experiment(validation_df=exp_data_csource,
                                                                              flux_df=flux_df_wide,
                                                                              rxns_to_validate=[rxn for rxn in rxns_to_validate_extracell if rxn in flux_df_wide.columns],
                                                                              substr_rxn=substrate_id)
        for i, rxn in enumerate(rxns_intracell+rxns_extracell):
            errors_for_models_rows.append({'model': model,
                                           'substrate_id': substrate_id,
                                           'rxn_id':rxn,
                                           'r_squared': r_squared[i],
                                           'mean_rsquared_intra': np.mean(r_squared_intracell),
                                           'mean_rsquared_extra': np.mean(r_squared_extracell),
                                           'pearson_intra': pearson_intra,
                                           'pearson_extra': pearson_extra,
                                           'pval_intra': pval_intra,
                                           'pval_extra': pval_extra,
                                           })
    return pd.DataFrame(errors_for_models_rows)

def calculate_error_for_reactions(validation_df: pd.DataFrame,
                                   flux_df: pd.DataFrame,
                                  substrate_id: str,
                                  rxns_to_validate: list) -> Tuple[List[float], List[str]]:
    # calculate error for different exchange rates
    error = []
    rxns = []
    for rxn in rxns_to_validate:
        # only select the rows which are filled with data
        validation_data = validation_df.dropna(axis=0, subset=rxn)
        # if there are no reference data points, continue to the next reaction
        if len(validation_data) == 0:
            continue
        rsquareds = []
        for i, validation in validation_data.iterrows():
            rate = validation[substrate_id]

            r_squared = calculate_r_squared_for_reaction(rxn, validation_df, substrate_id,
                                                               flux_df[flux_df.substrate == rate])
            rsquareds.append(r_squared)
        rxns += [rxn]
        error += [np.mean(rsquareds)]
    return error, rxns

def get_fluxomics_data_multiple_csources():

    flux_csources = pd.read_excel(fluxomics_data_file,
                                  sheet_name='Fluxes_Csources',
                                  index_col=1)
    flux_csources_df = flux_csources.drop(['Flux (publication)', 'Reversibility'], axis=1)
    fluxes_to_simulate = flux_csources_df.copy()
    new_indices = []
    for i, row in fluxes_to_simulate.iterrows():
        if isinstance(row.name, str):
            if row.name[-2:] == '_b':
                new_indices.append(row.name[:-2])
                fluxes_to_simulate.loc[row.name] = -row
            else:
                new_indices.append(row.name)
        else:
            new_indices.append(row.name)

    fluxes_to_simulate.index = new_indices
    fluxes_to_simulate_parsed = fluxes_to_simulate[fluxes_to_simulate.index.notnull()]
    fluxes_to_simulate_parsed = fluxes_to_simulate_parsed.rename(
        index={'BIOMASS_Ec_iML1515_WT_75p37M': 'BIOMASS_Ecoli_core_w_GAM'}
    ).drop('Glucose (flux ratio Glc)', axis=1)    # substrate_rates = [

    flux_mapper = {col: fluxes_to_simulate_parsed.index[i] for i, col in enumerate(fluxes_to_simulate_parsed.columns)}
    fluxes_to_save = []
    # Get the fluxes we want to save
    for flux in fluxes_to_simulate_parsed.index:
        if isinstance(flux, str):
            fluxes_to_save += [f for f in flux.split(', ')]

    # parse the fluxes such that we can run and validate simulations easily
    validation_df = pd.DataFrame(columns=list(fluxes_to_simulate_parsed.index))
    substrate_ids = []
    substrate_uptake = []
    for substrate, fluxes in fluxes_to_simulate_parsed.items():
        substrate_ids += [flux_mapper[substrate]]
        substrate_uptake += [fluxes.loc[flux_mapper[substrate]]]
        validation_df = pd.concat([validation_df, fluxes.to_frame().T], ignore_index=True)

    validation_df.index = list(flux_mapper.values())
    return validation_df, substrate_ids, substrate_uptake, fluxes_to_save

if __name__ == '__main__':
    mfa_data = pd.read_excel(MFA_DATA_FILE_PATH)
    mfa_data_glc = (mfa_data[mfa_data.condition == 'Glucose'][['reaction', 'measured', 'condition']]
                    .pivot(columns=['reaction'], values='measured', index='condition')
                    .reset_index(drop=True))
    fluxomics_data_file = os.path.join('Data', 'Ecoli_phenotypes', 'Ecoli_phenotypes_py.xls')
    flux_data_glc = get_fluxomics_data(fluxomics_data_file)
    rxns_to_save_glc, valid_df_glucose = get_reactions_to_save(flux_data_glc)
    valid_df_glucose['EX_glc__D_e'] = -valid_df_glucose['EX_glc__D_e']

    exchange_fluxes = pd.read_excel(EXCHANGE_FLUX_FILE_PATH, sheet_name='Yields')
    exchange_fluxes = exchange_fluxes[exchange_fluxes.Strain == 'MG1655'][[col for col in exchange_fluxes.columns
                                                                           if not any([substr in col for substr in ['Yield', 'Strain', 'Reference']])]]

    all_exp_data_glc = pd.concat([valid_df_glucose, exchange_fluxes], ignore_index=True).assign(substrate_id = 'EX_glc__D_e')

    substrate_uptake_rates_glc = list(all_exp_data_glc['EX_glc__D_e'].dropna().drop_duplicates())

    valid_df_other_csources, substrate_ids, substrate_uptake_rates_csource, fluxes_to_save_csources = get_fluxomics_data_multiple_csources()
    # for substrate, rate in zip(substrate_ids, substrate_uptake_rates_csource):
    #     perform_and_save_simulations(num_alternative_models=10,
    #                                       substrate_uptake_rates=[rate],
    #                                       substrate_id=substrate,
    #                                       rxns_to_save=fluxes_to_save_csources,)
    #
    #
    # perform_and_save_simulations(num_alternative_models=10,
    #                                           substrate_uptake_rates=substrate_uptake_rates_glc,
    #                                           rxns_to_save=rxns_to_save_glc,)
    #
    simulated_fluxes = pd.read_excel(RESULT_FLUX_PATH, sheet_name=None)
    #
    # for substrate, fluxes in simulated_fluxes.copy().items():
    #     if substrate != 'EX_glc__D_e':
    #         fluxes['substrate'] = valid_df_other_csources.at[substrate, substrate]
    #     simulated_fluxes[substrate] = fluxes

    all_fluxes = pd.concat([fluxes.assign(substrate_id=substrate) for substrate, fluxes in simulated_fluxes.items()], ignore_index=True)
    all_fluxes = all_fluxes.melt(id_vars=['substrate_id', 'substrate', 'model'], var_name='rxn_id', value_name='flux')
    # all_fluxes = all_fluxes[all_fluxes.rxn_id.isin(rxns_to_validate)]
    all_exp_data = pd.concat([all_exp_data_glc, valid_df_other_csources.reset_index(names = 'substrate_id')], ignore_index=True)

    corr_coeff_df = determine_coeff_of_corr_sim_vs_measurements(all_exp_data, all_fluxes=all_fluxes)
    corr_coeff_df_glc = corr_coeff_df[corr_coeff_df.substrate_id == 'EX_glc__D_e']
    corr_coeff_df_other_csource = corr_coeff_df[corr_coeff_df.substrate_id != 'EX_glc__D_e']
    for col in ['pearson_intra', 'pearson_extra']:
        column = {}

        for model, df in corr_coeff_df_other_csource.groupby('model'):

            column[model] = df[col].mean()
            # print(model, df.mean_rsquared_intra.mean(), df.mean_rsquared_extra.mean())
            # corr_coeff_df_glc[]
        corr_coeff_df_glc = pd.merge(corr_coeff_df_glc, pd.Series(column, name=col+'_csource').to_frame().reset_index(names = 'model'), on='model')

    for csource in ['', '_csource']:
        for type in ['pearson_intra', 'pearson_extra']:
            ref_val = corr_coeff_df_glc.loc[
                corr_coeff_df_glc["model"] == "GotEnzymes", f"{type}{csource}"
            ].iloc[0]
            corr_coeff_df_glc[f'diff_{type}{csource}'] = corr_coeff_df_glc[type] - ref_val
    pd.set_option('display.float_format', '{:.2g}'.format)
    print(corr_coeff_df_glc.drop_duplicates(['model', 'substrate_id'])[
              ['model'] + [col for col in corr_coeff_df_glc.columns if ('diff_pearson' in col) or ('pval' in col)]
          ].to_latex(index=False, na_rep='-')
          )
    print(corr_coeff_df_glc.drop_duplicates(['model', 'substrate_id'])[
              ['model'] + [col for col in corr_coeff_df_glc.columns if col[:3] == 'pea']
              ].to_latex(index=False, na_rep='-')
          )

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


