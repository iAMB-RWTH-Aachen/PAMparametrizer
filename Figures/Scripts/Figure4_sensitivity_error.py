import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
from cobra.io.sbml import read_sbml_model
from matplotlib import gridspec
import matplotlib.colors as mcolors
import seaborn as sns
from scipy.cluster.hierarchy import leaves_list
from scipy.stats import mannwhitneyu
from typing import List, Dict, Tuple, Iterable, Union

from PAModelpy.utils import set_up_pam
from PAModelpy import PAModel

from Scripts.pam_generation import setup_ecoli_pam as set_up_ecoli_pam_curated
from Scripts.i3_analysis.metabolic_flux_distribution_vs_exp import RXNS_TO_VALIDATE, get_reactions2plot_pathway_mapping
from Modules.utils.sector_config_functions import change_sector_parameters_with_config_dict

from Modules.utils.pamparametrizer_analysis import (get_results_from_simulations,
                                                    calculate_error_for_reactions, calculate_r_squared_for_reaction,
                                                    calculate_difference_simulation_experiment)
from Modules.utils.pam_generation import (create_pamodel_from_diagnostics_file,
                                          _extract_reaction_id_from_catalytic_reaction_id)

# from Modules.utils import calculate_r_squared_for_reaction
# from Scripts.Visualization.PAMparametrizer_progress_cleaned_figure import run_simulations

def make_simulation_error_boxplot(gotenz_param_file:str,
                                  preprocessed_param_file:str,
                                  diagnostic_files: List[str],
                                  model_file_path:str,
                                  fluxomics_data_file:str,
                                  ax:plt.Axes,
                                  substrate_uptake_ids: List[str] = ['EX_glc__D_e'],
                                  fontsize:int = 16):
    ecoli_pams = set_up_pams_different_parameters(gotenz_param_file,
                                     preprocessed_param_file,
                                     diagnostic_files,
                                     model_file_path)
    ecoli_pams['iML1515'] = read_sbml_model(os.path.join('Models', 'iML1515.xml'))
    flux_data = get_fluxomics_data(fluxomics_data_file)
    rxns_to_save, valid_df = get_reactions_to_save(flux_data)
    substrate_rates = [
        [-f for f in flux_data.loc[substr_id, :].values] for substr_id in substrate_uptake_ids
    ]
    simulated_fluxes = get_simulation_results_for_models(ecoli_pams,
                                                         substrate_rates = substrate_rates,
                                                         substrate_uptake_ids = substrate_uptake_ids,
                                                         rxns_to_save = rxns_to_save)
    calculate_error_for_simulations(simulated_fluxes, valid_df, rxns_to_save)
    return ecoli_pams, create_simulation_error_boxplot(simulated_fluxes,
                                           valid_df,
                                           rxns_to_save,
                                           ax, fontsize)



def set_up_pams_different_parameters(gotenz_param_file:str,
                                     preprocessed_param_file:str,
                                     diagnostic_files: List[str],
                                     model_file_path:str):
    # setup the model
    ecoli_pam_wt = set_up_pam(gotenz_param_file,
                              model=model_file_path,
                              sensitivity=False)  # not curation for reference
    ecoli_pam_curated = set_up_ecoli_pam_curated(
        pam_data_file_path=os.path.join('Data', 'proteinAllocationModel_iML1515_EnzymaticData_py.xls'),
        sensitivity=False)  # curated for reference

    pam = set_up_pam(preprocessed_param_file,
                     model=model_file_path,
                     sensitivity=False)
    pam.change_reaction_bounds('EX_glc__D_e', lower_bound=0)

    new_ecoli_pams = {f"Alternative {alt+1}": create_pamodel_from_diagnostics_file(file,
                                                                    pam.copy(copy_with_pickle=True)) for alt, file in
                      enumerate(diagnostic_files)}

    return {'GotEnzymes': ecoli_pam_wt, 'Curated': ecoli_pam_curated, **new_ecoli_pams}

def get_fluxomics_data(fluxomics_data_file_path: str)->pd.DataFrame:
    def process_index(name):
        if isinstance(name, str):
            if name.endswith('_b'):
                return name[:-2]  # Remove '_b' from index
            elif name.endswith('_f'):
                return name[:-2]  # Remove '_f' from index
        return name  # Keep unchanged if no suffix
    flux_df = pd.read_excel(fluxomics_data_file_path,
                            sheet_name='Fluxes',
                            engine='openpyxl',
                            index_col=1)

    flux_df.index = flux_df.index.map(process_index)
    return flux_df.drop(columns=['Description', 'Reaction_position'])

def get_reactions_to_save(flux_data: pd.DataFrame) -> Tuple[List[str], pd.DataFrame]:
    fluxes_to_save = None
    flux_information = {}
    validation_df = pd.DataFrame(columns=flux_data.index)
    studies = []
    substrate_uptake = []
    for study, fluxes in flux_data.items():
        studies += [study]
        substrate_uptake += [-fluxes.loc['EX_glc__D_e']]
        flux_information[study] = fluxes.to_dict()
        validation_df = pd.concat([validation_df, fluxes.to_frame().T], ignore_index=True)

        # store the names of the fluxes to save
        if fluxes_to_save is None:
            fluxes_to_save = list(fluxes.index)
    #parse validation_data in proper format for calculating differences
    #validation_df = validation_df.T.reset_index()
    #print(validation_df)
    # validation_df['index'] = validation_df['index'].str.split(', ')
    # validation_df = validation_df.explode('index').set_index('index').T
    #validation_df = validation_df.set_index('Reaction_ID').T

    return fluxes_to_save, validation_df

def get_simulation_results_for_models(ecoli_models:Dict[str, Union['PAModel', 'Model']],
                                      substrate_rates:list[Iterable],
                                      substrate_uptake_ids:list[str],
                                      rxns_to_save: List[str]) -> Dict[str, pd.DataFrame]:
    kwargs = {'substrate_ids': substrate_uptake_ids,
          'substrate_rates': substrate_rates,
          'fluxes_to_save' : rxns_to_save}
    simulated_fluxes = {}
    for label, model in ecoli_models.items():
        if isinstance(model, PAModel): model.change_reaction_bounds('EX_glc__D_e', 0, 1e3)
        else: model.reactions.EX_glc__D_e.bounds = (0,1e3)
        simulated_fluxes[label] = get_results_from_simulations(model, **kwargs)['fluxes']
    return simulated_fluxes

def calculate_error_for_simulations(simulation_results:Dict[str, pd.DataFrame],
                                    validation_df: pd.DataFrame,
                                    fluxes_to_save: List[str]):

    error_new_dict = {alt: calculate_error_for_reactions(validation_df,
                                                         fluxes,
                                                         fluxes_to_save[1:]) for alt, fluxes in simulation_results.items()}
    for alt, error_list in error_new_dict.items():
        print(f'R^2 values for alternative model {alt} with the optimized parameters: ', np.nanmean(error_list))

def create_simulation_error_boxplot(simulated_fluxes,
                                    valid_df: pd.DataFrame,
                                    fluxes_to_save,
                                    ax: plt.Axes,
                                    fontsize: int,
                                    other_colors={'GotEnzymes': 'grey', 'After preprocessing': 'black',
                                                  'Curated': 'chocolate', 'iML1515': 'white'}
                                    ) -> plt.Axes:
    models = [label for label in simulated_fluxes if label not in other_colors.keys()]
    other_colors = {k: v for k, v in other_colors.items() if k in simulated_fluxes}
    model_colors = sns.color_palette("winter", n_colors=len(simulated_fluxes) - len(other_colors))
    cmap = {**dict(zip(models, model_colors)), **other_colors}

    # Combine data into a DataFrame
    all_differences = pd.DataFrame()
    curated_differences = None  # Placeholder for Curated errors
    df_list = []
    for model, flux_df in simulated_fluxes.items():

        differences = []
        for _, row in flux_df.iterrows():
            substrate_id = row['substrate_id']
            difference = calculate_difference_simulation_experiment(
                valid_df, row, fluxes_to_save[1:], substrate_id)
            differences += difference
        df_list.append({
            'model': model,
            'median' :np.median(differences),
            'mean' : np.mean(differences),
            'stdev' : np.std(differences)
        })

        temp_df = pd.DataFrame({'Model': [model] * len(differences), 'Difference': differences})
        if model == 'GotEnzymes':
            curated_differences = differences
        elif curated_differences is not None:
            # Statistical test
            stat, p = mannwhitneyu(curated_differences, differences, alternative='greater')
            print(
                f"{model}: U-statistic = {stat}, p-value = {p}, median = {np.median(differences)}, mean = {np.mean(differences)}")

        # Append to the main DataFrame
        all_differences = pd.concat([all_differences, temp_df], ignore_index=True)

    # Boxplot or Violin Plot
    sns.boxplot(x='Model', y='Difference', data=all_differences, ax=ax, palette=cmap, showfliers=False)

    # Adjust y-axis to focus on bulk data and add space for annotations
    #     ax.set_ylim([all_differences['Difference'].quantile(0.05), all_differences['Difference'].quantile(0.96)])

    # Set labels and title
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right', fontsize=fontsize)
    ax.set_xlabel('Model', fontsize=fontsize)
    ax.set_ylabel(r'Difference (exp - sim) mmol/$\text{g}_{CDW}$/h', fontsize=fontsize)
    df = pd.DataFrame(df_list, columns=['model', 'median', 'mean', 'stdev'])
    print(df.to_latex())
    return ax

def create_sensitivity_heatmap(ecoli_pams: dict, gs, fig, fontsize):
    etc_reactions_proteins = get_reactions_from_etc(ecoli_pams['Alternative 1'], ecoli_pams['Curated'])
    enzyme_sensitivities = get_most_sensitive_enzymes(ecoli_pams, etc_reactions_proteins)
    plot_split_clustermap(enzyme_sensitivities, etc_reactions_proteins,gs, fig, fontsize)

def get_reactions_from_etc(pam, pam_curated):
    etc_reactions = pd.DataFrame({'Reactions':
                                      ['SUCDi',  # succinate dehydrogenase
                                       'NADH16pp', 'NADH17pp', 'NADH18pp', 'NADH10', 'NADH9', 'NADH5',
                                       # NADH dehydrogenase
                                       'ATPS4rpp',  # ATP synthase
                                       'CYTBDpp', 'CYTBD2pp']})  # cytochrome B oxidase

    etc_reactions_proteins = etc_reactions.copy()
    etc_reactions_proteins['enzyme_id'] = etc_reactions.Reactions.apply(
        lambda x: (pam.get_enzymes_with_reaction_id(x)  # uniprot ids, current pam version
                   + pam_curated.get_enzymes_with_reaction_id(x))  # ec numbers, previous pam version, curated
    )
    etc_reactions_proteins = etc_reactions_proteins.explode('enzyme_id')
    return etc_reactions_proteins

def get_most_sensitive_enzymes(ecoli_pams:dict, etc_reactions_proteins:pd.DataFrame):
    # get the most sensitive enzymes
    top_enzymes = pd.DataFrame(columns=['rxn_id', 'enzyme_id', 'coefficient', 'model'])
    escs = {}
    top_esc = {}

    for label, pam in ecoli_pams.items():
        if not isinstance(pam, PAModel): continue
        pam.change_reaction_bounds('EX_glc__D_e', -10, 0)
        pam.sensitivity = True
        pam.optimize()
        print('Optimizing the following model: ', label)
        print('Growth rate is ', pam.objective.value)

        # save top escs
        esc = pam.enzyme_sensitivity_coefficients.copy()
        esc_df = pam.enzyme_sensitivity_coefficients.sort_values('coefficient', ascending=False).iloc[:6]
        top_esc[label] = esc_df[esc_df.coefficient > 0.05]
        #parse esc
        esc['rxn_id'] = esc.rxn_id.str.split(',')
        esc = esc.explode('rxn_id')
        # extract actual reaction ids
        esc['Reactions'] = esc.rxn_id.apply(lambda x: _extract_reaction_id_from_catalytic_reaction_id(x))

    top_enzymes = pd.DataFrame(columns=['rxn_id', 'enzyme_id', 'coefficient', 'model'])
    for label, esc in top_esc.items():
        esc['model'] = label
        top_enzymes = pd.concat([top_enzymes, esc])
    top_enzymes = top_enzymes.drop_duplicates(['rxn_id', 'enzyme_id']).reset_index()

    for label, esc in escs.items():
        #get only interesting proteins
        esc = pd.merge(esc, etc_reactions_proteins, how='inner', on=['Reactions', 'enzyme_id'])
        # only get the sensitivities for the enzymes associated with energy generation
        esc = pd.merge(esc[['enzyme_id']],
                       pam.enzyme_sensitivity_coefficients,
                       on='enzyme_id', how='inner'
                       ).drop_duplicates()
        esc['model'] = label
        escs[label] = esc
        top_enzymes = pd.concat([top_enzymes, esc])

        # get all the enzymes which are important in the models
        esc_enzymes = pd.merge(pam.enzyme_sensitivity_coefficients,
                               top_enzymes[['enzyme_id']],
                               on='enzyme_id', how='inner')

        esc_enzymes['model'] = label
        top_enzymes = pd.concat([top_enzymes, esc_enzymes])

    #parse for plotting
    enzyme_sensitivities = top_enzymes.drop_duplicates(['model', 'enzyme_id'])
    enzyme_sensitivities['rxn_id'] = enzyme_sensitivities['rxn_id'].str.split(',')
    enzyme_sensitivities = enzyme_sensitivities.explode('rxn_id')
    enzyme_sensitivities['reaction'] = enzyme_sensitivities['rxn_id'].apply(
        lambda x: _extract_reaction_id_from_catalytic_reaction_id(
            x,
            #         default_enzyme_id_pattern=r'E[0-9][0-9]*|Enzyme_*|\d{1,2}(\.\d{0,2}){0,3}'#for ec numbers
        )
    )
    return enzyme_sensitivities


def plot_split_clustermap(enzyme_sensitivities, genes_reactions_etc, gs, fig, fontsize):
    # Split data into two subsets
    mask = enzyme_sensitivities.reaction.isin(genes_reactions_etc.Reactions)
    heatmap_1 = enzyme_sensitivities[mask]
    heatmap_2 = enzyme_sensitivities[~mask]

    # Pivot data to create heatmap format
    heatmap_1 = pd.pivot_table(
        heatmap_1.drop_duplicates(['reaction', 'enzyme_id', 'model']),
        index='reaction', columns='model', values='coefficient', aggfunc="sum"
    ).dropna(how='all').fillna(0)

    heatmap_2 = pd.pivot_table(
        heatmap_2.drop_duplicates(['reaction', 'enzyme_id', 'model']),
        index='enzyme_id', columns='model', values='coefficient', aggfunc="sum"
    ).dropna(how='all').fillna(0)

    # Replace enzyme IDs with reaction IDs
    enzyme_to_reaction = enzyme_sensitivities.drop_duplicates(
        'enzyme_id'
    ).set_index('enzyme_id')[['reaction']]
    heatmap_2 = pd.merge(
        heatmap_2, enzyme_to_reaction, left_index=True, right_index=True
    ).reset_index().set_index('reaction').drop('enzyme_id', axis=1)

    # Define custom colormap
    colors_pos = plt.cm.coolwarm(np.linspace(0, 1, 256))
    colors_zero = np.array([[0.549, 0.5725, 1, 0.9]])
    colors = np.vstack((colors_zero, colors_pos))
    combined_cmap = mcolors.ListedColormap(colors, name='custom_cmap')

    vmin, vmax = 0, 1  # round(max(heatmap_1.max().max(), heatmap_2.max().max()))
    bounds = np.linspace(vmin, vmax, len(colors))
    norm = mcolors.BoundaryNorm(bounds, combined_cmap.N)

    # 2 columns: heatmap and cmap
    gs_main = gridspec.GridSpecFromSubplotSpec(1, 2, width_ratios=[10, 1],
                                               wspace=0.4, subplot_spec=gs)
    # 2 rows: ETC heatmap and top sensitivities heatmap
    gs_inner = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs_main[0],
                                                height_ratios=[heatmap_1.shape[0], heatmap_2.shape[0]],
                                                hspace=0.05)

    # Heatmap 1 (ETC)
    ax1 = fig.add_subplot(gs_inner[0, 0])
    sns.heatmap(heatmap_1, cmap=combined_cmap, norm=norm, cbar=False, ax=ax1)
    ax1.xaxis.set_visible(False)
    ax1.set_xlabel("")
    ax1.set_ylabel("ETC", fontsize=fontsize)

    # Step 2: Extract reordered row indices
    rxns_to_plot = get_reactions2plot_pathway_mapping(heatmap_2.T)
    # Flatten all reaction lists into a set for fast lookup
    all_rxns = set(rxn for rxn_list in rxns_to_plot.values() for rxn in rxn_list)

    # Loop over indices and add to 'anabolism' if not present in any of the lists
    for idx in heatmap_2.index:
        if idx not in all_rxns:
            rxns_to_plot.setdefault('Other', []).append(idx)

    row_order = []  # leaves_list(cluster_map.dendrogram_row.linkage)  # Get sorted row indices
    for rxns in rxns_to_plot.values(): row_order += rxns
    data_reordered = heatmap_2.loc[row_order]  # Reorder data

    # Step 4: Add heatmap to GridSpec
    ax2 = fig.add_subplot(gs_inner[1, 0])
    sns.heatmap(data_reordered, cmap=combined_cmap, norm=norm, cbar=False, ax=ax2)
    ax2.set_ylabel("Reaction", fontsize=fontsize)

    ax2.set_xlabel("")
    ax2.set_xticklabels(ax2.get_xticklabels(), rotation=45, ha='right')

    # add pathway annotation
    group_labels = []
    group_positions = []
    start = 0
    for group_name, row in rxns_to_plot.items():
        length = len(row)
        center = start + length / 2
        group_labels.append(group_name)
        group_positions.append(center)
        start += length
        # drawing line between the groups
        if start < heatmap_2.shape[0]:  # Avoid drawing a line at the far right
            ax2.axhline(y=start, color='black', linewidth=2)

    # Add group labels as a second y-axis (on the left)
    ax_top = ax2.twinx()
    # Match the ticks and limits
    ax_top.set_ylim(ax2.get_ylim())
    ax_top.set_yticks(group_positions)
    ax_top.set_yticklabels(group_labels, fontsize=fontsize, fontweight='bold')
    ax_top.tick_params(axis='y', bottom=False, top=True, labelbottom=False, labeltop=True)
    ax_top.spines['bottom'].set_visible(False)

    # Colorbar spanning both heatmaps
    cbar_ax = fig.add_subplot(gs_main[0, 1])
    sm = plt.cm.ScalarMappable(cmap=combined_cmap, norm=norm)
    cbar = fig.colorbar(sm, cax=cbar_ax)
    # Define tick positions at whole and half numbers
    tick_positions = np.arange(0, vmax + 0.5, 0.5)  # Adjust `vmax` as needed
    valid_ticks = [t for t in tick_positions if vmin <= t <= vmax]
    cbar.ax.yaxis.set_minor_locator(plt.NullLocator())  # Remove minor ticks
    cbar.set_ticks(valid_ticks)

    cbar.ax.set_ylabel("Sensitivity Coefficient", fontsize=fontsize)



def main():
    N_ALT_MODELS = 10
    FONTSIZE = 16

    ECOLI_PHENOTYPE_DATA_PATH = os.path.join('Data', 'Ecoli_phenotypes', 'Ecoli_phenotypes_py.xls')

    MODEL_FILE_PATH = os.path.join('Models', 'iML1515.xml')

    PARAM_FILE_GOTENZ = os.path.join('Results', '1_preprocessing', 'proteinAllocationModel_iML1515_EnzymaticData_250912.xlsx')
    PARAM_FILE_PREPROC = os.path.join('Results', '2_parametrization',
                                     'proteinAllocationModel_iML1515_EnzymaticData_multi.xlsx')

    BEST_INDIV_RESULT_FILES = [os.path.join('Results', '2_parametrization', 'diagnostics',
                                        f'pam_parametrizer_diagnostics_{i}.xlsx') for i in range(1, N_ALT_MODELS + 1)]

    fig = plt.figure(figsize=(15, 40))

    gs_main = gridspec.GridSpec(2, 1, height_ratios=[5, 6], wspace=0)
    ax1 = fig.add_subplot(gs_main[0])

    ecoli_pams, ax1 = make_simulation_error_boxplot(PARAM_FILE_GOTENZ,
                                        PARAM_FILE_PREPROC,
                                        BEST_INDIV_RESULT_FILES,
                                        MODEL_FILE_PATH,
                                        ECOLI_PHENOTYPE_DATA_PATH,
                                        ax1, fontsize= FONTSIZE)
    create_sensitivity_heatmap(ecoli_pams, gs_main[1], fig, FONTSIZE)

    for ax in fig.axes:
        ax.tick_params(axis='both', labelsize=FONTSIZE)

    annotations = ["A", "", "B"]
    fontsize = FONTSIZE  # Adjust as needed

    for ax, label in zip(fig.axes, annotations):
        ax.annotate(label, xy=(0, 1), xycoords="axes fraction",
                    fontsize=fontsize, fontweight='bold',
                    xytext=(-5, 5), textcoords="offset points",
                    ha="right", va="bottom")

    fig.savefig(os.path.join('Figures', 'Figure2_sensitivity_simerror.png'))

if __name__ == '__main__':
    main()