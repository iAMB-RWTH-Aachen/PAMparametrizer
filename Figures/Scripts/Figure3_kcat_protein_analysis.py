import matplotlib.pyplot as plt
from cobra.io.sbml import read_sbml_model
from matplotlib.colors import to_hex
from matplotlib import gridspec
import seaborn as sns
import numpy as np
import pandas as pd
import os
from typing import List, Dict, Tuple, Callable
from PAModelpy.utils import set_up_pam

from Scripts.i2_parametrization.pam_parametrizer_iML1515 import set_up_pamparametrizer

from Modules.PAMparametrizer.utils.pam_generation import (get_rxn2kcat_protein2gene_dict,
                                          _extract_reaction_id_from_catalytic_reaction_id,
                                          _get_rxn2kcat_as_series,
                                          create_pamodel_from_diagnostics_file
                                         )
from Modules.PAMparametrizer.utils.pamparametrizer_analysis import (calculate_kcat_differences,
                                                    plot_histogram_logspace,
                                                    select_clustered_rows_by_variation,
                                                    get_clusters_from_clustermap,
                                                    set_up_pam_parametrizer_and_get_substrate_uptake_rates,
                                                                    get_results_from_simulations_fixed_mu,
                                                    parse_enzyme_complex_id,
                                                    convert_peptide_to_enzyme_concentrations,
                                                    normalize_simulated_protein_concentrations
                                                   )

from Modules.PAMparametrizer.utils.pamparametrizer_visualization import plot_valid_data, plot_simulation, plot_flux_vs_experiment
from Modules.PAMparametrizer.utils.pamparametrizer_setup import set_up_sector_config_from_diagnostic_file

COG_MAPPER = {'Amino acid transport and metabolism': 'Amino acid metabolism',
       'Carbohydrate transport and metabolism': 'Carbon metabolism',
       'Cell cycle control, cell division, chromosome partitioning': 'Cell cycle',
       'Cell wall/membrane/envelope biogenesis':"Cell membrane biogenesis",
       'Coenzyme transport and metabolism': 'Coenzyme metabolism',
              'Defense mechanisms':'Defense mechanisms',
       'Energy production and conversion': 'Energy generation', 'Function unknown': 'Function unknown',
       'General function prediction only': 'General function',
       'Inorganic ion transport and metabolism': 'Inorganic ion metabolism',
       'Intracellular trafficking, secretion, and vesicular transport': 'Intracellular transport',
       'Lipid transport and metabolism': 'Lipid metabolism',
       'Nucleotide transport and metabolism': 'Nucleotide metabolism',
       'Posttranslational modification, protein turnover, chaperones':'Protein modification',
       'Replication, recombination and repair': 'Replication',
       'Secondary metabolites biosynthesis, transport and catabolism': 'Secondary metabolites metabolism',
       'Signal transduction mechanisms': 'Signal transduction', 'Transcription': 'Transcription',
       'Translation, ribosomal structure and biogenesis': 'Translation'}


def gaussian(x, mean, amplitude, standard_deviation):
    return amplitude * np.exp( - (x - mean)**2 / (2*standard_deviation ** 2))

def calculate_distribution_statistics(bin_heigths: list[float],
                                      bin_borders:list[float]) -> None:
    peak_bin_index = np.argmax(bin_heigths)
    peak_bin_value = (bin_borders[peak_bin_index] + bin_borders[peak_bin_index + 1]) / 2
    area_under_curve = sum(
        [(abs(bin_borders[i]) - abs(bin_borders[i + 1])) * bin_heigths[i] for i in range(len(bin_heigths))])

    print(
        f'\tMost frequent flux:\t\t{peak_bin_value} mmol/gCDW/h\n\tArea under the curve:\t\t{area_under_curve} mmol/gCDW/h\n')

def create_kcat_histogram_old_vs_new(data_file_paths: list[pd.DataFrame],
                                     ax:plt.Axes,
                                     label_names:list[str],
                                     cumulative: bool = False,
                                     other_colors = {'GotEnzymes': 'grey', 'After preprocessing': 'black'},
                                     legend = True,
                                     fontsize =16,
                                     cmap: dict=None):
    # fig, ax = plt.subplots()
    n_bins = 50
    i = 0
    if cmap is None: cmap = plt.get_cmap("coolwarm")

    for label, data_file_path in zip(label_names, data_file_paths):
        aes_parameter_df = pd.read_excel(data_file_path, sheet_name='ActiveEnzymes')
        kcat_values = aes_parameter_df.kcat_values.dropna()
        print('------------------------------------------------------------------------')
        print(f'The kcat set from {label} has:\n \tMedian:\t\t\t{np.median(kcat_values)} \n \tMean:\t\t\t{np.mean(kcat_values)}')


        hist, bins = np.histogram(kcat_values, bins=n_bins)
        logbins = np.logspace(np.log10(bins[0]), np.log10(bins[-1]), len(bins))

        if label in other_colors.keys() and not isinstance(cmap, dict):
            color = other_colors[label]
        else:
            i += 1
            if not isinstance(cmap, dict):
                color = to_hex(cmap(i/ (len(data_file_paths)-len(other_colors))))
            else:
                color = cmap[label]
        bin_heights, bin_borders, _ = ax.hist(kcat_values, bins = logbins, histtype='step',
                                               fill=False, label= label, color=color, cumulative=cumulative)
        calculate_distribution_statistics(bin_heights, bin_borders)

    ax.vlines([13.7], 0, 1e4, linestyles='dotted')
    ax.tick_params(axis='both', labelsize=fontsize)

    plt.yscale('log')
    plt.ylabel('Frequency', fontsize=fontsize)
    plt.xscale('log')
    plt.xlabel('Kcat value [s-1]', fontsize=fontsize)
    if legend:
        plt.legend()
    return ax

def create_kcat_change_per_cog_barplot(original_pam_kcat_file: str,
                                       model_file: str,
                                       diagnostic_files: List[str],
                                       ax: plt.Axes,
                                       legend=True,
                                       fontsize = 16,
                                       cmap = None):
    individual_info_with_cog = get_kcat_changes_per_cog(original_pam_kcat_file,
                                       model_file,
                                       diagnostic_files)
    return create_cog_barplot(individual_info_with_cog, ax, legend=legend, fontsize=fontsize, cmap=cmap)

def get_kcat_changes_per_cog(original_pam_kcat_file:str,
                             model_file:str,
                             diagnostic_files:List[str])->pd.DataFrame:
    #extract kcat at starting position of optimization
    rxn2kcat, _ = get_rxn2kcat_protein2gene_dict(original_pam_kcat_file, model_file)

    best_individual_df = merge_all_diagnostic_files(diagnostic_files)
    best_individual_df = calculate_change_in_kcat(best_individual_df,
                                                  rxn2kcat)
    best_individual_with_cog = add_pathway_annotation_to_kcat_changes(best_individual_df)

    return summarize_and_pivot_cog_info_df_to_long(best_individual_with_cog)


def merge_all_diagnostic_files(diagnostic_files: List[str]) -> pd.DataFrame:
    best_individual_df = None

    for i, file in enumerate(diagnostic_files):
        pam_param_results = pd.read_excel(file, sheet_name='Best_Individuals')
        pam_param_error = pd.read_excel(file, sheet_name='Final_Errors')
        pam_param_results = pam_param_results.merge(pam_param_error, on='run_id', how='left')
        pam_param_results['alternative'] = i + 1

        if best_individual_df is None:
            best_individual_df = pam_param_results
        else:
            best_individual_df = pd.concat([best_individual_df, pam_param_results])
    return best_individual_df

def calculate_change_in_kcat(best_individual_df: pd.DataFrame,
                             rxn2kcat: Dict[str,dict]):
    best_individual_df['rxn_id'] = best_individual_df['rxn_id'].apply(lambda row:
                                                                      _extract_reaction_id_from_catalytic_reaction_id(
                                                                          row
                                                                      ))
    return best_individual_df.groupby(
        'alternative', group_keys=True).apply(calculate_kcat_differences, rxn2kcat).reset_index(drop=True)

def add_pathway_annotation_to_kcat_changes(best_individual_df: pd.DataFrame,
                              pathway_annotation_file:str = os.path.join('Data', 'GeneList_ecoli.xlsx'),
                              gene_to_protein_sheet:str = 'GeneList', rxn_to_gene_sheet: str = 'gene2rxn'
                                           )->pd.DataFrame:
    # get pathway information from genelist file
    gene_list_info = pd.read_excel(pathway_annotation_file, sheet_name=gene_to_protein_sheet)
    rxn2gene_info = pd.read_excel(pathway_annotation_file, sheet_name=rxn_to_gene_sheet)

    # map cog information to reactions using the gene id
    gene2rxn2cog = pd.merge(gene_list_info[['bnumber', 'COG description']], rxn2gene_info,
                            left_on='bnumber', right_on='Gene', how='inner')
    gene2rxn2cog.Reactions = gene2rxn2cog.Reactions.str.split(', ')
    gene2rxn2cog['COG description'] = gene2rxn2cog['COG description'].str.split(';')
    gene2rxn2cog = gene2rxn2cog.explode('Reactions')
    gene2rxn2cog = gene2rxn2cog.explode('COG description').reset_index(drop=True)

    best_kcats_with_cog = pd.merge(best_individual_df, gene2rxn2cog[['Reactions', 'COG description']],
                                   left_on='rxn_id', right_on='Reactions', how='left')
    best_kcats_with_cog = best_kcats_with_cog.drop_duplicates(
        ['rxn_id', 'enzyme_id', 'r_squared', 'COG description']).reset_index(drop=True)
    best_kcats_with_cog['enz_kcat_change']  = best_kcats_with_cog['kcat_change']
    total_change_per_iteration = best_kcats_with_cog.groupby(['alternative'])['kcat_change'].sum().abs().reset_index()
    return pd.merge(best_kcats_with_cog.drop('kcat_change', axis=1), total_change_per_iteration,
                             on=['alternative'])

def summarize_and_pivot_cog_info_df_to_long(cog_info_relative: pd.DataFrame):
    # Step 1: Compute the sums of positive and negative changes
    cog_info_relative['positive_change'] = cog_info_relative['enz_kcat_change'].apply(lambda x: x if x > 0 else 0)
    cog_info_relative['negative_change'] = cog_info_relative['enz_kcat_change'].apply(lambda x: x if x < 0 else 0)
    # Sum of changes grouped by alternative and COG description
    cog_summary = cog_info_relative.groupby(['COG description', 'alternative']).agg({
        'positive_change': 'sum',
        'negative_change': 'sum'
    }).reset_index()

    # Reshape the data for plotting
    cog_summary_long = pd.melt(
        cog_summary,
        id_vars=['COG description', 'alternative'],
        value_vars=['positive_change', 'negative_change'],
        var_name='Change Type',
        value_name='Change'
    )
    return cog_summary_long

def create_cog_barplot(cog_summary_long:pd.DataFrame,
                       ax:plt.Axes,
                       plotting_threshold=5e4,
                       bar_width=0.2, spacing_factor = 3,  # Increase spacing
                       fontsize = 15,
                       legend = True,
                       xlabel=r'$\sum \frac{k_{cat,new}-k_{cat,old}}{k_{cat,old}}$',
                       other_colors:dict={},
                       cmap:dict = None):
    # only plot the most important cogs
    cog_summary_long = cog_summary_long.groupby('COG description').filter(
        lambda x: x['Change'].abs().sum() > plotting_threshold
    )

    # Sort COG descriptions by total absolute change to get a nice descending plot
    cog_summary_long['abs_change'] = cog_summary_long['Change'].abs()  # Compute absolute changes
    # cog_summary_long = cog_summary_long.sort_values('abs_change').iloc[:num_pathways_to_plot,:]
    cog_order = sorted(
        set(cog_summary_long['COG description']),
        key=lambda x: abs(cog_summary_long[cog_summary_long['COG description'] == x]['Change']).sum(),
        reverse=False
    )

    #create the plot
    # fig, ax = plt.subplots(figsize=(12, len(cog_summary_long['COG description'].unique()) * 0.5), ax)

    # Define sorted y positions
    y_positions = np.arange(0, len(cog_order) * spacing_factor, spacing_factor)

    # Get unique alternatives and define a color mapping using coolwarm
    if cmap is None:
        alternatives = cog_summary_long['alternative'].unique()
        model_colors = sns.color_palette("viridis", n_colors=len(alternatives)-len(other_colors))
        cmap = {
            **{l: c for l, c in zip([f'Alternative {i + 1}' for i in range(len(alternatives)-len(other_colors))], model_colors)},
            **other_colors}

        print(cmap)
        cmap['Alternative 10'] = model_colors[0]


    # Plot each alternative separately
    for i, alt in enumerate(cog_summary_long['alternative'].unique()[::-1]):
        alt_data = cog_summary_long[cog_summary_long['alternative'] == alt]

        # Get y positions for existing data points
        y_pos_alt = [y_positions[cog_order.index(cog)] for cog in alt_data['COG description'] if cog in cog_order]

        # Separate positive and negative changes
        positive_data = alt_data[alt_data['Change'] > 0].reset_index()
        negative_data = alt_data[alt_data['Change'] < 0].reset_index()

        # Apply dodge effect by shifting bars horizontally
        # Ensure the number of y-positions matches the data length
        y_pos_positive = np.array(y_pos_alt[:len(positive_data)]) + (i - len(cog_summary_long['alternative'].unique()) / 2) * bar_width
        y_pos_negative = np.array(y_pos_alt[:len(negative_data)]) + (i - len(cog_summary_long['alternative'].unique()) / 2) * bar_width

        # Plot positive changes
        if isinstance(alt, str): label =alt
        else: label = f'Alternative {alt}'
        ax.barh(y_pos_positive, positive_data['Change'].values,
                color=cmap[label], label=label,
                height=bar_width, align='center')

        # Plot negative changes
        ax.barh(y_pos_negative, negative_data['Change'].values,
                color=cmap[label], height=bar_width, align='center')

    # Add a vertical reference line at 0
    ax.axvline(0, color='black', linestyle='--', linewidth=1)

    # make x-axis logaritmic
    # ax.set_xscale('symlog')
    # ax.set_xlim([-2.5 * 1e9, 2.5 * 1e9])

    # Adjust y-axis labels
    ax.tick_params(axis='both', labelsize=fontsize)
    ax.set_yticks(y_positions)
    ax.set_yticklabels([COG_MAPPER[cog] for cog in cog_order],
                       fontsize=fontsize)  # Set labels in sorted order

    # Label axes
    ax.set_xlabel(xlabel, fontsize=fontsize * 1.5)
    # ax.set_ylabel('COG Description', fontsize=fontsize * 1.5)
    ax.grid(visible=True, alpha=0.2, linewidth=0.7)
    ax.set_axisbelow(True)

    # Add legend (ensuring all alternatives are included)
    if legend:
        handles, labels = ax.get_legend_handles_labels()
        unique_labels = dict(zip(labels, handles))
        ax.legend(unique_labels.values(), unique_labels.keys(),
                  bbox_to_anchor=(1.05, 1), loc='upper left', fontsize = fontsize)

    # Show plot
    # plt.tight_layout()
    return ax

def get_simulated_protein_concentrations_for_alternatives(got_enzymes_file: str,
                                                          multi_file: str,
                                                          diagnostic_files: List[str],
                                                          model_file_path: str
                                                          ):
    substrate_rates = [0.67, 0.5, 0.35]
    pam_ge = set_up_pam(got_enzymes_file, model_file_path, sensitivity=False)
    pam = set_up_pam(multi_file, model_file_path, sensitivity=False)

    for model, label in zip(
        [pam_ge, pam], ['GotEnzymes', 'After preprocessing']
    ):
        enzymes = [enz.id for enz in pam_ge.enzyme_variables if enz._model is not None]

        proteomics_results = {label: get_results_from_simulations_fixed_mu(pamodel=model,
                                                                                  growth_rates=substrate_rates,
                                                                                  proteins_to_save=enzymes,
                                                                                  method_ids=['Batch', 'mu_5', 'mu_35']
                                                                                  )['proteins']}


    for file, label in zip(diagnostic_files,
                           [f'Alternative {i}' for i in range(1, len(diagnostic_files) + 1)]
                           ):
        print('\n')
        print(label)
        model = create_pamodel_from_diagnostics_file(file,
                                                     pam.copy(copy_with_pickle=True))
        enzymes = [enz.id for enz in model.enzyme_variables if enz._model is not None]
        proteomics_results[label] = get_results_from_simulations_fixed_mu(pamodel=model,
                                                                          growth_rates=substrate_rates,
                                                                          proteins_to_save=enzymes,
                                                                          method_ids=['Batch', 'mu_5', 'mu_35']
                                                                          )['proteins']
    enzyme_db = pd.read_excel(got_enzymes_file, sheet_name = 'ActiveEnzymes')
    masses_per_enzyme = enzyme_db[['enzyme_id', 'molMass']].drop_duplicates()

    prot_predicted_df_list = []
    for model_name, protein_df in proteomics_results.items():
        print(f"Processing {model_name}...")
        protein_df = pd.merge(protein_df, masses_per_enzyme, on='enzyme_id')
        # this function also corrects the concentration units from mmol to g
        protein_df = normalize_simulated_protein_concentrations(protein_df,
                                                                enzyme_db,
                                                                pam_ge.sectors.get_by_id('UnusedEnzymeSector'))
        protein_df['model'] = model_name
        prot_predicted_df_list.append(protein_df)

    # Merge all at once
    prot_predicted_long = pd.concat(prot_predicted_df_list, ignore_index=True)
    return prot_predicted_long

def parse_input_proteomics_data(ref_proteome_data_file:str,
                                uniprot_ref_file:str,
                                got_enzymes_file:str
                                ):
    locustag_regex = r'\b([b|s]\d{4})\b'
    # average values from Bakken and Olsen (1984); https://pubmed.ncbi.nlm.nih.gov/16346263/
    buoyant_volume = 1.09  # g_wetcell/cm^3
    dryweight_per_wetweight = 0.3
    DRYWEIGHT_PER_LITER = buoyant_volume * 1e3 * dryweight_per_wetweight  # 1e3 to convert cm^-3 to dm^-3
    proteome_df_metadata = pd.read_excel(ref_proteome_data_file,
                                         sheet_name='GrowthRates',
                                         engine='openpyxl',
                                         index_col=0)

    enzyme_db = pd.read_excel(got_enzymes_file, sheet_name = 'ActiveEnzymes')
    masses_per_enzyme = enzyme_db[['enzyme_id', 'molMass']].drop_duplicates()

    ref_proteomics_data = pd.read_excel(ref_proteome_data_file,
                                        sheet_name='ProteinMasses',
                                        engine='openpyxl',
                                        index_col=0)[['Glucose', 'Chemostat µ=0.5', 'Chemostat µ=0.35']]

    ref_proteomics_data.columns = ['Batch', 'mu_5', 'mu_35']

    # use information from uniprot to map the locus tag ids to protein identifiers
    uniprot_info_df = pd.read_excel(UNIPROT_INFO_FILE)

    # get the gene id from the gene names
    uniprot_info_df['b_number'] = uniprot_info_df['Gene Names'].str.extract(locustag_regex)
    print(uniprot_info_df)
    uniprot_df = uniprot_info_df[['b_number', 'Entry', 'Mass']]
    ref_proteomics_data_mapped = pd.merge(ref_proteomics_data, uniprot_df, how = 'right',
                               left_on = 'Bnumber', right_on='b_number')
    #make sure units are comparable
    proteomics_per_enzyme = convert_peptide_to_enzyme_concentrations(
        ref_proteomics_data_mapped.rename({'Entry': 'enzyme_id'}, axis=1),
        enzyme_db,
        concentration_columns=list(ref_proteomics_data.columns)
    ).dropna()
    proteomics_per_enzyme = (proteomics_per_enzyme.
                             merge(masses_per_enzyme, on='enzyme_id').
                             merge(uniprot_df[['Entry','b_number']], left_on='enzyme_id', right_on='Entry')
                             )
    proteomics_per_enzyme_mmol = proteomics_per_enzyme.copy()
    proteomics_per_enzyme_normalized = proteomics_per_enzyme.copy()
    for exp, col in zip(['Glucose', 'Chemostat µ=0.5', 'Chemostat µ=0.35'], ref_proteomics_data.columns):
        fl_per_cell = proteome_df_metadata.Cell_volume.loc[exp]
        # units are in fg/cell, cell * fl/cell * gcdw/l -> fgcdw
        proteomics_per_enzyme_mmol[col] = (proteomics_per_enzyme[col] /
                                  (fl_per_cell * DRYWEIGHT_PER_LITER))  # cell * fl/cell * gcdw/l -> fgcdw

        # need to sum all measured proteins while ignoring double entries because of isozymes
        total_protein_content = proteomics_per_enzyme_mmol.drop_duplicates('rxn_id')[col].sum()
        proteomics_per_enzyme_normalized[col] = proteomics_per_enzyme_mmol[col].div(total_protein_content)

    ref_proteomics_normalized = pd.melt(proteomics_per_enzyme_normalized, value_vars=['Batch', 'mu_5', 'mu_35'],
                                id_vars=['enzyme_id', 'rxn_id', 'b_number'],
                                var_name='experiment', value_name='fraction').rename({'experiment': 'method'}, axis = 1)
    return ref_proteomics_normalized

def get_changes_in_protein_abundance_per_cog(protein_conc_vs_cog_df):
    measured_proteins = protein_conc_vs_cog_df.drop_duplicates(['enzyme_id', 'method'])
    measured_proteins['model'] = 'Measurements'
    measured_proteins['fraction'] = measured_proteins['normalized_fraction']
    total_proteome_df_cog_sum = pd.concat([protein_conc_vs_cog_df, measured_proteins])
    total_proteome_df_cog_sum = total_proteome_df_cog_sum.drop_duplicates()

    # Step 1: Compute the sums of positive and negative changes
    total_proteome_df_cog_sum['positive_change'] = total_proteome_df_cog_sum['fraction'].apply(
        lambda x: x if x > 0 else 0)
    total_proteome_df_cog_sum['negative_change'] = total_proteome_df_cog_sum['fraction'].apply(
        lambda x: x if x < 0 else 0)

    # Compute the sum of 'fraction' for each model
    total_proteome_df_cog_sum['fraction_sums'] = total_proteome_df_cog_sum.groupby(['model'])['fraction'].transform(
        'sum')
    total_proteome_df_cog_sum['positive_change'] = total_proteome_df_cog_sum['positive_change'] / \
                                                   total_proteome_df_cog_sum['fraction_sums']
    # Aggregate the data and divide by the fraction sum per model

    cog_summary_sum = total_proteome_df_cog_sum.groupby(['COG Name', 'model']).agg({
        'positive_change': 'sum',
        'negative_change': 'sum'
    }).reset_index()

    cog_summary_long_sum = pd.melt(
        cog_summary_sum,
        id_vars=['COG Name', 'model'],
        value_vars=['positive_change', 'negative_change'],
        var_name='Change Type',
        value_name='Change'
    ).rename({'COG Name': 'COG description', 'model': 'alternative'}, axis=1)

    return cog_summary_long_sum


def create_proteomics_bargraph(got_enzymes_file: str,
                               multi_file: str,
                               diagnostic_files: List[str],
                               model_file_path: str,
                               ref_proteome_data_file: str,
                               uniprot_ref_file:str,
                               ax:plt.Axes,
                               ):

    protein2cog_df = pd.read_excel(ref_proteome_data_file,
                                    sheet_name='Gene2COG',
                                    engine='openpyxl',
                                    index_col=0)
    # simulated_protein_concentrations = get_simulated_protein_concentrations_for_alternatives(got_enzymes_file=got_enzymes_file,
    #                                                                            multi_file=multi_file,
    #                                                                            diagnostic_files=diagnostic_files,
    #                                                                            model_file_path=model_file_path)
    # simulated_protein_concentrations.to_excel('Results/3_analysis/predicted_protein_concentrations_normalized_iML1515.xlsx')
    simulated_protein_concentrations = pd.read_excel('Results/3_analysis/predicted_protein_concentrations_normalized_iML1515.xlsx')
    measured_protein_concentrations = parse_input_proteomics_data(ref_proteome_data_file=ref_proteome_data_file,
                                                                  uniprot_ref_file=uniprot_ref_file,
                                                                  got_enzymes_file=got_enzymes_file)
    all_protein_concentrations = pd.merge(simulated_protein_concentrations, measured_protein_concentrations,
                                 on=['method', 'enzyme_id', 'rxn_id'], how='inner')

    protein_conc_vs_cog_df = all_protein_concentrations.merge(protein2cog_df[['COG Name']],
                                                        left_on='b_number',
                                                        right_index=True).drop_duplicates()
    cog_changed_protein_conc = get_changes_in_protein_abundance_per_cog(protein_conc_vs_cog_df,)

    create_cog_barplot(cog_changed_protein_conc, ax, plotting_threshold=0.2, legend=True,
                       xlabel=r'$\sum\frac{g_{prot,COG}}{g_{protein}}$',
                       other_colors={'GotEnzymes': 'grey', 'Measurements': 'grey', 'After preprocessing': 'black'})

def create_proteomics_linegraphs(got_enzymes_file: str,
                               multi_file: str,
                               diagnostic_files: List[str],
                               model_file_path: str,
                               ref_proteome_data_file: str,
                               uniprot_ref_file:str,
                               ax:plt.Axes,
                               ):

    protein2cog_df = pd.read_excel(ref_proteome_data_file,
                                    sheet_name='Gene2COG',
                                    engine='openpyxl',
                                    index_col=0)
    # simulated_protein_concentrations = get_simulated_protein_concentrations_for_alternatives(got_enzymes_file=got_enzymes_file,
    #                                                                            multi_file=multi_file,
    #                                                                            diagnostic_files=diagnostic_files,
    #                                                                            model_file_path=model_file_path)
    # simulated_protein_concentrations.to_excel('Results/3_analysis/predicted_protein_concentrations_normalized_iML1515.xlsx')
    simulated_protein_concentrations = pd.read_excel('Results/3_analysis/predicted_protein_concentrations_normalized_iML1515.xlsx')
    measured_protein_concentrations = parse_input_proteomics_data(ref_proteome_data_file=ref_proteome_data_file,
                                                                  uniprot_ref_file=uniprot_ref_file,
                                                                  got_enzymes_file=got_enzymes_file)
    all_protein_concentrations = pd.merge(simulated_protein_concentrations, measured_protein_concentrations,
                                 on=['method', 'enzyme_id', 'rxn_id'], how='inner')

    protein_conc_vs_cog_df = (all_protein_concentrations
                              .merge(protein2cog_df[['COG Name']],
                                                        left_on='b_number',
                                                        right_index=True)
                              .drop_duplicates()[['enzyme_id', 'rxn_id', 'model', 'COG Name','method', 'normalized_fraction', 'fraction']]
                              .groupby(['model', 'COG Name','method'])
                              .sum()
                              .reset_index()
                              )
    top_cog = (protein_conc_vs_cog_df.groupby('COG Name')
               .sum()
               .sort_values(by=['normalized_fraction'], ascending=False)
               .head(5)
               .reset_index()
               )['COG Name']

    model_colors = sns.color_palette("viridis", n_colors=10)
    cmap = {
        **{l: c for l, c in
           zip([f'Alternative {i + 1}' for i in range(10)], model_colors)},
        **{'After preprocessing': 'black', 'measurements': 'black'}}

    print(cmap)
    method2mu = {'Batch': 0.67, 'mu_5': 0.5, 'mu_35':0.35}
    protein_conc_vs_cog_df['growth_rate'] = protein_conc_vs_cog_df['method'].apply(
        lambda x: method2mu[x])
    protein_conc_vs_cog_df = protein_conc_vs_cog_df.sort_values(by=['growth_rate'])
    fig, ax = plt.subplots(nrows=len(top_cog),figsize = (10,15))
    cog2index = {cog:i for i, cog in enumerate(top_cog)}
    for (model, cog), df in protein_conc_vs_cog_df.groupby(['model', 'COG Name']):
        if not cog in top_cog.to_list():continue
        i = cog2index[cog]
        ax[i].set_title(cog)
        ax[i].plot(df.growth_rate, df.normalized_fraction, label = model, color =cmap[model])
        ax[i].scatter(df.growth_rate, df.fraction, label = 'measurements', color =cmap['measurements'])
        ax[i].set_ylabel('Normalized Fraction [gProtein/gCDW]')
    h, l = ax[0].get_legend_handles_labels()
    handles, labels = [], []
    for handle, label in zip(h,l):
        if label in labels: continue
        handles.append(handle)
        labels.append(label)

    fig.legend(handles,labels, loc='upper right')
    fig.tight_layout()
    fig.subplots_adjust(right=0.8)
    fig.savefig('Figures/tests.png')

def main():
    NUM_ALT_MODELS = 10
    FONTSIZE = 14
    REF_PROTEOMICS_FILE = os.path.join('Data', 'Ecoli_phenotypes', 'proteome_data_extract_schmidt2016.xlsx')
    UNIPROT_INFO_FILE = os.path.join('Data', 'Databases', 'uniprotkb_ecolik12_240726.xlsx')
    PARAM_FILE_ORI = os.path.join('Results', '1_preprocessing',
                                  'proteinAllocationModel_iML1515_EnzymaticData_250912.xlsx')
    PARAM_FILE_PREPROC = os.path.join('Results', '2_parametrization',
                                      'proteinAllocationModel_iML1515_EnzymaticData_multi.xlsx')

    MODEL_FILE = os.path.join('Models', 'iML1515.xml')

    diagnostic_files = [os.path.join('Results', '2_parametrization', 'diagnostics',
                                     f'pam_parametrizer_diagnostics_{file_nmbr}.xlsx') for file_nmbr in
                        range(1, NUM_ALT_MODELS + 1)]
    parameter_files = [os.path.join('Results', '3_analysis', 'parameter_files',
                                    f'proteinAllocationModel_EnzymaticData_iML1515_{model}.xlsx') for model in
                       range(1, NUM_ALT_MODELS + 1)]

    model_colors = sns.color_palette("viridis", n_colors=NUM_ALT_MODELS)
    other_colors = {'GotEnzymes': 'grey', 'After preprocessing': 'black', 'iML1515': 'purple'}
    cmap = {
        **{l: c for l, c in
           zip([f'Alternative {i + 1}' for i in range(NUM_ALT_MODELS)], model_colors)},
        **other_colors}

    #create a pretty figure
    fig = plt.figure(figsize=(21/2.58, 15/2.58))
    plt.rcParams.update({'font.size': FONTSIZE})

    # Outer GridSpec: 2 rows (80% heatmaps, 20% colorbar)
    gs_main = gridspec.GridSpec(1, 2,wspace=0.6,
                                width_ratios=[10,6])
    gs_inner_l = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs_main[1],
                                                hspace=0.2,height_ratios=[5,3]
                                                )

    gs_inner_r = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs_main[0],
                                                    hspace=0.5,height_ratios=[10,1])

    protein_bar_ax = fig.add_subplot(gs_inner_l[0])
    create_proteomics_bargraph(got_enzymes_file=PARAM_FILE_ORI,
                               multi_file=PARAM_FILE_PREPROC,
                               diagnostic_files=diagnostic_files,
                               model_file_path=MODEL_FILE,
                               ref_proteome_data_file=REF_PROTEOMICS_FILE,
                               uniprot_ref_file=UNIPROT_INFO_FILE,
                               ax=protein_bar_ax)

    ax2 = fig.add_subplot(gs_inner_l[1])
    hist_ax = create_kcat_histogram_old_vs_new([PARAM_FILE_ORI,
                                                      PARAM_FILE_PREPROC] + parameter_files,
                                                     ax2,
                                                     label_names=['GotEnzymes', 'After preprocessing'] \
                                                                 + [f'Alternative {i}' for i in
                                                                    range(1, NUM_ALT_MODELS + 1)],
                                                     legend=False, fontsize=FONTSIZE, cmap = cmap)
    hist_ax.grid(visible=True, alpha=0.2, linewidth=0.7)
    hist_ax.set_axisbelow(True)

    # create a legend
    legend_ax = fig.add_subplot(gs_inner_r[1])
    handles, labels = [], []
    # for ax in [line_axs[0], ax2]:
    legend_ax.axis("off")  # Hide axes
    hist_ax.plot([],[],label='GotEnzymes', color='grey', linewidth=5)#dummy line for complete legend
    hist_ax.plot([],[],label='After preprocessing', color='black', linewidth=5)#dummy line for complete legend

    hist_ax.plot([],[],label='iML1515', color='black', linestyle ='--',linewidth=5)#dummy line for complete legend

    for ax in [hist_ax, ax2]:
        h, l = ax.get_legend_handles_labels()
        for label, handle in zip(l, h):
            if label not in labels:
                handles.extend([handle])
                labels.extend([label])
    legend_ax.legend(handles, labels, loc="center",
                     fontsize=FONTSIZE,
                     bbox_to_anchor = (0.2,0),
                     ncol=round(len(labels)/5), frameon=False)
    # ax2.legend(handles, labels, loc="lower center", ncol=round(len(labels)/2), frameon=False)

    ax3 = fig.add_subplot(gs_inner_r[0])
    bar_ax = create_kcat_change_per_cog_barplot(PARAM_FILE_PREPROC,
                                                MODEL_FILE,
                                                diagnostic_files,
                                                ax3, legend=False,
                                                cmap = cmap)

    #Add alphabet annotations
    #first 4 axes are added in reverse, and must ignore additional axis for x/y label config
    annotations = ["D", "C", "B", "A","","E","", "F"]
    fontsize = FONTSIZE  # Adjust as needed

    for ax, label in zip(fig.axes, annotations):
        ax.annotate(label, xy=(0, 1), xycoords="axes fraction",
                    fontsize=fontsize, fontweight='bold',
                    xytext=(-5, 5), textcoords="offset points",
                    ha="right", va="bottom")
    fig.tight_layout()
    fig.savefig(os.path.join('Figures', 'Figure3_parametrization_results_analysis.png'))

if __name__ == '__main__':
    REF_PROTEOMICS_FILE = os.path.join('Data', 'Ecoli_phenotypes', 'proteome_data_extract_schmidt2016.xlsx')
    UNIPROT_INFO_FILE = os.path.join('Data', 'Databases', 'uniprotkb_ecolik12_240726.xlsx')
    PARAM_FILE_ORI = os.path.join('Results', '1_preprocessing',
                                  'proteinAllocationModel_iML1515_EnzymaticData_250912.xlsx')
    PARAM_FILE_PREPROC = os.path.join('Results', '2_parametrization',
                                      'proteinAllocationModel_iML1515_EnzymaticData_multi.xlsx')

    MODEL_FILE = os.path.join('Models', 'iML1515.xml')
    NUM_ALT_MODELS = 10

    diagnostic_files = [os.path.join('Results', '2_parametrization', 'diagnostics',
                                     f'pam_parametrizer_diagnostics_{file_nmbr}.xlsx') for file_nmbr in
                        range(1, NUM_ALT_MODELS + 1)]
    create_proteomics_linegraphs(got_enzymes_file=PARAM_FILE_ORI,
                               multi_file=PARAM_FILE_PREPROC,
                               diagnostic_files=diagnostic_files,
                               model_file_path=MODEL_FILE,
                               ref_proteome_data_file=REF_PROTEOMICS_FILE,
                               uniprot_ref_file=UNIPROT_INFO_FILE,
                                 ax = ''
                               )
    # main()