import matplotlib.pyplot as plt
from matplotlib.colors import to_hex
from matplotlib import gridspec
import seaborn as sns
import numpy as np
import pandas as pd
import os
from typing import List, Dict, Tuple, Callable
from PAModelpy.utils import set_up_pam

from Scripts.i2_parametrization.pam_parametrizer_iML1515 import set_up_pamparametrizer

from Modules.utils.pam_generation import (get_rxn2kcat_protein2gene_dict,
                                          _extract_reaction_id_from_catalytic_reaction_id,
                                          _get_rxn2kcat_as_series,
                                          create_pamodel_from_diagnostics_file
                                         )
from Modules.utils.pamparametrizer_analysis import (calculate_kcat_differences,
                                                    plot_histogram_logspace,
                                                    select_clustered_rows_by_variation,
                                                    get_clusters_from_clustermap,
                                                    set_up_pam_parametrizer_and_get_substrate_uptake_rates
                                                   )

from Modules.utils.pamparametrizer_visualization import plot_valid_data, plot_simulation, plot_flux_vs_experiment


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
                                     fontsize =16):
    # fig, ax = plt.subplots()
    n_bins = 50
    i = 0
    cmap = plt.get_cmap("coolwarm")

    for label, data_file_path in zip(label_names, data_file_paths):
        aes_parameter_df = pd.read_excel(data_file_path, sheet_name='ActiveEnzymes')
        kcat_values = aes_parameter_df.kcat_values.dropna()
        print('------------------------------------------------------------------------')
        print(f'The kcat set from {label} has:\n \tMedian:\t\t\t{np.median(kcat_values)} \n \tMean:\t\t\t{np.mean(kcat_values)}')


        hist, bins = np.histogram(kcat_values, bins=n_bins)
        logbins = np.logspace(np.log10(bins[0]), np.log10(bins[-1]), len(bins))

        if label in other_colors.keys():
            color = other_colors[label]
        else:
            i += 1
            color = to_hex(cmap(i/ (len(data_file_paths)-len(other_colors))))
        bin_heights, bin_borders, _ = ax.hist(kcat_values, bins = logbins, histtype='step',
                                              stacked=True, fill=False, label= label, color=color, cumulative=cumulative)
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
                                       fontsize = 16):
    individual_info_with_cog = get_kcat_changes_per_cog(original_pam_kcat_file,
                                       model_file,
                                       diagnostic_files)
    return create_cog_barplot(individual_info_with_cog, ax, legend=legend, fontsize=fontsize)

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
    best_kcats_with_cog['relative_change']  = best_kcats_with_cog['kcat_change']
    total_change_per_iteration = best_kcats_with_cog.groupby(['alternative'])['kcat_change'].sum().abs().reset_index()
    return pd.merge(best_kcats_with_cog.drop('kcat_change', axis=1), total_change_per_iteration,
                             on=['alternative'])

def summarize_and_pivot_cog_info_df_to_long(cog_info_relative: pd.DataFrame):
    # Step 1: Compute the sums of positive and negative changes
    cog_info_relative['positive_change'] = cog_info_relative['kcat_change'].apply(lambda x: x if x > 0 else 0)
    cog_info_relative['negative_change'] = cog_info_relative['kcat_change'].apply(lambda x: x if x < 0 else 0)

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
                       plotting_threshold=5*1e11,
                       bar_width=0.2, spacing_factor = 3,  # Increase spacing
                       fontsize = 15,
                       legend = True,
                       xlabel=r'$\sum \frac{k_{cat,new}-k_{cat,old}}{k_{cat,old}}$',
                       other_colors:dict={}):
    # only plot the most important cogs
    cog_summary_long = cog_summary_long.groupby('COG description').filter(
        lambda x: x['Change'].abs().sum() > plotting_threshold
    )

    # Sort COG descriptions by total absolute change to get a nice descending plot
    cog_summary_long['abs_change'] = cog_summary_long['Change'].abs()  # Compute absolute changes
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
    alternatives = cog_summary_long['alternative'].unique()
    model_colors = sns.color_palette("coolwarm", n_colors=len(alternatives)-len(other_colors))
    alternative_colors = {
        **{l: c for l, c in zip([f'Alternative {i + 1}' for i in range(len(alternatives)-len(other_colors))], model_colors)},
        **other_colors}


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
        y_pos_positive = np.array(y_pos_alt[:len(positive_data)]) + (i - len(alternatives) / 2) * bar_width
        y_pos_negative = np.array(y_pos_alt[:len(negative_data)]) + (i - len(alternatives) / 2) * bar_width

        # Plot positive changes
        if isinstance(alt, str): label =alt
        else: label = f'Alternative {alt}'
        ax.barh(y_pos_positive, positive_data['Change'].values,
                color=alternative_colors[alt], label=label,
                height=bar_width, align='center')

        # Plot negative changes
        ax.barh(y_pos_negative, negative_data['Change'].values,
                color=alternative_colors[alt], height=bar_width, align='center')

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

    # Add legend (ensuring all alternatives are included)
    if legend:
        handles, labels = ax.get_legend_handles_labels()
        unique_labels = dict(zip(labels, handles))
        ax.legend(unique_labels.values(), unique_labels.keys(),
                  bbox_to_anchor=(1.05, 1), loc='upper left', fontsize = fontsize)

    # Show plot
    # plt.tight_layout()
    return ax

def recreate_progress_plot(best_indiv_files:list[str],
                           labels:list[str], fig,
                           axs,
                           legend:bool = True,
                           fontsize=20,
                           pamparam_setup: Callable= None,
                           pamparam_kwargs: dict = {'max_substrate_uptake_rate':-0.1},
                           rxns_to_plot: List[str] = None,
                           substrate_uptake_id:str = 'EX_glc__D_e',
                           other_measurements: bool = False
                           ):
    j=0

    if pamparam_setup is None:
        parametrizer, substrate_rates = set_up_ecoli_pam_parametrizer_and_get_substrate_uptake_rates()
    else:
        parametrizer, substrate_rates = set_up_pam_parametrizer_and_get_substrate_uptake_rates(pamparam_setup,
                                                                                               pamparam_kwargs)
    if rxns_to_plot is not None:
        parametrizer.validation_data.get_by_id(substrate_uptake_id)._reactions_to_plot = rxns_to_plot

    substrate_rates = sorted(substrate_rates)
    fig, axs = plot_valid_data(parametrizer,axs, fig, fontsize=fontsize)
    print('Run reference simulations')
    fluxes, _ = parametrizer.run_simulations_to_plot(substrate_uptake_id='EX_glc__D_e',
                                                                   substrate_rates = substrate_rates,
                                                                   sensitivity = False)
    fig, axs = plot_simulation(fig, axs, fluxes, [abs(rate) for rate in substrate_rates],
                               parametrizer.validation_data.get_by_id('EX_glc__D_e')._reactions_to_plot,
                               iteration=0, color='black')

    for file, label in zip(best_indiv_files, labels):
        j +=1
        print('\nAlternative ', label, ' from file ', file)
        parametrizer.pamodel = create_pamodel_from_diagnostics_file(file,
                                                                    parametrizer._pamodel.copy(copy_with_pickle = True))
        fluxes, _ = parametrizer.run_simulations_to_plot(substrate_uptake_id='EX_glc__D_e',
                                                                       substrate_rates=substrate_rates,
                                                                       sensitivity=False)
        fig, axs, color = plot_simulation(fig, axs, fluxes, [abs(rate) for rate in substrate_rates],
                                   parametrizer.validation_data.get_by_id('EX_glc__D_e')._reactions_to_plot,
                                   iteration=j + 1, max_iteration=len(best_indiv_files), label = label,
                                   return_color=True)
        if other_measurements:
            print('plotting other carbon sources')
            plot_flux_vs_experiment(axs[len(rxns_to_plot)], parametrizer,
                                    color)


    lines, labels = fig.axes[1].get_legend_handles_labels()

    if legend:
        fig.legend(lines, labels, loc='upper left', bbox_to_anchor=(0.1, 0.99), frameon=False,
                   fontsize=fontsize - 5)
    # fig.tight_layout()
    return fig, axs

def set_up_ecoli_pam_parametrizer_and_get_substrate_uptake_rates() -> Tuple:
    kwargs = {'min_substrate_uptake_rate':-12,
              'max_substrate_uptake_rate': -0.1,
              'kcat_increase_factor': 3}
    return set_up_pam_parametrizer_and_get_substrate_uptake_rates(set_up_pamparametrizer,
                                                           kwargs)

def main():
    NUM_ALT_MODELS = 8
    FONTSIZE = 16
    PARAM_FILE_ORI = os.path.join('Results', '1_preprocessing',
                                  'proteinAllocationModel_iML1515_EnzymaticData_250225.xlsx')
    PARAM_FILE_PREPROC = os.path.join('Results', '2_parametrization',
                                      'proteinAllocationModel_iML1515_EnzymaticData_multi.xlsx')

    MODEL_FILE = os.path.join('Models', 'iML1515.xml')

    diagnostic_files = [os.path.join('Results', '2_parametrization', 'diagnostics',
                                     f'pam_parametrizer_diagnostics_{file_nmbr}.xlsx') for file_nmbr in
                        range(1, NUM_ALT_MODELS + 1)]
    parameter_files = [os.path.join('Results', '3_analysis', 'parameter_files',
                                    f'proteinAllocationModel_EnzymaticData_iML1515_{model}.xlsx') for model in
                       range(1, NUM_ALT_MODELS + 1)]

    #create a pretty figure
    fig = plt.figure(figsize=(30, 30))

    # Outer GridSpec: 2 rows (80% heatmaps, 20% colorbar)
    gs_main = gridspec.GridSpec(1, 2,wspace=0.6,
                                width_ratios=[10,6])
    gs_inner_l = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs_main[0],
                                                hspace=0.2,height_ratios=[5,3]
                                                )

    gs_inner_top = gridspec.GridSpecFromSubplotSpec(2, 2, subplot_spec=gs_inner_l[0,0],
                                                    wspace = 0.4,
                                                    hspace=0.4)
    gs_inner_r = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs_main[1],
                                                    hspace=0.5,height_ratios=[10,1])
    axs1 =[None,None,None,None]
    for i in reversed(range(4)):  # Assuming line_axs[0] is a row of axes
        row, col = (i, 0)  # Determine row/column position in the 2x2 grid
        kwargs={}
        if i>1:
            row, col = (i-2, 1)
        if i==0 or i==1:
            kwargs={'sharex':axs1[i+1]}
        axs1[i]= fig.add_subplot(gs_inner_top[row, col], **kwargs)
    line, line_axs = recreate_progress_plot(diagnostic_files,
                                           labels=[f'Alternative {i}' for i in range(1, NUM_ALT_MODELS + 1)],
                                            fig=fig, axs=axs1, legend=False, fontsize=FONTSIZE)
    #share x and y axis labels for progress plot
    ax_group = fig.add_subplot(gs_inner_l[0])
    ax_group.set_xticks([])
    ax_group.set_yticks([])
    ax_group.set_frame_on(False)
    ax_group.set_ylabel(r"Flux [mmol/$\text{g}_{\text{CDW}}$/h]",
                        labelpad=20, fontsize=FONTSIZE)
    ax_group.yaxis.set_label_coords(-0.1, 0.5)  # Move ylabel to the right
    ax_group.set_xlabel(r"Glucose uptake [mmol/$\text{g}_{\text{CDW}}$/h]",
                        labelpad=20, fontsize=FONTSIZE)

    ax2 = fig.add_subplot(gs_inner_l[1])
    hist_ax = create_kcat_histogram_old_vs_new([PARAM_FILE_ORI,
                                                      PARAM_FILE_PREPROC] + parameter_files,
                                                     ax2,
                                                     label_names=['GotEnzymes', 'After preprocessing'] \
                                                                 + [f'Alternative {i}' for i in
                                                                    range(1, NUM_ALT_MODELS + 1)],
                                                     legend=False, fontsize=FONTSIZE)

    # create a legend
    legend_ax = fig.add_subplot(gs_inner_r[1])
    handles, labels = [], []
    # for ax in [line_axs[0], ax2]:
    legend_ax.axis("off")  # Hide axes
    line_axs[0].plot([],[],label='GotEnzymes', color='grey')#dummy line for complete legend
    for ax in [line_axs[0], ax2]:
        h, l = ax.get_legend_handles_labels()
        for label, handle in zip(l, h):
            if label not in labels:
                handles.extend([handle])
                labels.extend([label])
    legend_ax.legend(handles, labels, loc="center",
                     fontsize=FONTSIZE,
                     ncol=round(len(labels)/4), frameon=False)
    # ax2.legend(handles, labels, loc="lower center", ncol=round(len(labels)/2), frameon=False)

    ax3 = fig.add_subplot(gs_inner_r[0])
    bar_ax = create_kcat_change_per_cog_barplot(PARAM_FILE_PREPROC,
                                                MODEL_FILE,
                                                diagnostic_files,
                                                ax3, legend=False)
    # Get the current position of ax3
    # pos = ax3.get_position()

    # Shrink the axis to fit inside the GridSpec boundaries
    # ax3.set_position([pos.x0 + 0.2, pos.y0 + 0.2, pos.width - 0.2, pos.height - 0.2])

    # # create a legend
    # legend_ax = fig.add_subplot(gs_main[2])
    # legend_ax.axis("off")  # Hide axes
    #
    # # Manually create a legend based on the entries of the histogram
    # handles, labels = [], []
    # for ax in [line_axs[0], ax2]:
    #     h, l = ax.get_legend_handles_labels()
    #     print(h,l)
    #     if l not in labels:
    #         handles.extend(h)
    #         labels.extend(l)
    # legend_ax.legend(handles, labels, loc="center", ncol=len(labels), frameon=False)

    #Add alphabet annotations
    #first 4 axes are added in reverse, and must ignore additional axis for x/y label config
    annotations = ["D", "C", "B", "A","","E","",  "F"]
    fontsize = FONTSIZE  # Adjust as needed

    for ax, label in zip(fig.axes, annotations):
        ax.annotate(label, xy=(0, 1), xycoords="axes fraction",
                    fontsize=fontsize, fontweight='bold',
                    xytext=(-5, 5), textcoords="offset points",
                    ha="right", va="bottom")
    # fig.set_constrained_layout(True)
    fig.tight_layout()
    # plt.show()
    # fig.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1, hspace=0.3, wspace=0.2)
    fig.savefig(os.path.join('Figures', 'Figure1_parametrization_results_analysis.png'))

if __name__ == '__main__':
    main()