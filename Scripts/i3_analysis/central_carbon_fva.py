from cobra.flux_analysis import flux_variability_analysis
import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

from typing import Dict, List

from cobra.io.sbml import read_sbml_model
from PAModelpy.utils import set_up_pam
from Modules.utils.pam_generation import setup_pputida_pam


def run_and_save_fva(models_to_check: Dict[str, 'Model'],
                     fva_result_file: str
                     ) -> None:
    for name, model in models_to_check.items():
        result = flux_variability_analysis(model = model)
        if not os.path.exists(fva_result_file): kwargs = {'mode':'w'}
        else:kwargs = {'mode':'a', 'if_sheet_exists':'replace', 'engine':'openpyxl'}
        with pd.ExcelWriter(fva_result_file, **kwargs) as writer:
            result.to_excel(writer, sheet_name=name)

def plot_fva_bargraph(fva_result: pd.DataFrame,
                      rxns_to_plot: List[str],
                      fig_result_file: str,
                      other_colors: Dict[str, str] = {'GotEnzymes': 'grey', 'iML1515': 'purple'},
                      legend = True
                      ) -> None:
    alternatives = fva_result['model'].unique()
    model_colors = sns.color_palette("viridis", n_colors=len(alternatives) - len(other_colors))
    cmap = {
            **{l: c for l, c in
               zip([f'Alternative_{i + 1}' for i in range(len(alternatives) - len(other_colors))], model_colors)},
            **other_colors}
    spacing_factor = 0.2
    bar_width = spacing_factor/(len(alternatives)+2)
    fontsize = 11
    xlabel = 'Flux [mmol/gCDW/h]'
    y_positions = np.arange(0, len(rxns_to_plot) * spacing_factor, spacing_factor)
    fig, ax = plt.subplots(figsize=(20/2.56, 20/2.56))

    # Plot each alternative separately
    for i, (model_id, model_df) in enumerate(fva_result.groupby('model')):
        model_df = (
            model_df.set_index('rxn')
            .reindex(rxns_to_plot)
            .reset_index()
        )

        offset = (list(alternatives).index(model_id) - (len(alternatives) - 1) / 2) * bar_width
        y_pos_alt = np.array(y_positions) + offset
        ax.barh(y_pos_alt, width = model_df['maximum'].values - model_df['minimum'].values,
                left = model_df['minimum'].values,
                color=cmap[model_id], label=model_id,
                height=bar_width, align='center')
        ax.scatter((model_df.maximum +model_df.minimum)/2,y_pos_alt,  marker ='|', color = cmap[model_id])
    # Add a vertical reference line at 0
    ax.axvline(0, color='black', linestyle='--', linewidth=1)

    # Adjust y-axis labels
    ax.tick_params(axis='both', labelsize=fontsize)
    ax.set_yticks(y_positions)
    ax.set_yticklabels(rxns_to_plot,
                       fontsize=fontsize)  # Set labels in sorted order

    # Label axes
    ax.set_xlabel(xlabel, fontsize=fontsize * 1.5)
    ax.set_xlim([-25,25])
    # ax.set_ylabel('COG Description', fontsize=fontsize * 1.5)
    ax.grid()

    # Add legend (ensuring all alternatives are included)
    if legend:
        handles, labels = ax.get_legend_handles_labels()
        unique_labels = dict(zip(labels, handles))
        ax.legend(unique_labels.values(), unique_labels.keys(),
                  bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=fontsize)
    plt.subplots_adjust(right = 0.7)
    plt.savefig(fig_result_file)

def main_ecoli():
    NUM_MODELS = 10
    FVA_RESULT_FILE = os.path.join('Results', '3_analysis', 'fva_central_carbon_metabolism_ecoli.xlsx')
    FIG_RESULT_FILE = os.path.join('Results', '3_analysis','fva_central_carbon_metabolism_iML1515.png')
    rxns_to_plot = ['GLCptspp', 'HEX1', 'PGI', 'PFK', 'F6PA', 'TALA', 'FBA', 'GAPD']
    models_to_check = {f'Alternative_{i}': os.path.join('Results', '3_analysis', 'parameter_files',
                                    f'proteinAllocationModel_EnzymaticData_iML1515_{i}.xlsx') for i in range(1, NUM_MODELS+1)}
    # models_to_check = {}
    models_to_check['GotEnzymes'] = (os.path.join('Results','2_parametrization', 'proteinAllocationModel_iML1515_EnzymaticData_multi.xlsx'))
    models_to_check = {name: set_up_pam(file, sensitivity=False) for name, file in models_to_check.items()}
    models_to_check['iML1515'] = read_sbml_model(os.path.join('Models', 'iML1515.xml'))

    run_and_save_fva(models_to_check = models_to_check, fva_result_file=FVA_RESULT_FILE)
    fva_result = pd.concat(
        [df.reset_index(names='rxn').assign(model=name) for name, df in
         pd.read_excel(FVA_RESULT_FILE, sheet_name=None, index_col=0).items()],
        ignore_index=True
    )
    fva_result = fva_result[fva_result['rxn'].isin(rxns_to_plot)]
    print(fva_result.to_markdown())
    plot_fva_bargraph(fva_result=fva_result, rxns_to_plot=rxns_to_plot, fig_result_file = FIG_RESULT_FILE)

def main_pputida():
    NUM_MODELS = 5
    GLC_UPTAKE_FLUXOMICS = 6.14
    FVA_RESULT_FILE = os.path.join('Results', '3_analysis', 'fva_central_carbon_metabolism_pputida.xlsx')
    FIG_RESULT_FILE = os.path.join('Results', '3_analysis','fva_central_carbon_metabolism_iJN1463.png')
    rxns_to_plot = ['GLCptspp', 'HEX1', 'PGI', 'PFK', 'F6PA', 'TALA', 'FBA', 'GAPD', 'PDH', 'CS', 'OAADC', 'ACONTa', 'ACONTb', 'ICDHyr', 'SUCOAS', 'SUCDi', 'FUM', 'MDH', 'ME2']
    models_to_check = {f'Alternative_{i}': os.path.join('Results', '3_analysis', 'parameter_files',
                                                        f'proteinAllocationModel_EnzymaticData_iJN1463_{i}.xlsx') for i
                       in range(1, NUM_MODELS + 1)}
    # models_to_check = {}
    models_to_check['GotEnzymes'] = (
        os.path.join('Results', '2_parametrization', 'proteinAllocationModel_iJN1463_EnzymaticData_multi.xlsx'))
    models_to_check = {
        name: setup_pputida_pam(file, sensitivity=False)
                       for name, file in models_to_check.items()
    }
    for model in models_to_check.values():
        model.change_reaction_bounds('EX_glc__D_e', -GLC_UPTAKE_FLUXOMICS, 0)

    models_to_check['iJN1463'] = read_sbml_model(os.path.join('Models', 'iJN1463.xml'))
    models_to_check['iJN1463'].reactions.EX_glc__D_e.lower_bound = -GLC_UPTAKE_FLUXOMICS

    run_and_save_fva(models_to_check=models_to_check, fva_result_file=FVA_RESULT_FILE)
    fva_result = pd.concat(
        [df.reset_index(names='rxn').assign(model=name) for name, df in
         pd.read_excel(FVA_RESULT_FILE, sheet_name=None, index_col=0).items()],
        ignore_index=True
    )
    fva_result = fva_result[fva_result['rxn'].isin(rxns_to_plot)]
    print(fva_result.to_markdown())
    plot_fva_bargraph(fva_result=fva_result,
                      rxns_to_plot=rxns_to_plot,
                      fig_result_file=FIG_RESULT_FILE,
                      other_colors = {'GotEnzymes': 'grey', 'iJN1463': 'purple'})



if __name__ == '__main__':
    main_pputida()


