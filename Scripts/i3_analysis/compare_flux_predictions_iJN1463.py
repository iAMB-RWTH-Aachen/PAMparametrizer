import os
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from cobra import DictList
from PAModelpy import PAModel

from Scripts.i2_parametrization.pam_parametrizer_iJN1463 import set_up_validation_data
from Scripts.pam_generation_uniprot_id import setup_pputida_pam
from Modules.utils.pam_generation import create_pamodel_from_diagnostics_file
from Modules.utils.pamparametrizer_analysis import get_results_from_simulations

def perform_simulation_experiments(model: PAModel,
                                   validation_data_list: DictList,
                                   translational_sector: dict[str, dict[str,float]]) -> pd.DataFrame:
    substrate_ids = []
    substrate_rates= []
    fluxes_to_save = set()
    for vd in validation_data_list:
        substrate_ids += [vd.id]
        fluxes_to_save.update(vd.valid_data.columns)
        substrate_rates.append(list(vd.valid_data[vd.id+'_ub']))

    return get_results_from_simulations(model,
                                 substrate_ids = substrate_ids,
                                 substrate_rates = substrate_rates,
                                 fluxes_to_save=list(fluxes_to_save),
                                        transl_sector_config = translational_sector)['fluxes']

def plot_growth_rate_vs_experiments(simulation_results: pd.DataFrame,
                                    validation_data: DictList,
                                    axs: plt.Axes = None,
                                    model: str = 'GotEnzymes',
                                    fontsize:int = 16):
    if axs is None:
        fig, axs = plt.subplot()

    axs.plot([0,1], [0,1], linestyle = 'dashed')
    model_colors = sns.color_palette("coolwarm", n_colors=len(validation_data))
    cmap = dict(zip([vd.id for vd in validation_data], model_colors))

    for vd in validation_data:
        results = simulation_results[simulation_results.substrate_id == vd.id].sort_values('substrate')
        if len(results) == 0:
            continue
        axs.scatter(vd.valid_data.sort_values(vd.id).BIOMASS_KT2440_WT3, results.BIOMASS_KT2440_WT3,
                    color = cmap[vd.id],
                    label = vd.id)

    axs.set_title(f'Alternative model {model}', fontsize=fontsize)
    axs.set_xlabel('simulated growth rate (1/h)', fontsize=fontsize)
    # axs.set_ylabel('experimental growth rate (1/h)')
    return axs


def plot_simulations_vs_experiments(simulation_results: pd.DataFrame, validation_data: DictList,
                                    axs: plt.Axes = None, fontsize:int = 16):
    if axs is None:
        fig, axs = plt.subplot()

    axs.plot([0,15], [0,15], linestyle = 'dashed')

    model_colors = sns.color_palette("coolwarm", n_colors=len(validation_data))
    cmap = dict(zip([vd.id for vd in validation_data], model_colors))

    for vd in validation_data:
        results = simulation_results[simulation_results.substrate_id == vd.id].sort_values('substrate')
        if len(results) == 0:
            continue
        validation = vd.valid_data.sort_values([vd.id+'_ub'])
        for rxn in results.columns[2:]:
            if rxn in validation.columns and rxn in results.columns:
                axs.scatter(validation[rxn].abs(), results[rxn].abs(), color = cmap[vd.id], label=vd.id)

    axs.set_xlabel('simulated fluxes (mmol/gCDW/h)', fontsize = fontsize)
    # axs.set_ylabel('experimental fluxes (mmol/gCDW/h)')
    return axs


if __name__ == '__main__':
    RESULT_FIGURE_FILE_PATH = os.path.join('Results', '3_analysis', 'flux_comparison_iJN1463.png')
    NUM_MODELS = 2
    DIAGNOSTIC_FILES = [os.path.join(
        'Results', '2_parametrization', 'diagnostics', f'pam_parametrizer_diagnostics_iJN1463_{i}.xlsx')
        for i in range(1, NUM_MODELS+1)
    ]
    fig, axs = plt.subplots(ncols=2,nrows=NUM_MODELS+1,figsize = [7,10])

    fontsize = 15

    pam = setup_pputida_pam(sensitivity = False)
    pam.change_reaction_bounds('EX_glc__D_e', 0,1e3)
    validation_data = set_up_validation_data()
    for i, file in enumerate(DIAGNOSTIC_FILES):
        new_pam = create_pamodel_from_diagnostics_file(file, pam.copy(copy_with_pickle=True))
        translational_sector = pd.read_excel(file, sheet_name='translational_sector').set_index('substrate_uptake_id').T.to_dict()
        simulation_results = perform_simulation_experiments(pam, validation_data, translational_sector)
        plot_growth_rate_vs_experiments(simulation_results, validation_data,axs[i,0], model = i+1, fontsize = fontsize)
        plot_simulations_vs_experiments(simulation_results, validation_data, axs[i,1],fontsize = fontsize)


    pam = setup_pputida_pam(os.path.join(
        'Results', '1_preprocessing','proteinAllocationModel_iJN1463_EnzymaticData_250117.xlsx'),
        sensitivity = False
    )
    simulation_results = perform_simulation_experiments(pam, validation_data, translational_sector)

    plot_growth_rate_vs_experiments(simulation_results, validation_data, axs[i+1, 0])
    plot_simulations_vs_experiments(simulation_results, validation_data, axs[i+1, 1])

    # # Shrink current axis's height by 10% on the bottom
    box = axs[NUM_MODELS,1].get_position()
    axs[NUM_MODELS,1].set_position([box.x0, box.y0 + box.height * 0.1,
                            box.width, box.height * 0.9])
    handles, labels = axs[0,0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='lower center', bbox_to_anchor=(0.5, -0.05), ncols=8, fontsize=fontsize)
    fig.supylabel('experimental growth rate (1/h)', fontsize = fontsize)
    # plt.legend()
    plt.tight_layout()
    plt.savefig(RESULT_FIGURE_FILE_PATH)



