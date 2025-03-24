import os
from matplotlib import gridspec
import matplotlib.pyplot as plt

from Figures.Scripts.Figure1_iml1515_kcat_analysis import recreate_progress_plot
from Scripts.i2_parametrization.pam_parametrizer_iCGB21FR import set_up_pamparametrizer as pamparam_setup_icgb21fr
from Scripts.i2_parametrization.pam_parametrizer_iJN1463 import set_up_pamparametrizer as pamparam_setup_ijn1463


NUM_MODELS = 5
PAM_KCAT_FILES_ICG = [os.path.join('Results', '2_parametrization', 'diagnostics',
                               f'pam_parametrizer_diagnostics_iCGB21FR_{file_nmbr}.xlsx') for file_nmbr in
                  range(1, NUM_MODELS + 1)]

PAM_KCAT_FILES_IJN = [os.path.join('Results', '2_parametrization', 'diagnostics',
                               f'pam_parametrizer_diagnostics_iJN1463_{file_nmbr}.xlsx') for file_nmbr in
                  range(1, NUM_MODELS + 1)]
FONTSIZE=16

labels = [f'Alternative {i}' for i in range(1, NUM_MODELS+1)]

def main():
    fig, axs = plt.subplots(2,4)
    i=0
    for pamparamsetup, kcat_file_list, kwargs, rxns_to_plot in zip([pamparam_setup_icgb21fr, pamparam_setup_ijn1463],
                                             [PAM_KCAT_FILES_ICG, PAM_KCAT_FILES_IJN],
                                             [{'max_substrate_uptake_rate':-0.1,
                                               'c_sources': ['Glucose', 'Fructose', 'Succinate','Gluconate', 'Acetate']
                                               },
                                              {'max_substrate_uptake_rate':-0.1,
                                               'min_substrate_uptake_rate':-15,
                                               'c_sources': ['Glycerol', 'Glucose','Octanoate',
                                                            'm-Xylene','Succinate', 'Benzoate',
                                                            'Fructose']
                                               }
                                              ],[['Growth', 'EX_co2_e', 'EX_o2_e'],['BIOMASS_KT2440_WT3']]
                                             ):

         recreate_progress_plot(kcat_file_list,
                                labels, fig, axs[i],
                                fontsize=FONTSIZE,
                                pamparam_setup=pamparamsetup,
                                pamparam_kwargs = kwargs,
                                rxns_to_plot = rxns_to_plot,
                                other_measurements = True)
         i+=1

    fig.savefig(os.path.join('Figures', 'Figure3_cglutanicum_pputida.png'))

if __name__ == '__main__':
    main()