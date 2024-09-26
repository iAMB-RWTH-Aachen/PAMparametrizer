from matplotlib import pyplot as plt
import pandas as pd
import os
import sys
import numpy as np
if os.path.split(os.getcwd())[1] == 'Figures':
    os.chdir(os.path.split(os.getcwd())[0])
sys.path.append('C:\\Users\\claud\\Documents\\iamb-student-folders\\iamb-folder-template\\mcPAM_package')

from mcPAModelpy.configuration import Config
from Scripts.mcpam_generation_uniprot_id  import setup_ecolicore_pam, set_up_ecoli_pam

# print(mcpam_core.solver.shadow_prices)
config = Config()
config.reset()
def change_kcats_for_an_enyzme_set(model, enzyme_set):
    for enzyme, kcats in enzyme_set.items():
        model.change_kcat_value(enzyme, kcats)


mcpam_core = set_up_ecoli_pam(sensitivity=True)
mcpam_core.tolerance = 1e-9

enzyme_set1 = {'P06959_P0A9P0_P0AFG8': {'CE_PDH_P06959_P0A9P0_P0AFG8': {'f': 1e8, 'b': 0}},
              'P0A6P9': {'CE_ENO_P0A6P9': {'f': 1e8, 'b': 15.49830396049956}},
              'P0AB71': {'CE_FBA_P0AB71': {'f': 1e8, 'b': 362.7443342149361}},
              'P33570': {'CE_TKT1_P33570': {'f': 1e8, 'b': 0}, 'CE_TKT2_P33570': {'f': 0.5, 'b': 2.345}, 'TKT2': {'f': 0.5, 'b': 2.345}},
              'P52697': {'CE_PGL_P52697': {'f': 1e8, 'b': 221.220906006779}},
              'P69054_P0AC44_P0AC41_P07014': {'CE_SUCDi_P69054_P0AC44_P0AC41_P07014': {'f': 1e8, 'b': 0}}}


enzyme_set2 = {'P00864': {'CE_PPC_P00864': {'f': 1e8, 'b': 1.941478463907463}},
              'P06999': {'CE_PFK_P06999': {'f': 1e8, 'b': 0}},
              'P08200': {'CE_ICDHyr_P08200': {'f': 1e8, 'b': 31.73895571943552}},
              'P0A799': {'CE_PGK_P0A799': {'f': 1e8, 'b': 23.3638}},
              'P25516': {'CE_ACONTa_P25516': {'f': 1e8, 'b': 81.96157312716353},
                      'CE_ACONTb_P25516': {'f': 1e8, 'b': 1.255995794640362},
                      'ACONTb': {'f': 1e8, 'b': 1.255995794640362}},
              'P33570': {'CE_TKT1_P33570': {'f': 1e8, 'b': 0},
                         'CE_TKT2_P33570': {'f': 1e8, 'b': 2.345}, 'TKT2': {'f': 2.345, 'b': 2.345}},
              'P36683': {'CE_ACONTa_P36683': {'f': 1e8, 'b': 308105.6127807793},
                         'CE_ACONTb_P36683': {'f': 1e8, 'b': 4.4574}},
              'P61889': {'CE_MDH_P61889': {'f': 1e8, 'b': 5.6591}}}

enzyme_set3 = {'P00864': {'CE_PPC_P00864': {'f': 1e8, 'b': 1.941478463907463}},
              'P06999': {'CE_PFK_P06999': {'f': 1e8, 'b': 0}},
              'P0A799': {'CE_PGK_P0A799': {'f': 1e8, 'b': 23.3638}},
              'P25516': {'CE_ACONTa_P25516': {'f': 1e8, 'b': 81.96157312716353},
                      'CE_ACONTb_P25516': {'f': 1e8, 'b': 1.255995794640362},
                      'ACONTb': {'f': 1e8, 'b': 1.255995794640362}},
              'P36683': {'CE_ACONTa_P36683': {'f': 1e8, 'b': 308105.6127807793},
                         'CE_ACONTb_P36683': {'f': 1e8, 'b': 4.4574}},
              'P61889': {'CE_MDH_P61889': {'f': 1e8, 'b': 5.6591}}}

enzyme_set4 = {'P00864': {'CE_PPC_P00864': {'f': 1e8, 'b': 1.941478463907463}},
              'P25516': {'CE_ACONTa_P25516': {'f': 1e8, 'b': 81.96157312716353},
                      'CE_ACONTb_P25516': {'f': 1e8, 'b': 1.255995794640362},
                      'ACONTb': {'f': 1e8, 'b': 1.255995794640362}}}

# Changing the kcats
change_kcats_for_an_enyzme_set(mcpam_core, enzyme_set1)
change_kcats_for_an_enyzme_set(mcpam_core, enzyme_set4)
#
# for enz_complex in mcpam_core.enzyme_variables:
#     print(enz_complex, enz_complex.kcats)

# BIOMASS_REACTION = 'BIOMASS_Ecoli_core_w_GAM'
#
# mcpam_core.objective = BIOMASS_REACTION
mcpam_core.optimize()

print('mcpam objective value:', mcpam_core.objective.value)

occupied_area, available_area = mcpam_core.calculate_occupied_membrane()
print('membrane constraint', mcpam_core.constraints['membrane'])
print(f'occupied area {occupied_area/available_area*100}%')

#### Fluxes simulations for different glc uptake rates

# load phenotype data from excel file
pt_data = pd.read_excel(os.path.join('Data', 'Ecoli_phenotypes','Ecoli_phenotypes_py_rev.xls'),
                        sheet_name='Yields', index_col=None)

# extract reaction specific data
rxn_to_pt = {}
rxn_transform = {
    'EX_ac_e': 'EX_ac_e',
    'EX_co2_e': 'EX_co2_e',
    'EX_o2_e': 'EX_o2_e',
    'BIOMASS_Ecoli_core_w_GAM':'BIOMASS_Ec_iML1515_core_75p37M'
    # 'BIOMASS_Ecoli_core_w_GAM':'BIOMASS_Ecoli_core_w_GAM'
}
for rxn_id, pt_id in rxn_transform.items():
    rxn_to_pt[rxn_id] = pt_data[['EX_glc__D_e', pt_id]].dropna().rename(columns={pt_id: rxn_id})


# Optimizing mcpam core
with mcpam_core:
    # change glucose uptake rate
    mcpam_core.reactions.EX_glc__D_e.lower_bound = -6.0
    # solve the model
    sol_mcpam = mcpam_core.optimize()
    # print(pamodel.summary())
    # with pd.option_context('display.max_rows', None):
    #     print(sol_pam.fluxes)
mcpam_core.optimize()

glc_uptake_rates = np.linspace(0.5, 11.5, 20)

# Initializing fluxes and concentrations for pam core
fluxes_pam_core = []
concentrations_pam_core = [0]

# Initializing fluxes and concentrations for mcpam core
fluxes_mcpam_core = []
concentrations_mcpam_core = [0]

# Simulating mcpam core
for glc in glc_uptake_rates:
    with mcpam_core:
        # change glucose uptake rate
        mcpam_core.reactions.EX_glc__D_e.lower_bound = -glc
        # disable pyruvate formate lyase (inhibited by oxygen)
        mcpam_core.reactions.PFL.upper_bound = 0
        # solve the model
        sol_mcpam = mcpam_core.optimize()
        # save data
        fluxes_mcpam_core.append(sol_mcpam.fluxes)  # flux distributions
        concentration = 0
        for enz_var in mcpam_core.enzyme_variables:
            concentration += enz_var.concentration
        concentrations_mcpam_core.append(concentration)

# plot flux changes with glucose uptake
rxn_id = ['EX_ac_e', 'EX_co2_e', 'EX_o2_e', 'BIOMASS_Ecoli_core_w_GAM']
fig, axs = plt.subplots(2, 2, dpi=90)
for r, ax in zip(rxn_id, axs.flatten()):
    # plot data
    if r in rxn_to_pt.keys():
        ax.scatter(abs(rxn_to_pt[r]['EX_glc__D_e']), abs(rxn_to_pt[r][r]),
                   color='firebrick', marker='o', s=30, linewidths=1.3,
                   facecolors=None, zorder=0,
                   label='Data')

    # plot simulation for mcpam core
    ax.plot(glc_uptake_rates, [abs(f[r]) for f in fluxes_mcpam_core],
            label='Simulation mcpam core', linewidth=2.5,
            zorder=5)

    # options
    ax.set_xlabel('glc uptake rate [mmol/gDW/h]')
    ax.set_ylabel('flux [mmol/gDW/h]')
    ax.set_title(r)
    # set grid
    ax.grid(True, axis='both', linestyle='--', linewidth=0.5, alpha=0.6)
    ax.set_axisbelow(True)
    # show legend
    ax.legend(fontsize=8, edgecolor='white', facecolor='white', framealpha=1)

# Add parameter box
param_text = (
    f"mcPAM/PAM_core model \n"
    f"Total protein: {mcpam_core.total_protein_fraction} g/g DW \n"
    f"Max membrane area = {mcpam_core.membrane_sector.max_membrane_area * 100}% \n"
    f"mcPAM occupied area = {occupied_area/available_area*100}%"
)
# Position the text box in the upper left corner of each subplot
fig.text(0.6, 2.75, param_text, transform=ax.transAxes, fontsize=8,
        verticalalignment='top', bbox=dict(facecolor='white', alpha=0.5, boxstyle="round,pad=0.3"))

plt.tight_layout()
plt.show()
