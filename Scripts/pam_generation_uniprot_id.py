import cobra
import pandas as pd
import numpy as np
import os
from typing import Union
import json


# load PAMpy modules
from PAModelpy.PAModel import PAModel
from PAModelpy.EnzymeSectors import ActiveEnzymeSector, UnusedEnzymeSector, TransEnzymeSector, CustomSector
from PAModelpy.configuration import Config
from PAModelpy.utils import set_up_pam


'Function library for making Protein Allocation Models as described in the publication'

def setup_ecolicore_pam(total_protein:bool = True,
                         active_enzymes: bool = True,
                         translational_enzymes:bool = True,
                         unused_enzymes:bool = True,
                         sensitivity:bool =True):
    # Setting the relative paths
    PAM_DATA_FILE_PATH = os.path.join('Data', 'proteinAllocationModel_iML1515_EnzymaticData_240730.xlsx')

    config = Config()
    config.reset()

    # some other constants
    BIOMASS_REACTION = 'BIOMASS_Ecoli_core_w_GAM'
    TOTAL_PROTEIN_CONCENTRATION = 0.16995  # [g_prot/g_cdw]
    config.BIOMASS_REACTION = BIOMASS_REACTION

    # load the genome-scale information
    model = cobra.io.load_json_model(os.path.join('Models', 'e_coli_core.json'))

    #load example data for the E.coli iML1515 model
    if active_enzymes:
        # load active enzyme sector information
        enzyme_db = parse_enzyme_db(pd.read_excel(PAM_DATA_FILE_PATH, sheet_name='ActiveEnzymes'))
        # create enzyme objects for each gene-associated reaction
        rxn2protein, protein2gene = parse_reaction2protein(enzyme_db, model)

        active_enzyme_sector = ActiveEnzymeSector(rxn2protein=rxn2protein, protein2gene=protein2gene,
                                                  configuration=config)

    else:
        active_enzyme_sector = None

    if translational_enzymes:
        # translational protein sector parameter (substrate dependent)
        id_list_tps = ['EX_glc__D_e']
        tps_0 = [0.04992]  # g/gDW
        tps_mu = [-0.002944]  # g h/gDW -> transformed to match glucose uptake variable
        molmass_tps = [405903.94]  # g/mol

        # translational protein sector
        translation_enzyme_sector = TransEnzymeSector(
            id_list=id_list_tps,
            tps_0=tps_0,
            tps_mu=tps_mu,
            mol_mass=molmass_tps,
        )
    else:
        translation_enzyme_sector = None

    if unused_enzymes:
        id_list_ups = [BIOMASS_REACTION]
        ups_0 = [0.0407]  # g/gDW
        ups_mu = [-0.0214]  # g h/gDW -> negative relation with growth rate
        molmass_ups = [405903.94]  # g/mol

        unused_enzyme_sector = UnusedEnzymeSector(
            id_list=id_list_ups,
            ups_0=ups_0,
            ups_mu=ups_mu,
            mol_mass=molmass_ups,
        )
    else:
        unused_enzyme_sector = None

    if total_protein: total_protein = TOTAL_PROTEIN_CONCENTRATION

    pa_model = PAModel(id_or_model=model, p_tot=total_protein, sensitivity=sensitivity,
                       active_sector=active_enzyme_sector, translational_sector=translation_enzyme_sector,
                       unused_sector=unused_enzyme_sector, configuration=config)
    return pa_model

def set_up_ecoli_pam(pam_info_file:str= os.path.join('Data', 'proteinAllocationModel_iML1515_EnzymaticData_240730.xlsx'),
                     model:str = 'iML1515.xml', config:Config = None,
                     total_protein: Union[bool, float] = True, active_enzymes: bool = True,
                     translational_enzymes: bool = True, unused_enzymes: bool = True, sensitivity = True,
                     enzyme_db:pd.DataFrame = None):


    if config is None:
        config = Config()
        config.reset()

    # some other constants
    TOTAL_PROTEIN_CONCENTRATION = 0.258  # [g_prot/g_cdw]

    #setup the gem ecoli iML1515 model
    model = cobra.io.read_sbml_model(os.path.join('Models', model))

    #check if a different total protein concentration is given
    if isinstance(total_protein, float):
        TOTAL_PROTEIN_CONCENTRATION = total_protein

    # load example data for the E.coli iML1515 model
    if active_enzymes:
        # load active enzyme sector information
        if enzyme_db is None:
            enzyme_db = parse_enzyme_db(pd.read_excel(pam_info_file, sheet_name='ActiveEnzymes'))
        # create enzyme objects for each gene-associated reaction
        rxn2protein, protein2gene = parse_reaction2protein(enzyme_db, model)

        active_enzyme_sector = ActiveEnzymeSector(rxn2protein=rxn2protein, protein2gene=protein2gene,
                                                  configuration=config)

    else:
        active_enzyme_sector = None

    if translational_enzymes:
        translational_info = pd.read_excel(pam_info_file, sheet_name='Translational')
        translation_enzyme_sector = TransEnzymeSector(
            id_list=[translational_info[translational_info.Parameter == 'id_list'].loc[0, 'Value']],
            tps_0=[translational_info[translational_info.Parameter == 'tps_0'].loc[1, 'Value']],
            tps_mu=[translational_info[translational_info.Parameter == 'tps_mu'].loc[2, 'Value']],
            mol_mass=[translational_info[translational_info.Parameter == 'mol_mass'].loc[3, 'Value']],
            configuration = config)
    else:
        translation_enzyme_sector = None

    if unused_enzymes:
        unused_protein_info = pd.read_excel(pam_info_file, sheet_name='UnusedEnzyme').set_index('Parameter')

        ups_0 = unused_protein_info.at['ups_0', 'Value']
        ups_mu = unused_protein_info.at['ups_mu', 'Value']

        unused_protein_sector = UnusedEnzymeSector(
            id_list=[unused_protein_info.at['id_list', 'Value']],
            ups_mu=[ups_mu],
            ups_0=[ups_0],
            mol_mass=[unused_protein_info.at['mol_mass', 'Value']],
            configuration = config)
    else:
        unused_protein_sector = None

    if total_protein: total_protein = TOTAL_PROTEIN_CONCENTRATION

    pamodel = PAModel(id_or_model=model, p_tot=total_protein,
                       active_sector=active_enzyme_sector, translational_sector=translation_enzyme_sector,
                       unused_sector=unused_protein_sector, sensitivity=sensitivity, configuration = config
                      )
    return pamodel


def _parse_enzyme_information_from_file(file_path:str):
    # load active enzyme sector information
    enzyme_db = pd.read_excel(file_path, sheet_name='enzyme-gene-reaction')
    # enzyme_db = enzyme_db.set_index('rxnID')
    # correct reaction IDS
    for idx in enzyme_db.rxnID.to_list():
        # transport reactions
        if 'pp' in idx:
            idx_new = idx.replace('pp', '')
            if idx_new not in enzyme_db.index:
                enzyme_db.rename(index={idx: idx_new}, inplace=True)
        if 'ex' in idx:
            idx_new = idx.replace('ex', '')
            if idx_new not in enzyme_db.index:
                enzyme_db.rename(index={idx: idx_new}, inplace=True)

                # replace NaN values with unique identifiers
    # replace NaN enzyme ids with a dummy enzyme identifier
    # select the NaN values
    nan_values = enzyme_db['EC_nmbr'].isnull()
    # make a list with unique ids
    nan_ids = [f'E{i}' for i in range(nan_values.sum())]
    # replace nan values by unique id
    enzyme_db.loc[nan_values, 'EC_nmbr'] = nan_ids

    return enzyme_db

def parse_enzyme_db(enzyme_db: pd.DataFrame) -> None:
    try:
        enzyme_db = enzyme_db.drop(['kegg_id', 'Reactants', 'Products', 'EC', 'Length'],axis=1).pivot_table(
                          index=['rxn_id', 'GPR', 'gene', 'uniprot_id', 'molMass'],
                          columns='direction',
                          values='kcat_values').reset_index().rename({'f': 'kcat_f', 'b':'kcat_b'}, axis=1)
    except:
        enzyme_db = enzyme_db.pivot_table(
            index=['rxn_id', 'GPR', 'gene', 'uniprot_id', 'molMass'],
            columns='direction',
            values='kcat_values').reset_index().rename({'f': 'kcat_f', 'b': 'kcat_b'}, axis=1)
    return enzyme_db


def parse_reaction2protein(enzyme_db: pd.DataFrame, model:cobra.Model) -> dict:
    # Initialize dictionaries
    rxn2protein = {}
    protein2gpr = {}

    # replace NaN values with unique identifiers
    # select the NaN values
    nan_values = enzyme_db['uniprot_id'].isnull()
    # make a list with unique ids
    nan_ids = [f'E{i}' for i in range(nan_values.sum())]
    # replace nan values by unique id
    enzyme_db.loc[nan_values, 'uniprot_id'] = nan_ids


    protein2gene, gene2protein = _get_genes_for_proteins(enzyme_db, model)

    # Iterate over each row in the DataFrame
    for index, row in enzyme_db.iterrows():
        # Parse data from the row
        rxn_id = row['rxn_id']
        # only parse those reactions which are in the model
        if rxn_id not in model.reactions: continue
        kcat_f_b = [row['kcat_f'], row['kcat_b']]
        kcat_f_b = [kcat if not np.isnan(kcat) else 0 for kcat in kcat_f_b]
        # if all([np.isnan(kcat) for kcat in kcat_f_b]): continue

        rxn = model.reactions.get_by_id(rxn_id)
        # get the identifiers and replace nan values by dummy placeholders
        enzyme_id = row['uniprot_id']
        gene_id = row['gene']

        # check if there are genes associates with the reaction #TODO move these if statements
        if len(rxn.genes) > 0 or isinstance(gene_id, str): #TODO if not: stop the loop and continue
            if not isinstance(enzyme_id, str):
                enzyme_id = 'Enzyme_' + rxn_id
                row['molMass'] = 39959.4825  # default molmass
            if not isinstance(gene_id, str):
                gene_id = [gene.id for gene in rxn.genes][0]  # TODO

            # get the gene-protein-reaction-associations for this specific enzyme
            gr, pr = parse_gpr_information_for_rxn2protein(row['GPR'],
                                                           gene2protein, protein2gene,enzyme_id)

            if pr is None: pr = [[enzyme_id]]

            if enzyme_id not in protein2gpr:
                protein2gpr[enzyme_id] = gr
            else:
                protein2gpr[enzyme_id].append(gr)

            # Create rxn2protein dictionary
            if rxn_id not in rxn2protein:
                rxn2protein[rxn_id] = {}
            if enzyme_id not in rxn2protein[rxn_id]:
                rxn2protein[rxn_id][enzyme_id] = {
                    'f': kcat_f_b[0],  # Forward kcat
                    'b': kcat_f_b[1],  # Backward kcat
                    'molmass': row['molMass'],
                    'genes': gr,
                    'protein_reaction_association': pr
                }
            else:
                rxn2protein[rxn_id][enzyme_id]['genes'].append(gene_id)

    # if no enzyme info is found, add dummy enzyme with median kcat and molmass
    for rxn in model.reactions:
        if (
                rxn.id not in rxn2protein.keys()
                and 'EX'.lower() not in rxn.id.lower()
                and 'BIOMASS' not in rxn.id
                and len(rxn._genes) > 0
                and list(rxn._genes)[0].id != 's0001'
        ): #TODO make if not, continue
            rev = _check_reaction_reversibility(rxn)
            if rev == 0:
                kcat_dict = {'f': 22}
            elif rev == 1:
                kcat_dict = {'b': 22}
            else:
                kcat_dict = {'f': 22, 'b': 22} #TODO put in dict?
            # no enzyme information found
            print('No enzyme information found for reaction: ' + rxn.id)
            enzyme_id = 'Enzyme_' + rxn.id
            gpr_info = parse_gpr_information_for_protein2genes(rxn.gpr)

            rxn2protein[rxn.id] = {enzyme_id: {
                **kcat_dict,
                'molmass': 3.947778784340140e04,
                'protein_reaction_association': [[enzyme_id]],
                'genes': gpr_info
            }}
            # add geneinfo for unknown enzymes
            protein2gpr[enzyme_id] = gpr_info

    return rxn2protein, protein2gpr

def _get_genes_for_proteins(enzyme_db: pd.DataFrame, model) -> dict:
    protein2gene = {}
    gene2protein = {}
    for index, row in enzyme_db.iterrows():
        # Parse data from the row
        rxn_id = row['rxn_id']
        if rxn_id not in model.reactions: continue
        rxn = model.reactions.get_by_id(rxn_id)
        # get the identifiers and replace nan values by dummy placeholders
        enzyme_id = row['uniprot_id']
        gene_id = row['gene']

        # check if there are genes associates with the reaction
        if len(rxn.genes) > 0 or isinstance(gene_id, str):
            if not isinstance(enzyme_id, str):
                enzyme_id = 'Enzyme_' + rxn_id
                row['molMass'] = 39959.4825  # default molmass
            if not isinstance(gene_id, str):
                gene_id = [gene.id for gene in rxn.genes][0]  # TODO

            gene2protein[gene_id] = enzyme_id

            # Create gene-protein-reaction associations
            if enzyme_id not in protein2gene:
                protein2gene[enzyme_id] = [gene_id]
            elif enzyme_id in protein2gene:
                # assume that there is a single protein if the previous relation was not assigned a gene
                if 'gene_' in protein2gene[enzyme_id][0]:
                    protein2gene[enzyme_id] = [gene_id]
                else:
                    protein2gene[enzyme_id].append(gene_id)

    return protein2gene, gene2protein


def _check_rxn_identifier_format(rxn_id:str) -> str:
    if 'pp' in rxn_id:
        idx_new = rxn_id.replace('pp', '')
    elif 'ex' in rxn_id:
        idx_new = rxn_id.replace('ex', '')
    else:
        idx_new = rxn_id
    return idx_new

def _get_fwd_bckw_kcat(rxn_id: str, kcat:float, model:PAModel) -> Union[list, None]:
    #skip exchange and biomass reaction
    if 'EX' in rxn_id or 'BIOMASS' in rxn_id:
        return None

    # Extract the base identifier without any suffixes
    base_id = rxn_id.split('_')[0]

    # Iterate over each identifier in the input
    if base_id in model.reactions:
        # Determine the form of the identifier
        if rxn_id.endswith('_f'):
            kcat_fwd = kcat
            kcat_rev = 0
        elif rxn_id.endswith('_b'):
            kcat_fwd = 0
            kcat_rev = kcat
        # identifier has no suffix
        elif base_id == rxn_id:
            kcat_fwd = kcat
            kcat_rev = kcat
        else:
            return None
    elif rxn_id in model.reactions:
        kcat_fwd = kcat
        kcat_rev = kcat
    else:
        return None
    return [kcat_fwd, kcat_rev]

def _check_reaction_reversibility(reaction):
    if reaction.lower_bound >= 0:
        # irreversible reaction (forward direction)
        rev = 0
    elif reaction.upper_bound <= 0:
        # irreversible reaction (reverse direction)
        rev = 1
    else:
        rev = 2
        # reversible r
    return rev

def parse_gpr_information_for_protein2genes(gpr_info:str):
    #filter out nan entries
    if isinstance(gpr_info, str):
        gpr_list = parse_gpr(gpr_info)
        # Extracting the inner lists and removing parentheses
        # gpr_list = [[[item.strip("(')")] + sublist[1:] for item in sublist[0].split(" or ")] for sublist in nested_list]
        return gpr_list
    else:
        return [['gene_dummy']]

def parse_gpr_information_for_rxn2protein(gpr_info:str, gene2protein: dict,protein2gene:dict, enzyme_id:str):
    #filter out nan entries
    if isinstance(gpr_info, str):
        gpr_list = parse_gpr(gpr_info)
        # Extracting the inner lists and removing parentheses
        gene = protein2gene[enzyme_id]

        if isinstance(protein2gene[enzyme_id], list):
            gene = gene[0]
        #only get the genes associated with this enzyme
        gpr_list = filter_sublists(gpr_list, gene)

        #convert the genes to the associated proteins
        enzyme_relations = []
        for sublist in gpr_list:
            enz_sublist = []
            for item in sublist:
                if item in gene2protein.keys():
                    enz_sublist.append(gene2protein[item])
            enzyme_relations += [enz_sublist]
        # only get the enzyme relations which involves this enzyme
        enzyme_relations = filter_sublists(enzyme_relations, enzyme_id)
        return gpr_list, enzyme_relations
    else:
        return [['gene_dummy']], None

def parse_gpr(gpr_info):
    # Split the string by 'or' first
    or_groups = gpr_info.split(' or ')
    parsed_gpr = []
    for group in or_groups:
        # Within each 'or' group, split by 'and'
        and_genes = group.split(' and ')
        # Clean up each gene string
        and_genes = [gene.strip(' ()') for gene in and_genes]
        parsed_gpr.append(and_genes)
    return parsed_gpr

def filter_sublists(nested_list, target_string):
    """
    Filters out all sublists from a nested list that contain the target string.

    Args:
        nested_list (list of list of str): The nested list to filter.
        target_string (str): The string to filter out sublists that contain it.

    Returns:
        list of list of str: A new nested list with the filtered sublists.
    """
    return [sublist for sublist in nested_list if target_string in sublist]

def setup_yeast_pam(pam_info_file:str= os.path.join('Data', 'proteinAllocationModel_yeast9_EnzymaticData_TurnUp.xlsx'),
                     model:str = 'Models/yeast9.xml', config:Config = None,
                     total_protein: Union[bool, float] = 0.28697423725932236, active_enzymes: bool = True,
                    translational_enzymes: bool = True, unused_enzymes: bool = True, sensitivity = True) -> PAModel:
    if config is None:
        config = set_up_yeast_config()
    yeast_pam = set_up_pam(pam_info_file, model, config,
                     total_protein, active_enzymes, translational_enzymes,
                     unused_enzymes, sensitivity = sensitivity)
    if total_protein:
        #the data we used to calibrate the model includes also housekeeping proteins in the UP section. Thus we
        #can assume include this section in the total protein constraint (no correction needed)
        protein_growth_relation = CustomSector(
            id_list= ['r_2111'], #biomass formation
            name='total_protein_growth_relation',
            cps_0=[0],
            cps_s=[-0.45431074007034344],
            mol_mass=1e6
        )
        yeast_pam.add_sector(protein_growth_relation)
    return yeast_pam


def set_up_yeast_config():
    config = Config()
    config.TOTAL_PROTEIN_CONSTRAINT_ID = "TotalProteinConstraint"
    config.P_TOT_DEFAULT = 0.388  # g_protein/g_cdw
    config.CO2_EXHANGE_RXNID = "r_1672"
    config.GLUCOSE_EXCHANGE_RXNID = "r_1714"
    config.BIOMASS_REACTION = "r_2111"
    config.OXYGEN_UPTAKE_RXNID = "r_1992"
    config.ACETATE_EXCRETION_RXNID = "r_1634"
    config.PHYS_RXN_IDS = [
    config.BIOMASS_REACTION,
    config.GLUCOSE_EXCHANGE_RXNID,
    config.ACETATE_EXCRETION_RXNID,
    config.CO2_EXHANGE_RXNID,
    config.OXYGEN_UPTAKE_RXNID]
    config.ENZYME_ID_REGEX = r'(Y[A-P][LR][0-9]{3}[CW])'
    return config

def setup_pputida_pam(pam_info_file:str= os.path.join(
                                             'Results', '1_preprocessing',
                                             'proteinAllocationModel_iJN1463_EnzymaticData_250110.xlsx'),
                     model:str = 'Models/iJN1463.xml',
                     total_protein: Union[bool, float] = 0.3, active_enzymes: bool = True,
                    translational_enzymes: bool = True, unused_enzymes: bool = True, sensitivity = True):
    config = Config()
    config.reset()
    config.BIOMASS_REACTION = 'BIOMASS_KT2440_WT3'
    pputida_pam = set_up_pam(pam_info_file, model, config,
                     total_protein, active_enzymes, translational_enzymes,
                     unused_enzymes, sensitivity = sensitivity)
    return pputida_pam

def get_rxn2kcat_protein2gene_dict(param_file):
    # create enzyme objects for each gene-associated reaction
    ecoli_pam = set_up_ecoli_pam(param_file, sensitivity=False)
    enzyme_db = pd.read_excel(param_file, sheet_name='ActiveEnzymes').iloc[:,1:]
    rxn2protein, protein2gene = parse_reaction2protein(enzyme_db, ecoli_pam)

    ae_sector = ecoli_pam.sectors.ActiveEnzymeSector
    new_rxn2prot = ae_sector.rxn2protein.copy()
    for rxn, enz_dict in ae_sector.rxn2protein.items():
        if rxn[:2]!='CE':
            for enzyme_id, enzyme_dict in enz_dict.items():
                protein_reaction = enzyme_dict['protein_reaction_association']
                if ae_sector._enzyme_is_enzyme_complex(protein_reaction, enzyme_id):
                    for pr in protein_reaction:
                        if len(pr) > 1:
                            enzyme_complex_id = '_'.join(pr)
                            new_rxn2prot[rxn] = {**new_rxn2prot[rxn],
                                                 **{enzyme_complex_id: enzyme_dict}}
    return new_rxn2prot, protein2gene

def pickle_rxn2kcat_protein2gene_dicts(pam_info_file: str = os.path.join('Data', 'proteinAllocationModel_iML1515_EnzymaticData_240730.xlsx')):
    # setup the gem ecoli iML1515 model
    model = cobra.io.read_sbml_model(os.path.join('Models', 'iML1515.xml'))
    enzyme_db = pd.read_excel(pam_info_file, sheet_name='ActiveEnzymes')

    # create enzyme objects for each gene-associated reaction
    rxn2kcat, protein2gene = get_rxn2kcat_protein2gene_dict(enzyme_db, model)

    # with open('Models/rxn2protein_iML1515_uniprot.json', 'w') as file:
    #     json.dump(rxn2protein, file, indent=4, sort_keys=True)
    with open('Models/rxn2kcat_iML1515_uniprot.json', 'w') as file:
        json.dump(rxn2kcat, file, indent=4, sort_keys=True)
    with open('Models/protein2gene_iML1515_uniprot.json', 'w') as file:
        json.dump(protein2gene, file, indent=4, sort_keys=True)


def increase_kcats_in_parameter_file(kcat_increase_factor: int,
                                     pam_info_file_path_ori: str,
                                     pam_info_file_path_out: str = os.path.join(
    'Data', 'proteinAllocationModel_iML1515_EnzymaticData_240730_multi.xlsx')):

    pam_info_file_ori = pd.read_excel(pam_info_file_path_ori, sheet_name='ActiveEnzymes')
    transprot = pd.read_excel(pam_info_file_path_ori, sheet_name='Translational')
    unusedenz = pd.read_excel(pam_info_file_path_ori, sheet_name='UnusedEnzyme')

    pam_info_file_new = pam_info_file_ori.rename({'kcat_values': 'kcat_values_ori'})
    pam_info_file_new['kcat_values'] = pam_info_file_ori['kcat_values']*kcat_increase_factor

    write_mode = 'w'
    kwargs = {}
    if os.path.isfile(pam_info_file_path_out):
        write_mode = 'a'
        kwargs = {'if_sheet_exists': 'replace'}

    with pd.ExcelWriter(pam_info_file_path_out, engine = 'openpyxl',
                        mode = write_mode, **kwargs) as writer:
        pam_info_file_new.to_excel(writer, sheet_name = 'ActiveEnzymes', index=False)
        transprot.to_excel(writer, sheet_name='Translational', index=False)
        unusedenz.to_excel(writer, sheet_name='UnusedEnzyme', index=False)




if __name__ == '__main__':
    # pam = setup_yeast_pam(os.path.join(
    #     'Results', '1_preprocessing',  'yeast9','proteinAllocationModel_yeast9_EnzymaticData_TurnUp_multi.xlsx'),
    #     sensitivity=True)
    # # pam.change_total_protein_constraint(1)
    # pam.change_reaction_bounds('r_1714', -1e3,0)
    # pam.optimize()
    # print(pam.summary())
    # print(pam.enzyme_sensitivity_coefficients.sort_values('coefficient', ascending=False))
    # print(pam.capacity_sensitivity_coefficients.sort_values('coefficient', ascending=False))

    pam_info_file: str = os.path.join('Results', '1_preprocessing', 'proteinAllocationModel_iML1515_EnzymaticData_241209.xlsx')
    model = cobra.io.read_sbml_model(os.path.join('Models', 'iML1515.xml'))

    pam = set_up_pam(pam_info_file, model)
    print(pam.catalytic_events.CE_ICYSDS.enzymes)
    # enzyme_db = pd.read_excel(pam_info_file, sheet_name='ActiveEnzymes')

    # create enzyme objects for each gene-associated reaction
    # rxn2kcat, protein2gene = get_rxn2kcat_protein2gene_dict(enzyme_db, model)


