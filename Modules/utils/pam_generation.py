import pandas as pd
import numpy as np
import cobra
from typing import TypedDict, Literal

from collections import defaultdict
from dataclasses import dataclass

DEFAULT_MOLMASS = 39959.4825 #kDa
DEFAULT_KCAT = 11 #s-1

#TODO: test functions,
class EnzymeInformation(TypedDict):
    enzyme_id:str
    f:float  # Forward kcat
    b:float  # Backward kcat
    genes: list[str]
    protein_reaction_association: str
    molmass: float

@dataclass
class ReactionInformation:
    rxnid: str
    enzymes: list[EnzymeInformation]
    model_reaction: cobra.Reaction = None

    def get_reaction_from_model(self, model):
        self.model_reaction = model.reactions.get_by_id(self.rxnid)

    def check_reaction_reversibility(self) -> Literal[-1,0,1]:
        if self.model_reaction.lower_bound >= 0:
            # irreversible reaction (forward direction)
            rev = 1
        elif self.model_reaction.upper_bound <= 0:
            # irreversible reaction (reverse direction)
            rev = -1
        else:
            rev = 0
            # reversible r
        return rev

def create_enzyme_information(rxn_id: str,
                              genes: list[str],
                              protein_reaction_association: list[str],
                              enzyme_id: str = None,
                              kcat_f: float= DEFAULT_KCAT,
                              kcat_b: float= DEFAULT_KCAT,
                              molmass:float = DEFAULT_MOLMASS) -> EnzymeInformation:
    if enzyme_id is None:
        enzyme_id = 'Enzyme_' + rxn_id
    return EnzymeInformation(enzyme_id,
                             check_kcat_value(kcat_f),
                             check_kcat_value(kcat_b),
                             genes,
                             protein_reaction_association,
                             molmass
                             )

def parse_gpr_information(gpr_info:str, gene:str,
                                          enzyme_id:str,
                                          gene2protein: dict[str, str]) -> tuple[list,list]:
    #filter out nan entries
    if not isinstance(gpr_info, str):
        return [['gene_dummy']], None
    gpr_list = _parse_gpr(gpr_info)

    #only get the genes associated with this enzyme
    gpr_list = _filter_sublists(gpr_list, gene)

    #convert the genes to the associated proteins
    enzyme_relations = []
    for sublist in gpr_list:
        enz_sublist = []
        for item in sublist:
            if item in gene2protein.keys():
                enz_sublist.append(gene2protein[item])
        enzyme_relations += [enz_sublist]
    # only get the enzyme relations which involves this enzyme
    enzyme_relations = _filter_sublists(enzyme_relations, enzyme_id)

    return gpr_list, enzyme_relations

def get_protein_gene_mapping(enzyme_db: pd.DataFrame, model) -> tuple[dict, dict]:
    protein2gene = defaultdict(list)
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
            if not isinstance(gene_id, str):
                gene_id = [gene.id for gene in rxn.genes][0]  # TODO

            gene2protein[gene_id] = enzyme_id

            # Create gene-protein-reaction associations
            if enzyme_id not in protein2gene:
                protein2gene[enzyme_id].append(gene_id)

    return dict(protein2gene), gene2protein

def _parse_gpr(gpr_info):
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

def _filter_sublists(nested_list, target_string):
    """
    Filters out all sublists from a nested list that contain the target string.

    Args:
        nested_list (list of list of str): The nested list to filter.
        target_string (str): The string to filter out sublists that contain it.

    Returns:
        list of list of str: A new nested list with the filtered sublists.
    """
    return [sublist for sublist in nested_list if target_string in sublist]


def check_kcat_value(kcat:float):
    return [kcat if not np.isnan(kcat) else 0][0]

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
        if rxn_id not in model.reactions:continue
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
