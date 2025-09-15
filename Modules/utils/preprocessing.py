import pandas as pd
import numpy as np
import re
from cobra import Model as CobraModel
from typing import List, Dict, Union

DEFAULT_KCAT = 13.7 #s-1, median from Brenda as determined by Bar-Evan et al. (2016)
DEFAULT_MOLMASS = 39959.4825 #g/mol
DEFAULT_PROT_LENGTH = 325 # amino acids, data from E. coli (Zhang J. Protein-length distributions for the three domains of life. Trends Genet. 2000 Mar16(3):107-9.PubMed ID10689349)

def extract_locus_tags(text:str,
                       locustag_regex:str) -> List[str]:
    if pd.isna(text):
        return text
    return re.findall(locustag_regex, text)


def create_id_mapper_from_model(model: CobraModel,
                                rxn_annotation_keys: list[str]=['kegg.reaction', 'ec-code', ],
                                met_annotation_key: str='kegg.compound',
                                exclude:List[str]=['Growth', 'ATPM', 'BIOMASS']
                                ) -> pd.DataFrame:
    """
    Create a mapping DataFrame from a COBRA model's reactions.

    Parameters
    ----------
    model : cobra.Model
        The COBRA model from which to extract reaction information.
    rxn_annotation_keys : list of str, optional
        Keys to extract from the reaction annotations.
    met_annotation_key : str, optional
        Key to extract from metabolite annotations (used for both reactants and products).
    exclude: str, optional
        Reactions not assoiated with a protein, for example growth and maintenance pseudo reactions

    Returns
    -------
    pd.DataFrame
        A DataFrame with columns for reaction ID, annotations, reactants, products, and GPR.
    """
    data = []
    mapped_to_kegg = 0

    for rxn in model.reactions:
        #check if reaction is passive transport
        if any(ex in rxn.id for ex in ['tpp',r't[0-9]pp']):
            if has_same_reactants_and_products(rxn):continue
        if any(ex in rxn.id for ex in ['tex', 'EX_','sink']) or rxn.id in exclude:
            continue

        entry = {
            'rxn_id': rxn.id,
            'Reactants': [met.annotation.get(met_annotation_key, np.nan) for met in rxn.reactants],
            'Products': [met.annotation.get(met_annotation_key, np.nan) for met in rxn.products],
            'GPR': rxn.gene_reaction_rule,
            'reversible': rxn.reversibility
        }

        for key in rxn_annotation_keys:
            value = rxn.annotation.get(key, np.nan)
            entry[key] = value
            if key == 'kegg.reaction' and value is not np.nan:
                mapped_to_kegg += 1

        data.append(entry)

    print(f"Mapped {mapped_to_kegg} of {len(data)} reactions to KEGG reaction IDs.")

    return pd.DataFrame(data)

def has_same_reactants_and_products(reaction):
    """
    Check if products and reactants of a COBRApy reaction
    are the same metabolites regardless of compartment.
    """
    # Get sets of metabolite "base names" (without compartment suffix)
    reactants = {m.formula for m in reaction.reactants}
    products  = {m.formula for m in reaction.products}

    return reactants == products

def create_genetokeggid_mapper(model:CobraModel) -> pd.DataFrame:
    data = []
    mapped_to_kegg = 0

    for gene in model.genes:
        if not 'kegg.genes' in gene.annotation: continue
        entry = {
            'gene_id': gene.id,
            'kegg_gene_id': gene.annotation.get('kegg.genes'),
        }
        mapped_to_kegg += 1

        data.append(entry)
    print(f'Mapped {mapped_to_kegg} out of {len(model.genes)} genes to a KEGG identifier')
    return pd.DataFrame(data)

def replace_locustags_in_text(text:str,
                              id_map:Dict[str, str]
                              )-> str:
    for old_id, new_id in id_map.items():
        text = text.replace(old_id, new_id)
    return text

def map_kcat_values_to_reaction_protein_association(id_mapper: pd.DataFrame,
                                                    gotenzymes_df: pd.DataFrame
                                                    ) -> pd.DataFrame:
    """Merges kcat data from GotEnzymes into model annotation based on gene ID or EC number.

    Args:
        id_mapper (pd.DataFrame): Model annotation with 'rxn_id', 'locus_tag', 'kegg.reaction', and 'EC'.
        gotenzymes_df (pd.DataFrame): GotEnzymes data with 'gene', 'reaction_id', 'ec_number', 'kcat_values'.

    Returns:
        pd.DataFrame: Annotated reactions with available kcat values and model gene IDs.
    """
    # Match based on gene ID and KEGG reaction
    merged_by_gene = pd.merge(
        id_mapper, gotenzymes_df,
        left_on=['locus_tag', 'kegg.reaction'],
        right_on=['gene', 'reaction_id'],
        how='left'
    )
    mapped_gene = merged_by_gene.dropna(subset=['kcat_values'])
    unmapped = merged_by_gene[merged_by_gene['kcat_values'].isna()]
    print(f'Mapped {len(mapped_gene )} '
          f'out of {len(merged_by_gene)} kcat values to reactions based on kegg gene id')

    # Match remaining by EC number and reaction ID
    merged_by_ec = pd.merge(
        unmapped[[col for col in unmapped if col not in gotenzymes_df.columns]],
        gotenzymes_df,
        left_on=['EC', 'kegg.reaction'],
        right_on=['ec_number', 'reaction_id'],
        how='left'
    )
    print(f'Mapped {len(merged_by_ec)} '
          f'out of {len(merged_by_gene)} kcat values to reactions based on ec number')

    # Keep the original gene/locus_tag from the model
    merged_by_ec['gene'] = merged_by_ec['locus_tag']

    # Remove duplicates already in mapped_gene
    already_mapped_rxns = set(mapped_gene['rxn_id'])
    merged_by_ec = merged_by_ec[~merged_by_ec['rxn_id'].isin(already_mapped_rxns)]

    # Combine and clean
    combined = pd.concat([mapped_gene, merged_by_ec], ignore_index=True)
    combined = combined.drop(['reaction_id', 'locus_tag'], axis=1)

    print(f'Final merged dataset has {len(combined.dropna(subset=["kcat_values"]))}'
          f' unique reactions with kcat values.')

    return combined


def assign_missing_gprs(df: pd.DataFrame, use_ec: bool = False) -> pd.DataFrame:
    """Assigns default gene, GPR, and enzyme IDs for unmappable proteins.

    This function updates the 'gene', 'GPR', and 'enzyme_id' columns for rows where
    'GPR' is missing (empty string). If `use_ec=True`, EC numbers are also used
    to generate gene and enzyme identifiers.

    Args:
        df (pd.DataFrame): A DataFrame containing enzyme reaction data,
            including 'GPR', 'rxn_id', 'gene', 'enzyme_id', and optionally 'EC'.
        use_ec (bool, optional): Whether to include EC numbers in the assigned values.
            Defaults to False.

    Returns:
        pd.DataFrame: The updated DataFrame with missing 'gene', 'GPR', and 'enzyme_id'
        values assigned where applicable.
    """
    # Ensure GPR is treated as a string and replace NaNs with empty strings
    df['GPR'] = df['GPR'].fillna('')

    # Condition: GPR is missing
    missing_gene_mask = (df['GPR'] == '')

    if use_ec:
        missing_gene_mask &= df['EC'].notna()
        id_pattern = lambda row: (f'Gene_{row.rxn_id}_{row.EC}',
                                  f'Gene_{row.rxn_id}_{row.EC}',
                                  f'Enzyme_{row.rxn_id}_{row.EC}')
    else:
        id_pattern = lambda row: (f'Gene_{row.rxn_id}',
                                  f'Gene_{row.rxn_id}',
                                  f'Enzyme_{row.rxn_id}')

    # Apply only to rows that meet the condition
    new_values = df.loc[missing_gene_mask].apply(lambda row: pd.Series(id_pattern(row)), axis=1)
    # Assign values
    if len(new_values)>0:
        df.loc[missing_gene_mask, ['gene', 'GPR', 'enzyme_id']] = new_values.to_numpy()

    return df

def assign_directionalities_for_kcat_relations(eco_enzymes: pd.DataFrame
                                               ) -> pd.DataFrame:
    """
        Assigns directionality ('f' for forward, 'b' for backward) to kcat entries and duplicates
        reversible reactions that do not already include a backward direction.

        The direction is inferred by checking whether the compound associated with the kcat value
        appears in the 'Products' list of the reaction. If so, the direction is 'b' (backward),
        otherwise 'f' (forward). For reversible reactions, if no 'b' direction is found, a
        backward copy is created.

        Args:
            eco_enzymes (pd.DataFrame): A DataFrame with at least the following columns:
                - 'kcat_values'
                - 'compound'
                - 'Products' (list of compounds)
                - 'Reactants' (list of compounds)
                - 'reversible' (bool)

        Returns:
            pd.DataFrame: An updated DataFrame including direction assignments and, if needed,
                          additional rows for missing backward directions.
        """
    # Ensure Reactants and Products are lists
    eco_enzymes = eco_enzymes.copy()
    eco_enzymes['Products'] = eco_enzymes['Products'].apply(lambda x: x if isinstance(x, list) else [])
    eco_enzymes['Reactants'] = eco_enzymes['Reactants'].apply(lambda x: x if isinstance(x, list) else [])

    # Assign a direction: 'b' if compound in Products, otherwise 'f'
    eco_enzymes['direction'] = np.where(
        (eco_enzymes['kcat_values'].notna() &
         [row.compound in row.Products for i, row in eco_enzymes.iterrows()]),
        'b', 'f'
    )

    # Find reversible reactions missing a backward direction
    reversible_groups = eco_enzymes[eco_enzymes['reversible']].groupby('rxn_id')
    rows_to_add = []

    for rxn_id, group in reversible_groups:
        if 'b' not in group['direction'].values:
            for _, row in group.iterrows():
                row_copy = row.copy()
                row_copy['direction'] = 'b'
                rows_to_add.append(row_copy)

    # Append missing backward reactions
    if rows_to_add:
        eco_enzymes = pd.concat([eco_enzymes, pd.DataFrame(rows_to_add)], ignore_index=True)

    return eco_enzymes

def assign_defaults_for_proteins_without_mapping(eco_enzymes: pd.DataFrame,
                                                 default_kcat:float = DEFAULT_KCAT,
                                                 default_molmass:Union[float, int] = DEFAULT_MOLMASS,
                                                 default_protein_length: int = DEFAULT_PROT_LENGTH
                                                 )-> pd.DataFrame:
    """
        Assigns default values for unmapped protein properties in a kcat annotation DataFrame.

        This function ensures that all enzymes have valid values for key fields such as
        `kcat_values`, `molMass`, `Length`, `gene`, and `enzyme_id`. Default values are filled in
        where data is missing. It also attempts to infer gene and enzyme identifiers based on
        available EC numbers or existing values.

        Special handling is included for unmapped gene annotations such as "s0001".

        Args:
            eco_enzymes (pd.DataFrame):
                A DataFrame containing enzyme annotations, including columns like
                'kcat_values', 'Mass', 'Length', 'gene', 'enzyme_id', 'GPR', 'rxn_id',
                'ec_number', and 'compound'.
            default_kcat (float, optional):
                Default kcat value to assign where missing. Defaults to `DEFAULT_KCAT`.
            default_molmass (float or int, optional):
                Default molecular mass (in Da) to assign where missing. Defaults to `DEFAULT_MOLMASS`.
            default_protein_length (int, optional):
                Default protein length to assign where missing. Defaults to `DEFAULT_PROT_LENGTH`.

        Returns:
            pd.DataFrame:
                The updated DataFrame with defaults filled in, GPR assignments updated,
                enzyme and gene identifiers inferred, and final cleanup of unused columns.

        Notes:
            - Rows with GPR value `"s0001"` are treated as special unmapped cases and assigned
              enzyme IDs using their `rxn_id`.
            - Relies on `assign_missing_gprs`, which must be defined in the namespace.
    """

    eco_enzymes['kcat_values'] = eco_enzymes['kcat_values'].fillna(default_kcat)
    eco_enzymes['Mass'] = eco_enzymes['Mass'].fillna(default_molmass)
    eco_enzymes['Length'] = eco_enzymes['Length'].fillna(default_protein_length)
    # Assign default genes for unmappable proteins
    # if there is an ec number available
    eco_enzymes = assign_missing_gprs(eco_enzymes, use_ec=True)
    # if there is no information available
    eco_enzymes = assign_missing_gprs(eco_enzymes, use_ec=False)

    # Assign enzyme IDs if missing but gene information is available
    eco_enzymes['enzyme_id'] = np.where(
        eco_enzymes['enzyme_id'].isna(),
        eco_enzymes['gene'].apply(lambda x: f'Enzyme_{x}' if isinstance(x, str) else None),
        eco_enzymes['enzyme_id']
    )

    # Assign genes if missing if enzyme information is available
    eco_enzymes['gene'] = eco_enzymes.apply(
        lambda row: f'gene_{row.enzyme_id}' if isinstance(row.gene, float) else row.gene, axis=1
    )

    # Handle the specific case for 's0001' (annotation is gene is not mapped to a genome annotation)
    eco_enzymes.loc[eco_enzymes['GPR'] == 's0001', 'gene'] = eco_enzymes['GPR']
    eco_enzymes.loc[eco_enzymes['GPR'] == 's0001', 'enzyme_id'] = 'Enzyme_s0001' + eco_enzymes[
        'rxn_id']

    # Clean up final columns
    eco_enzymes = eco_enzymes.drop(
        ['ec_number', 'compound'], axis=1
    ).rename(columns={'Mass': 'molMass'})
    return eco_enzymes