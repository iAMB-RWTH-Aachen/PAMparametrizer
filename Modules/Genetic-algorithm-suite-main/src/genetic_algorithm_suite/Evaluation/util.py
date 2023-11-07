"""
Utility functions

"""

from warnings import warn


# %%

def gene_to_reaction_knockouts(model, gpr_expressions: dict, gene_stat_dict: dict,
                               reactions: list = -1) -> list:
    """Determine all reactions which are affected by a set of gene knockouts
    
    Inputs:
        :param cobra.core.model model: metabolic model in COBRA format
        :param dict gpr_expressions: executable GPR expression for every model reaction
        :param dict gene_stat_dict: indicates gene deletions by False values
        :param list reactions: reactions to be evaluated
        
    Outputs:
        :param list reaction_knockouts: IDs of deleted reactions
    


    """
    
    reaction_knockouts = []
    
    # check all reactions if no reaction list was provided
    if reactions == -1:
        reactions = list(gpr_expressions.keys())

    # find dependent reaction knockouts according to GPR logical expressions
            
    # evaluate GPR expressions and determine reaction knockouts
    for r in reactions:
        expr = gpr_expressions[r]
        if not(eval(expr)):
            # reaction is disabled due to gene knockouts
            reaction_knockouts.append(r)
    
                                        
    return reaction_knockouts


def create_GPR_expressions(model):
    """Transform gene (protein) reaction rules (GPR) into machine readable and 
    executable expressions

    Inputs:
        :param cobra.core.model model: metabolic model in COBRA format
        
    Outputs:
        :param dict gpr_expressions: executable GPR expression for every model reaction
        :param dict gene_stat_dict: indicates gene deletions by False values
            

    """
    
    # helper functions
    def transform_gpr(gpr, gene_ids):
        # transform GPR into a machine readable form
        # set logical operators
        gpr_expr = gpr.replace("and", "&").replace("or", "|")
        # replace gene ids with calls to gene status array
        for g in gene_ids:
            gpr_expr = gpr_expr.replace(g, "gene_stat_dict['" + g + "']")

        return gpr_expr        
    
    # allocate gene status array
    gene_stat_dict = {g.id: True for g in model.genes}

    
    gpr_expressions = {}  
    for r in model.reactions:
        gpr = r.gene_reaction_rule
        if len(gpr) == 0: continue # no GPR available for reaction
        
        # determine IDs of genes involved in GPR
        gene_ids = [g.id for g in r.genes]
        # transform GPR
        gpr_expressions[r.id] = transform_gpr(gpr, gene_ids)
    
    
    
    return gpr_expressions, gene_stat_dict









