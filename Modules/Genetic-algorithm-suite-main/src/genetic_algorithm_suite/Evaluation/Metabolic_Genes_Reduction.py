"""
Simple evaluation class and operators for the reduction of metabolic genes from a metabolic model

"""
import cobra
import random
import pandas as pd
import numpy as np
from multiprocessing import Pool
from pathlib import Path
from os.path import dirname, abspath

from genetic_algorithm_suite.Evaluation.util import create_GPR_expressions, gene_to_reaction_knockouts

# set standard paths
FILE_PATH = Path(abspath(dirname(__file__)))
DATA_PATH = FILE_PATH.parents[0].joinpath("Data")

# seed random number generator
random.seed()

class FitnessEvaluation():
    
    def __init__(self, model=None, fixed_attr_list=[], processes=2):
        """Initialize fitness evaluation class for a genetic algorithm
        
        
        """
        # save processes
        self.processes = processes
        
        # save list of fixed attributes
        self.fixed_attr_list = fixed_attr_list
        
        
        # initialize general parameters
        self._init_general_parameters()
        

        # initialize the model
        self.init_model(model)
      
    

    ##########################################################################
    # NECESSARY CUSTOM FUNCTION
    ##########################################################################
    
    def init_model(self, model=None):
        """Initialize all necessary parameter for the fitness function
        
        """
        if not model:
            # load default model 
            self.model = cobra.io.load_json_model(str(FILE_PATH.joinpath("../Models/iML1515.json")))
        else:
            self.model = model
        
        # set objective
        self.model.objective = "BIOMASS_Ec_iML1515_WT_75p37M"
        
        # enable WT Biomass equation
        self.model.reactions.EX_adocbl_e.lower_bound = -1 # uptake of Adenosylcobalamin obligatory
        
        
        # set constraints
        self.model.reactions.EX_o2_e.lower_bound = -20 # oxygen uptake limit
        self.model.reactions.EX_co2_e.lower_bound = 0 # avoid CO2 uptake
        
        
        # create machine readable GPR expressions for linking gene to reaction deletions
        self.gpr_expressions, self.gene_stat_dict = create_GPR_expressions(self.model)
        
        # improve constraints
        # restrict lower/upper bounds to -100/100
        for r in self.model.reactions:
            if r.upper_bound > 100: r.upper_bound = 100
            if r.lower_bound < -100: r.lower_bound = -100
                    
        
        # add exchange reactions for lipid recycling
        fa_id = ["pe160_p", "pe161_p", "pe181_p", "pg160_p", "pg161_p", "pg181_p",
              "clpn160_p", "clpn161_p", "clpn181_p"]
        self.fatty_acids_exchanges = []
        for fa in fa_id:
            ex_id = "EX_"+fa
            rxn_ex = cobra.Reaction(id=ex_id, lower_bound=0, upper_bound=0)
            rxn_ex.add_metabolites({self.model.metabolites.get_by_id(fa): -1})
            self.model.add_reaction(rxn_ex)
            self.fatty_acids_exchanges.append(ex_id)
            
            
        # add exchanges for murein/peptidoglycan recycling
        pc_id = ["murein3p3p_p", "murein3px4p_p", "murein4p4p_p", "murein4px4p_p",
                 "murein4px4px4p_p"]
        self.peptidoglycans_exchanges = []
        for pc in pc_id:
            ex_id = "EX_"+pc
            rxn_ex = cobra.Reaction(id=ex_id, lower_bound=0, upper_bound=0)
            rxn_ex.add_metabolites({self.model.metabolites.get_by_id(pc): -1})
            self.model.add_reaction(rxn_ex)
            self.peptidoglycans_exchanges.append(ex_id)
            
            
        key_precursors = ["chor_c", "accoa_c", "malcoa_c"]
        # add exchange reactions for key product precursors
        for p in key_precursors:
            rxn_ex = cobra.Reaction(id="EX_"+p, lower_bound=0, upper_bound=0)
            rxn_ex.add_metabolites({self.model.metabolites.get_by_id(p): -1})
            self.model.add_reaction(rxn_ex)
            
            
        # create one exchange reaction for amino acids in cytosol (mimic protein production)
        aa = ["his__L_c", "leu__L_c", "ile__L_c", "tyr__L_c",
                "phe__L_c", "trp__L_c", "ala__L_c",
                "gly_c", "glu__L_c", "gln__L_c", "val__L_c", 
                "arg__L_c", "asp__L_c", "asn__L_c", "cys__L_c", 
                "lys__L_c", "met__L_c", "pro__L_c", "ser__L_c", 
                "thr__L_c",
            ]
        rxn_ex = cobra.Reaction(id="EX_generic_protein", lower_bound=0, upper_bound=0)
        for a in aa:
            # add amino acids with stoichiometries from the biomass equation
            stoich = self.model.reactions.get_by_id("BIOMASS_Ec_iML1515_WT_75p37M").metabolites[self.model.metabolites.get_by_id(a)]
            rxn_ex.add_metabolites({self.model.metabolites.get_by_id(a): stoich})
        self.model.add_reaction(rxn_ex)
            
        # compute connectivities of metabolites
        # only consider metabolite with a connectivity greater than 2 (greater than the minimally possible)
        self.metabolite_connectivity = {m.id: set([r.id for r in m.reactions]) for m in self.model.metabolites}
        self.std_connectivity_WT = np.std([len(r) for m, r in self.metabolite_connectivity.items()])


        # compute wildtype functionalities
        valid_fun, sol_dict = self._compute_metabolic_functionalities(self.model)

        # compute minimum acceptaple objective rates
        self.sol_ref_dict = {}
        for key, sol in sol_dict.items():
            self.sol_ref_dict[key] = sol * self.minimum_relative_rate



    def compute_individuals_properties(self, pop) -> list:
        """Compute custom properties of a population's individuals. For each 
        individual return a dictionary with properties as keys
        
        Inputs:
            :param list pop: Individuals of a population
         
        Outputs:
            :param list ind_props: Property values for each individual. Return
                                    an empty list if no properties should be computed
        
        """
        ind_props = []
        
        
        for ind in pop:
            props_dict = {}
            
            # number of deleted genes
            props_dict["number_deleted_genes"] = sum([i==0 for i in ind])
            # standard deviation metabolite connectivities
            props_dict["std_metabolite_connectivity"] = ind.fitness.values[1]
            
            # save properties
            ind_props.append(props_dict)
        
        
        return ind_props

    ##########################################################################
    # DEAP-RELATED FUNCTIONS
    ##########################################################################
    
    
    def attribute_generator(self, probability=0) -> int:
        """Generates an attribute of a DEAP individual
        Inputs:
            :param float probability: probability for returning a "0" as attr
        
        Outputs:
            :param int attr: representation of an attribute
        
        """
        
        if random.random() >= probability:
            attr = 1
        else:
            attr = 0
        
        return attr
    
    
    def init_fitness(self) -> dict:
        """Specifies parameters to set up a fitness function in DEAP
        - weights of fitness functions
        
        Outputs:
            :param dict param_fit: parameters to set up a DEAP fitness function
        
        """
        
        param_fit = {}
        # defines the weighting of the fitness functions
        param_fit["weights"] = (-1, -1)
        
        return param_fit
        
    
    def init_attribute(self) -> dict:
        """Specifies parameters to describe attributes of individuals in DEAP
        - number of attributes
        - type of an attribute (metabolic gene, regulator etc.)
        
        Inputs:
            
        Outputs:
            :param dict param_attr: 
        
        """
        
        
        param_attr = {}
        
        # get protected genes
        protected_genes = self._determine_protected_genes()
            
        # only consider non-protected genes
        gene_id = [g.id for g in self.model.genes if g.id not in protected_genes]
          
        # compute non-essential genes
        with Pool(processes=self.processes) as pool:
            results = pool.map(self._determine_gene_essentiality, gene_id)
        non_essential_genes = [gene for result in results for gene in result]    
            
        non_essential_id = [g["id"] for g in non_essential_genes]
        
        # remove fixed attributes not in non-essential genes
        for attr in self.fixed_attr_list:
            if attr["id"] not in non_essential_id:
                print("Fixed attribute", attr["id"], "not in target space")
                
        self.fixed_attr_list = [attr for attr in self.fixed_attr_list if attr["id"] in non_essential_id]
        
        
        # get fixed attribute IDs
        fixed_attr_list_id = [fa["id"] for fa in self.fixed_attr_list]

        # save attribute details
        self.individual_attr_list = [g for g in non_essential_genes if g["id"] not in fixed_attr_list_id]
        
        print("Number of attributes:", len(self.individual_attr_list))
                
        
        # number of attribute per individual
        param_attr["number_attributes"] = len(self.individual_attr_list)
        
        
        return param_attr
    
    
    
    
    def init_individual(self) -> dict:
        """Set up charactersitics and parameters of individuals in DEAP
        
        Outputs:
            :param dict param_ind: parameters for setting up an individual 

        """
        
        param_ind = {}
        

        # create individuals represnting metabolic genes as variables/targets
        # creator.create("Individual", list, fitness=creator.FitnessMin)
        param_ind["individual_type"] = list
        
        return param_ind
            

    
    
    
    ##########################################################################
    # FITNESS FUNCTION
    ##########################################################################
    
    def eval_fitness(self, individual) -> int:
        """Evaluate the fitness of an individual
        
        Inputs:
            :param list individual: a list of variables representing an individual (solution)
            :param cobra.core.model model: metabolic model in cobrapy format
            
        Outputs:
            :param int fitness: fitness of the current individual evaluated
        
        """

        # apply knockouts and compute metabolic functionalities        
        with self.model as model:
              
            # determine reaction knockouts from gene knockouts regarding GPR logical rules
            reactions_affected = []
            for i in range(len(individual)):
                ind_attr = self.individual_attr_list[i]
                if ind_attr["type"] == "m_gene":
                    if individual[i] == 1:
                        # gene is intact
                        self.gene_stat_dict[ind_attr["id"]] = True
                    else:
                        # gene is deleted
                        self.gene_stat_dict[ind_attr["id"]] = False
                        # determine possibly affected reactions
                        reactions_affected.extend([r.id for r in model.genes.get_by_id(ind_attr["id"]).reactions])
            
            # add fixed attributes
            for attr in self.fixed_attr_list:
                if attr["type"] == "m_gene":
                    self.gene_stat_dict[attr["id"]] = attr["value"]
                    if not attr["value"]:
                        reactions_affected.extend([r.id for r in model.genes.get_by_id(attr["id"]).reactions])
                 
           
            reactions_affected = list(np.unique(reactions_affected))
                      
            # get reaction knockouts
            reaction_knockouts = gene_to_reaction_knockouts(
                                    model,
                                    self.gpr_expressions,
                                    self.gene_stat_dict,
                                    reactions=reactions_affected)
                       
            # set reaction knockouts
            for r in reaction_knockouts:
                model.reactions.get_by_id(r).bounds = (0, 0)
            

            # compute protected metabolic functionalities 
            valid_fun, sol_dict = self._compute_metabolic_functionalities(model, self.sol_ref_dict)

            
        # calculate fitness
        if valid_fun:            
            # count number of intact genes
            fit_gene_number = sum(individual)
            # set of reaction knockouts           
            reaction_knockouts_set = set(reaction_knockouts)
            # get metabolite's connectivities
            connectivity = [len(r - reaction_knockouts_set) for m, r in self.metabolite_connectivity.items()]
            fit_connectivity = np.std(connectivity)

        else:
            # infeasible solution, assume worst fitness value
            fit_gene_number = len(individual)
            # assume wiltype connectivity
            fit_connectivity = self.std_connectivity_WT


        # return a tuple of one element
        return tuple([float(fit_gene_number), float(fit_connectivity)])
        
   
 
    ##########################################################################
    # CUSTOM HELPER FUNCTIONS
    ##########################################################################
      
    
    def _compute_metabolic_functionalities(self, model, sol_ref_dict={}):
        """Compute protected functionalities of the metabolic network and ensure
        their activity
        Return objective values for each functionality
        
        Inputs:
            :param cobra.core.model model: metabolic model
            :param dict sol_ref_dict: reference solutions for each functionality
            
        Outputs:
            :param bool valid: True if all metabolic functionalities are active
            
        """
        # 
        valid = True    
        
        # save objective values
        sol_dict = {}
        
        program = ""
        
        
        # test for very basic functionality
        if program == "simple":
            key_fun = "simple"
            
            sol = model.slim_optimize(error_value=0)
            # save solution
            sol_dict[key_fun] = sol
     
                
            if key_fun in sol_ref_dict:
                if sol < sol_ref_dict[key_fun]:
                    # invalid function
                    valid = False
        
        
            elif sol < self.minimum_growth_rate:
                # invalid function
                valid = False
                
            return valid, sol_dict
           
        ##########################################
        # MINIMAL MEDIUM, DIFFERENT CARBON SOURCES
        ##########################################
                 
        
        source_targets = ['EX_arab__L_e','EX_acgam_e',
                     'EX_succ_e',
                     'EX_gal_e',
                     'EX_asp__L_e',
                     'EX_pro__L_e',
                     'EX_ala__D_e',
                     'EX_tre_e',
                     'EX_man_e',
                     'EX_ser__D_e',
                     'EX_sbt__D_e',
                     'EX_glyc_e',
                     'EX_fuc__L_e',
                     'EX_glcur_e',
                     'EX_glyc3p_e',
                     'EX_xyl__D_e',
                     'EX_lac__L_e',
                     'EX_mnl_e',
                     'EX_g6p_e',
                     'EX_mal__L_e',
                     'EX_mal__D_e',
                     'EX_rib__D_e',
                     'EX_rmn_e',
                     'EX_fru_e',
                     'EX_ac_e',
                     'EX_glc__D_e',
                     'EX_malt_e',
                     'EX_melib_e',
                     'EX_thymd_e',
                     'EX_asn__L_e',
                     'EX_akg_e',
                     'EX_lcts_e',
                     # 'EX_sucr_e',
                     'EX_uri_e',
                     'EX_gln__L_e',
                     'EX_g1p_e',
                     'EX_f6p_e',
                     'EX_malttr_e',
                     'EX_dad_2_e',
                     'EX_adn_e',
                     'EX_fum_e',
                     'EX_ppa_e',
                     'EX_glyclt_e',
                     'EX_ins_e',
                     'EX_ser__L_e',
                     'EX_thr__L_e',
                     'EX_ala__L_e',
                     'EX_acac_e',
                     'EX_acmana_e',
                     'EX_lyx__L_e',
                     'EX_pyr_e',
                     'EX_galur_e',
                     'EX_all__D_e',
                     'EX_acnam_e',
                     'EX_gam_e',
                     'EX_dha_e'] # growth was checked on BioCyc (Biolog PM1/2)  
        
        # source_targets = ['EX_arab__L_e', 'EX_glc__D_e', 'EX_succ_e',]
        
        
        with model as model:
            
            # protect model
            with model:
                # precaution, disable glucose uptake
                model.reactions.EX_glc__D_e.lower_bound = 0
                # scan substrates 
                for s in source_targets:
                    # functionality key
                    key_fun = s
                    ex_rxn = model.reactions.get_by_id(s)
                    # define uptake rate (equivalent of 60 C-mmol/gDW/h)
                    uptake_rate = 60 / list(ex_rxn.metabolites.keys())[0].elements["C"]
                    ex_rxn.lower_bound = -uptake_rate
                    # solve model
                    
                    sol = model.slim_optimize(error_value=0)
                    # save solution
                    sol_dict[key_fun] = sol
                    # revert substrate uptake
                    ex_rxn.lower_bound = 0
                    
                    if key_fun in sol_ref_dict:
                        if sol < sol_ref_dict[key_fun]:
                            # invalid function
                            valid = False
                            break
                        
                    elif sol < self.minimum_growth_rate:
                        # invalid function
                        valid = False
                        break
                    
                    ###########
                    # break

                # return valid, sol_dict
            
            ##########################################
            # (Micro) ANAEROBIC GROWTH ON GLUCOSE MINIMAL MEDIUM
            ##########################################    
            if valid:
                with model:
                    # set key of function
                    key_fun = "microaerobic"
                    # enable uptake of glucose 
                    model.reactions.EX_glc__D_e.lower_bound = -10
                    # disable o2 uptake
                    wt_o2_uptake_lb = model.reactions.EX_o2_e.lower_bound
                    model.reactions.EX_o2_e.lower_bound = -1
                    # solve model
                    sol = model.slim_optimize(error_value=0)
                    # save solution
                    sol_dict[key_fun] = sol
                    
                    # revert oxygen uptake rate
                    model.reactions.EX_o2_e.lower_bound = wt_o2_uptake_lb
                    
                    if key_fun in sol_ref_dict:
                        if sol < sol_ref_dict[key_fun]:
                            # invalid function
                            valid = False
                    
                    elif sol < self.minimum_growth_rate:
                        # invalid function
                        valid = False
                    
                
            ##################################################################
            # Amino acid uptake and metabolism (without glucose as carbon source)
            ################################################################## 
            if valid:
                with model:
                    key_fun = "AA_medium"
                    
                    aa_ex = ["EX_his__L_e", "EX_leu__L_e", "EX_ile__L_e", "EX_tyr__L_e",
                             "EX_phe__L_e", "EX_trp__L_e", "EX_ala__L_e",
                             "EX_gly_e", "EX_glu__L_e", "EX_gln__L_e", "EX_val__L_e", 
                             "EX_arg__L_e", "EX_asp__L_e", "EX_asn__L_e", "EX_cys__L_e", 
                             "EX_lys__L_e", "EX_met__L_e", "EX_pro__L_e", "EX_ser__L_e", 
                             "EX_thr__L_e",
                      ]
                    
                    model.reactions.EX_glc__D_e.lower_bound = 0
                    model.reactions.EX_o2_e.lower_bound = -100
                    
                    # set amino acid uptake
                    for aa in aa_ex:
              
                        model.reactions.get_by_id(aa).lower_bound = -5
                        
                    # solve model
                    sol = model.slim_optimize(error_value=0)   
                    # save solution
                    sol_dict[key_fun] = sol
                    
                    
                    
                    if key_fun in sol_ref_dict:
                        if sol < sol_ref_dict[key_fun]:
                            # invalid function
                            valid = False
                    
                    elif sol < self.minimum_growth_rate:
                        # invalid function
                        valid = False   
                    
                
            ##################################################################
            # recycling of fatty acids 
            ################################################################## 
            if valid:
                with model:
                    key_fun = "fatty_acids_recycling"
                    # set low glucose uptake
                    model.reactions.EX_glc__D_e.lower_bound = 0
                    # increase o2 uptake
                    model.reactions.EX_o2_e.lower_bound = -100
                    # activate fatty acid exchanges
                    for fa in self.fatty_acids_exchanges:
                        model.reactions.get_by_id(fa).lower_bound = -5
                        model.reactions.get_by_id(fa).upper_bound = -1
                
                    # solve model
                    sol = model.slim_optimize(error_value=0) 
       
                    # save solution
                    sol_dict[key_fun] = sol
                    
                    if key_fun in sol_ref_dict:
                        if sol < sol_ref_dict[key_fun]:
                            # invalid function
                            valid = False
                    
                    elif sol < self.minimum_growth_rate:
                        # invalid function
                        valid = False 
           
            ##################################################################
            # recycling of murein/peptidoglycan 
            ################################################################## 
            if valid:
                with model:
                    key_fun = "peptidoglycan_recycling"
                    # set low glucose uptake
                    model.reactions.EX_glc__D_e.lower_bound = 0
                    # increase o2 uptake
                    model.reactions.EX_o2_e.lower_bound = -100
                    # activate fatty acid exchanges
                    for pc in self.peptidoglycans_exchanges:
                        model.reactions.get_by_id(pc).lower_bound = -0.5
                        model.reactions.get_by_id(pc).upper_bound = -0.1
                
                    # solve model
                    sol = model.slim_optimize(error_value=0) 
       
                    
                    # save solution
                    sol_dict[key_fun] = sol
                    
                    if key_fun in sol_ref_dict:
                        if sol < sol_ref_dict[key_fun]:
                            # invalid function
                            valid = False
                    
                    elif sol < self.minimum_growth_rate:
                        # invalid function
                        valid = False 
                        
                        
                        
            ##################################################################
            # recycling of nucleotides
            ################################################################## 
            if valid:
                with model:
                    key_fun = "nucleotide_recycling"
                    
                    # define nucleotide exchange reactions
                    ncl_ex = ["EX_amp_e", "EX_damp_e", "EX_gmp_e", "EX_dgmp_e", 
                              "EX_cmp_e", "EX_dcmp_e", "EX_dtmp_e", "EX_ump_e"]
        
                    # set low glucose uptake
                    model.reactions.EX_glc__D_e.lower_bound = 0

                    # activate fatty acid exchanges
                    for ncl in ncl_ex:
                        model.reactions.get_by_id(ncl).lower_bound = -1
                
                    # solve model
                    sol = model.slim_optimize(error_value=0) 
       
                    
                    # save solution
                    sol_dict[key_fun] = sol
                    
                    if key_fun in sol_ref_dict:
                        if sol < sol_ref_dict[key_fun]:
                            # invalid function
                            valid = False
                    
                    elif sol < self.minimum_growth_rate:
                        # invalid function
                        valid = False 
                 
            ##################################################################
            # precursor production of value-added chemicals
            ##################################################################             
            # pyruvate, succinate, chorismate, acetyl-coa, manoyl-coa
            
            exchange_products = ["EX_pyr_e", "EX_succ_e", "EX_chor_c", "EX_accoa_c",
                                 "EX_malcoa_c"]
            
            if valid:
               with model:
                   # glucose minimal medium
                   model.reactions.EX_glc__D_e.lower_bound = -10
                   
                   for ex in exchange_products:
                        key_fun = ex+"_production"
                        # set constraints
                        model.reactions.get_by_id(ex).upper_bound = 100
                        # set objective function
                        model.objective = ex
                        # optimize model
                        sol = model.slim_optimize(error_value=0)
                        # save solution
                        sol_dict[key_fun] = sol
                        
                        # evaluate solution
                        if key_fun in sol_ref_dict:
                            if sol < sol_ref_dict[key_fun]:
                                # invalid function
                                valid = False
                                break
                            
                            elif sol < self.minimum_growth_rate:
                                # invalid function
                                valid = False 
                                break
                            
            
            ##################################################################
            # Protein production in form of amino acids overproduction
            ################################################################## 
            if valid:
                with model:
                    key_fun = "AA_production"
                    
                    model.reactions.EX_glc__D_e.lower_bound = -10
                    model.reactions.EX_o2_e.lower_bound = -100
                    
                    # allow protein production
                    model.reactions.get_by_id("EX_generic_protein").upper_bound = 100
                    # enforce minimal growth
                    model.reactions.get_by_id("BIOMASS_Ec_iML1515_WT_75p37M").lower_bound = 0.2
                    
                    # set objective to generic protein production reaction
                    model.objective = "EX_generic_protein"
                   
                    # solve model
                    sol = model.slim_optimize(error_value=0)   
                    # save solution
                    sol_dict[key_fun] = sol
                    
                    
                    
                    if key_fun in sol_ref_dict:
                        if sol < sol_ref_dict[key_fun]:
                            # invalid function
                            valid = False
                    
                    elif sol < self.minimum_growth_rate:
                        # invalid function
                        valid = False   
     
        return valid, sol_dict
    
    
    
    def _determine_protected_genes(self):
        """Determine genes which are protected and should not be deleted
        
        """
        
        protected_genes = []
        
        # protect amino acid transport to cytosol
        aa = ["his__L", "leu__L", "ile__L", "tyr__L", "phe__L", "trp__L", "ala__L",
              "gly", "glu__L", "gln__L", "val__L", "arg__L", "asp__L", "asn__L",
              "cys__L", "lys__L", "met__L", "pro__L", "ser__L", "thr__L",
              ]
        
        aa_rxns = []
        for a in aa:
            aa_c = self.model.metabolites.get_by_id(a + "_c")
            aa_p = self.model.metabolites.get_by_id(a + "_p")
            aa_e = self.model.metabolites.get_by_id(a + "_e")
            # get all reactions acting on amino acids in the periplasma and extracellulaar space
            for r in self.model.reactions:
                if ((aa_c in r.metabolites) and (aa_p in r.metabolites)) \
                    or ((aa_e in r.metabolites) and (aa_p in r.metabolites)):
                    aa_rxns.append(r)

        
        # get all genes in amino acid transport
        for r in aa_rxns:
            for g in r.genes:
                if g.id not in protected_genes:
                    protected_genes.append(g.id)
         
                    
        # protect experimentally validated essential genes
        e_genes_frame = pd.read_excel(DATA_PATH.joinpath("Essential_genes_Ecoli_MG1655_uniform.xlsx"),
                                      sheet_name="Essential_Genes", engine="openpyxl")
        for g_id in e_genes_frame["bnumber"].tolist():
            if g_id not in protected_genes:
                protected_genes.append(g_id)
                
        # protect cofactor metabolism
        # ATP synthase
        protected_genes.extend([g.id for g in self.model.reactions.ATPS4rpp.genes])
        # Transhydrogenase
        protected_genes.extend([g.id for g in self.model.reactions.NADTRHD.genes])
        # NADH dehydrogenase
        protected_genes.extend([g.id for g in self.model.reactions.NADH16pp.genes])
        # Cytochrome oxidases
        protected_genes.extend([g.id for g in self.model.reactions.CYTBO3_4pp.genes])
        protected_genes.extend([g.id for g in self.model.reactions.CYTBDpp.genes])
        protected_genes.extend([g.id for g in self.model.reactions.CYTBD2pp.genes])
        
        # protect specific genes or pathways
        # oxidative PPP
        protected_genes.extend(["b0767", "b1852", "b2029"]) # pgl, zwf, gnd
        # Entner Doudoroff
        protected_genes.extend(["b1851", "b1850"]) # edd, eda
        
        # Pyruvate Metabolism
        protected_genes.extend([g.id for g in self.model.reactions.POX.genes]) # pox
        protected_genes.extend([g.id for g in self.model.reactions.PDH.genes])
        
        # fatty acid biosynthesis
        protected_genes.append("b2836") # aas
        

        protected_genes = list(np.unique(protected_genes))
        
        return protected_genes
    
        
    def _determine_gene_essentiality(self, g_id):
        non_essential_genes = []
        # for g_id in gene_id:
        with self.model as model:
            # knockout gene
            model.genes.get_by_id(g_id).knock_out()
            # solve model
            valid_fun, sol_dict = self._compute_metabolic_functionalities(model)
            # evaluate solution
            if valid_fun:
                non_essential_genes.append({"id": g_id, "type": "m_gene"})
            
        return non_essential_genes
    
    
    
    def _init_general_parameters(self):
        """Specify general, internal parameters for setting up a DEAP genetic
        algorithm and fitness function evaluations

        """
        
        self.minimum_growth_rate = 0.1
        
        self.minimum_relative_rate = 0.95
    
