from dataclasses import dataclass
import gurobipy as gp
import pandas as pd

@dataclass
class KcatLinearModel:
    substrate_rates: list
    esc_matrix: [list, list]
    right_hand_side: list
    identifiers: [str, str]
    valid_data: pd.DataFrame
    lmbda: float

    @property
    def gurobi_model(self):
        model = gp.Model('kcat_linear_optimization')
        objective_expression = gp.LinExpr()
        #loop over all simulations
        for id in self.identifiers:
            model.addVar(name=id)
        model.update()
        for subst_upt, esc_simulation in zip(self.substrate_rates,self.esc_matrix):
        #add variables and constraints
            for id, esc, rhs in zip(self.identifiers, esc_simulation, self.right_hand_side):
                var = model.getVarByName(name= id)
                model.addConstr(esc*var == rhs, name=f'{id}_{round(subst_upt, 3)}')

                rxn = id.split('_')[0]
                if rxn in self.valid_data.columns:
                    true_value = self.valid_data[rxn][self.valid_data['R1_ub']==subst_upt].iloc[0]
                    #make objective expression absolute by adding to other constraintsm, see https://lpsolve.sourceforge.net/5.1/absolute.htm
                    abs_var = self.create_absolute_variable(model,gp.LinExpr(var), f'{id}_{round(subst_upt, 3)}_var')
                    abs_error = self.create_absolute_variable(model, gp.LinExpr(esc*var - true_value), id+f'_error_{round(subst_upt, 3)}')
                    #lasso regression term
                    objective_expression.addTerms(self.lmbda, abs_var)
                    # objective_expression.addTerms(self.lmbda, var)

            #error termo
            objective_expression.addTerms(1.0, abs_error)
            # objective_expression.addTerms(esc,var)
            # objective_expression.addConstant( -rhs)

        model.setObjective(objective_expression, gp.GRB.MINIMIZE)
        model.update()
        for const in model.getConstrs():
            print(const.ConstrName,const.Sense, const.RHS)
        print(objective_expression)
        return model

    def create_absolute_variable(self, model:gp.Model, expression: gp.LinExpr, name:str)-> gp.Var:
        # make objective expression absolute by adding to other constraintsm, see https://lpsolve.sourceforge.net/5.1/absolute.htm
        abs_var = model.addVar(name = name)
        model.addConstr(expression <= abs_var, name= name+'_pos')
        neg_expression = self.invert_linear_expression(expression)
        model.addConstr(expression >= -abs_var, name= name+'_neg')
        return abs_var

    @staticmethod
    def invert_linear_expression(expression:gp.LinExpr):
        #invert coefficients
        new_expression = gp.LinExpr()
        for i in range(expression.size()):
            var = expression.getVar(i)
            new_coeff = -expression.getCoeff(i)
            new_expression.addTerms(new_coeff, var)

            try:
                var = expression.getVar(i)
                new_coeff = -expression.getCoeff(i)
                new_expression.addTerms(new_coeff, var)
            except:
                continue
        #invert constants
        new_constant = -expression.getConstant()
        new_expression.addConstant(new_constant)
        return new_expression
