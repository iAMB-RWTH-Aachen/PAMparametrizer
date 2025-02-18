import gurobipy
from optlang.gurobi_interface import Model, Configuration
import os
from tempfile import NamedTemporaryFile

class CustomGurobiModel(Model):
    """Custom Gurobi Model that overrides the optlang deserialization method."""
    def __init__(self, problem=None, *args, **kwargs):
        # Ensure that the problem is passed and initialized properly
        if problem is not None:
            super().__init__(problem=problem, *args, **kwargs)
            self.problem = problem
        else:
            super().__init__(*args, **kwargs)


    def __setstate__(self, repr_dict):
        """Override optlang's default deserialization to use MPS instead of LP."""
        # Ensure problem is initialized before proceeding
        if self.problem is None:
            raise ValueError("The problem attribute must be initialized before deserialization.")

        # Save the LP data in a temporary file
        with NamedTemporaryFile(suffix=".mps", delete=False) as tmp_file:
            tmp_file_name = tmp_file.name
            tmp_file.close()

        self.problem.write(tmp_file_name)  # Write the model to an MPS file

        # Read the LP file into Gurobi
        try:
            problem = gurobipy.read(tmp_file_name)
        except gurobipy.GurobiError as e:
            print(f"Error reading LP file with Gurobi: {e}")
            raise

        # Initialize the problem in optlang's Gurobi interface
        super().__init__(problem=problem)

        # Restore the model's configuration
        self.configuration = repr_dict["config"]
        self.configuration.problem = self
        self.configuration.__setstate__(repr_dict["config"])

        # Clean up the temporary LP file
        if os.path.exists(tmp_file_name):
            os.remove(tmp_file_name)


    @property
    def interface(self):
        """Ensure the model correctly identifies the Gurobi interface."""
        import optlang.gurobi_interface
        return optlang.gurobi_interface

if __name__ == '__main__':
    from cobra.io import read_sbml_model

    # Load the SBML model
    model = read_sbml_model("Models/iCGB21FR_annotated.xml")


    # Replace the solver with our custom version
    model._solver = CustomGurobiModel(problem=model.solver.problem)
    print(type(model.solver))

    # Now copy the model safely
    model_copy = model.copy()

