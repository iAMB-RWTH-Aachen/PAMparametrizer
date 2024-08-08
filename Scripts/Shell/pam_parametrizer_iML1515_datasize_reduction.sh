#!/usr/bin/zsh

##SBATCH section
#SBATCH --time=72:00:00 --cpus-per-task=4 --nodes=1 --job-name=pamparamertizer_datareduction
##SBATCH --output=pamparametrizer_datareduction_240808.txt
##SBATCH --mail-type=END,FAIL --mail-user=samira.vandenbogaard@rwth-aachen.de

##Loading the required modules and packages
module purge
module load GCCcore/.12.2.0 #Required to load Gurobi and Python
module load Python/3.10.8
module load Gurobi/10.0.2 #Required for the Gurobi license bindings
export CONDA_ROOT=$HOME/miniconda3
source $CONDA_ROOT/etc/profile.d/conda.sh
export PATH="$CONDA_ROOT/bin:$PATH"
conda activate PAMparametrizer

## Run the script and save the results
cd $HOME/Programs/PAM_parametrization

# Run Python script in parallel
echo "Starting script..."

# Define function to run Python script with specific configuration
run_script() {
    python3 -m Scripts.Scripts.Data_requirements.global_datasize_reduction > Results/pam_parameterizer_datasize_reduction_240808.txt
}

# Run the Python script in parallel for each configuration
run_script

wait

echo "Script completed."
