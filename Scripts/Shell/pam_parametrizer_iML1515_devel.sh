#!/usr/bin/zsh

##SBATCH section
#SBATCH --time=00:45:00 --cpus-per-task=8 --nodes=1 --job-name=pamparamertizer
#SBATCH --partition=devel

##Loading the required modules and packages
module purge
module load GCCcore/.12.2.0 #Required to load Gurobi
module load Python/3.10.8
module load Gurobi/10.0.2 #Required for the Gurobi license bindings
export CONDA_ROOT=$HOME/miniconda3
source $CONDA_ROOT/etc/profile.d/conda.sh
export PATH="$CONDA_ROOT/bin:$PATH"
conda activate PAMparametrizer

## Run the script and save the results
## Run the script and save the results
cd $HOME/Programs/PAM_parametrization

# Run Python script in parallel
echo "Starting script..."

# Define function to run Python script with specific configuration
run_script() {
    python3 -m Scripts.Testing.pam_parametrizer_performance_analysis --hyper_processes 2 --iterations 5 --configuration $1 > Results/pam_parameterizer_performance_analysis_$1_240502.txt
}

# Run the Python script in parallel for each configuration
for config in all before False; do
    run_script $config &
done

wait