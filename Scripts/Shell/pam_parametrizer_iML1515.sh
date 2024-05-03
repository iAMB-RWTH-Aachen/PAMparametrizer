#!/usr/bin/zsh

##SBATCH section
#SBATCH --time=05:00:00 --cpus-per-task=10 --nodes=1 --job-name=pamparamertizer

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
python3 -m Scripts.Testing.pam_parametrizer_performance_analysis --hyper_processes 3 --configuration 'all'> Results/pam_parameterizer_performance_analysis_all_240502.txt
python3 -m Scripts.Testing.pam_parametrizer_performance_analysis --hyper_processes 3 --configuration 'before'> Results/pam_parameterizer_performance_analysis_before_240502.txt
python3 -m Scripts.Testing.pam_parametrizer_performance_analysis --hyper_processes 3 --configuration 'False'> Results/pam_parameterizer_performance_analysis_false_240502.txt