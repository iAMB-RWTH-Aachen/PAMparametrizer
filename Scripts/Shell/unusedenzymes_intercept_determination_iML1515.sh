#!/usr/bin/zsh

##SBATCH section
#SBATCH --partition=c23ms
#SBATCH --time=24:00:00 --cpus-per-task=3 --nodes=1
#SBATCH --job-name=UE_determination_iML1515 --output=UE_intercept_iML1515_242609.txt
#SBATCH --mail-type=END,FAIL --mail-user=samira.vandenbogaard@rwth-aachen.de


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
cd $HOME/Programs/PAM_Parametrization

# Run Python script in parallel
echo "Starting script..."

python3 -m Scripts.Protein_sectors.1_scan_maxgrowthrate_UE

wait

echo "Script completed."
