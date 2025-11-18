#!/bin/bash

# Set error handling
set -e
set -u

# Check if the conda command is available
if [ ! -x "$(command -v conda)" ]; then
    echo "Error: conda is not installed or not in the PATH." >&2
    exit 1
fi

# Specify the conda environment to run the script in (assuming the virtual environment name is 'docking')
CONDA_ENV_NAME="docking"

# Use conda run to execute the script in the specified environment
run_in_conda_env() {
    conda run -n "$CONDA_ENV_NAME" "$@"
}

# Check if the Python version in the specified conda environment is >= 3.0
check_python_version() {
    local env_name="$1"
    local version
    version=$(conda run -n "$env_name" python --version 2>&1 | awk '{print $2}')
    min_version="3.0"

    if [[ $(echo -e "$version\n$min_version" | sort -V | head -n1) != "$min_version" ]]; then
		echo "Error: Python version in the $env_name environment must be at least 3.0+. Current version: $version" >&2
		echo "Please make sure you have installed Python 3 in the '$env_name' environment."
		echo "You can check the Python version in the environment using 'python -V'."
		echo "If you have previously configured Python 3 in the 'docking' environment, please ensure it is activated."
		echo "To activate the 'docking' environment, run 'conda deactivate' to exit all conda environments,"
		echo "then run 'conda activate docking'."
		exit 1
    fi
} 

# Check the Python version in the 'docking' environment
check_python_version "$CONDA_ENV_NAME"

# Check the path ./data_input/data_ligand/partition_1
ligand_path="./data_input/data_ligand/partition_1"
if [ ! -d "$ligand_path" ]; then
    echo "Path $ligand_path does not exist, please create the path."
    exit 1
fi

ligand_contents=$(ls -A "$ligand_path")
if [ -z "$ligand_contents" ]; then
    echo "No ligand detected in the path $ligand_path, please add ligand to this path."
    exit 1
else
    echo "Files in $ligand_path:"
    for file in $ligand_contents; do
        echo "$file"
		echo ""
    done
fi

# Check the path ./data_input/data_protein/partition_1
protein_path="./data_input/data_protein/partition_1"
if [ ! -d "$protein_path" ]; then
    echo "Path $protein_path does not exist, please create the path."
	echo ""
    exit 1
fi

protein_contents=$(ls -A "$protein_path")
if [ -z "$protein_contents" ]; then
    echo "No protein detected in the path $protein_path, please add protein to this path."
    exit 1
fi

# Check other paths
paths=("./data_output/1_PDBQT/ligand/file_1" "./data_output/1_PDBQT/protein/file_1" "./data_output/2_cavity_output/cavity_1" "./data_output/2_cavity_output/file_1" "./data_output/4_docking_results/file_1" "./data_output/5_statistics")

for path in "${paths[@]}"; do
    if [ -e "$path" ]; then
        echo "Files already exist in the path $path, please make sure these are the files you intend to analyze."
    fi
done

echo ""
echo "All paths have passed the checks, and no non-compliant conditions were detected."
echo ""
echo "Starting..."

##################################################################
##################################################################
# Execute Python scripts one by one in the 'docking' environment #
##################################################################
##################################################################

echo ""
echo ""
echo "######               PDBQT Processing              ######"
run_in_conda_env python 1_PDBQT_process.py -mol ligand -num 1 --timeout 600 >> reserve_docking.log 2>&1 
run_in_conda_env python 1_PDBQT_process.py -mol protein -num 1 --timeout 600 >> reserve_docking.log 2>&1
echo "PDBQT processing completed"

echo ""
echo ""
echo "######               Cavity Detection              ######"
run_in_conda_env python 2_cavity_detected.py -num 1 --timeout 2000 >> reserve_docking.log 2>&1
echo "Cavity detection completed"

echo ""
echo ""
echo "######        Pocket Parameters Calculation        ######"
##To utilize the ligand located at './data_output/1_PDBQT/ligand/file_1', please uncomment the 112th line. The default location is './data_input/data_protein/partition_1'.
#ligand_path="./data_output/1_PDBQT/ligand/file_1"
ligand_contents=$(ls -A "$ligand_path")
for file in $ligand_contents; do
    filename_no_extension="${file%.*}"
    run_in_conda_env python 3_pocket_parameters.py -num 1 -ligand "$filename_no_extension" >> reserve_docking.log 2>&1
done 
echo "Pocket parameter calculation completed"

echo ""
echo ""
echo "######               Reserve Docking               ######"
run_in_conda_env python 4_reserve_docking.py -num 1 --timeout 2000 >> reserve_docking.log 2>&1
echo "Reserve docking completed"

echo ""
echo ""
echo "######       Statistical Results Calculation       ######"
run_in_conda_env python 5_stastics_results.py >> reserve_docking.log 2>&1
echo "Statistical results calculation completed"