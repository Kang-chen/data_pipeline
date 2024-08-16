#!/bin/bash

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

data_pipeline="$script_dir"

echo "The project directory is located at: $data_pipeline"

# Check if at least one argument (source_id) is provided
if [ "$#" -lt 1 ]; then
    echo "Usage: $0 <source_id> [start_step]"
    exit 1
fi

source_id=$1
start_step=${2:-1}

check_error() {
    if [ $? -ne 0 ]; then
        echo "Error in step $1: $2"
        exit 1
    fi
}

# Step 1: Execute download_xlsx.sh
if [ "$start_step" -le 1 ]; then
    echo "Running download_xlsx.sh..."
    $data_pipeline/bash/download_xlsx.sh
    check_error 1 "download_xlsx.sh failed."
fi

# Step 2: Run Docker container and execute R script
if [ "$start_step" -le 2 ]; then
    echo "Running Docker container and executing R script..."
    docker run --rm \
		-v $data_pipeline:/home/rstudio/data_pipeline/ \
		-v "$(pwd)":/home/rstudio/run \
		-w /home/rstudio/run \
		kang/scfetch:1.61 \
		Rscript /home/rstudio/data_pipeline/R/convert2Ann.R $source_id
    check_error 2 "Docker container or R script execution failed."
fi

# Step 3: Activate conda environment
if [ "$start_step" -le 3 ] || { [ "$start_step" -ge 5 ] && [ "$start_step" -le 6 ]; }; then
    echo "Activating conda environment 'lamindb'..."
    conda activate lamindb
    check_error 3 "Conda environment activation failed."
fi

# Step 4: Log in to lamin
if [ "$start_step" -le 4 ] || { [ "$start_step" -ge 5 ] && [ "$start_step" -le 6 ]; }; then
    echo "Logging in to lamin..."
    lamin login fujingge
    check_error 4 "Lamin login failed."
fi

# Step 5: Run Python script 1-qc.py
if [ "$start_step" -le 5 ]; then
    echo "Running Python script 1-qc.py..."
    python $data_pipeline/python/1-qc.py --source_id $source_id
    check_error 5 "Python script 1-qc.py execution failed."
fi

# Step 6: Run Python script 2-lamindb-aws.py with up to 3 retries
if [ "$start_step" -le 6 ]; then
    max_retries=3
    attempt=0
    success=false

    while [ $attempt -lt $max_retries ]; do
        echo "Running Python script 2-lamindb-aws.py (Attempt $((attempt+1))/$max_retries)..."
        python $data_pipeline/python/2-lamindb-refactor.py --source_id $source_id
        if [ $? -eq 0 ]; then
            success=true
            break
        else
            echo "Error in step 6: Python script 2-lamindb-aws.py execution failed. Retrying..."
            attempt=$((attempt+1))
        fi
    done

    if ! $success; then
        echo "Error in step 6: Python script 2-lamindb-aws.py failed after $max_retries attempts."
        exit 1
    fi
fi

echo "Pipeline completed successfully."

