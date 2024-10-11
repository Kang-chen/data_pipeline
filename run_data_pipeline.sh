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

# Step 2: Activate conda environment
if [ "$start_step" -le 2 ] || { [ "$start_step" -ge 3 ] && [ "$start_step" -le 6 ]; }; then
    echo "Activating conda environment 'lamindb'..."
    conda activate lamindb
    check_error 2 "Conda environment activation failed."
fi

# Step 3: Run DataExtract.py
if [ "$start_step" -le 3 ]; then
    echo "Running Python script DataExtract.py..."
    python $data_pipeline/python/dataextract.py --source_id $source_id
    check_error 3 "Python script dataextract.py execution failed."
fi

# Step 4: Run Python script QualityControl.py
if [ "$start_step" -le 4 ]; then
    echo "Running Python script qualitycontrol.py..."
    python $data_pipeline/python/qualitycontrol.py --source_id $source_id
    check_error 4 "Python script qualitycontrol.py execution failed."
fi

# Step 5: Log in to lamin
if [ "$start_step" -le 5 ] || { [ "$start_step" -ge 5 ] && [ "$start_step" -le 6 ]; }; then
    echo "Logging in to lamin..."
    lamin login gefujing98@gmail.com --key yE_MyufO0Uud4cWwMV1TwwDfS63ROOo4OltuY_Rm
    lamin init --storage s3://cartabio/ai/data/fujing_test2 --schema bionty
    check_error 5 "Lamin login failed."
fi

# Step 6: Run Python script UploadLamindb.py with up to 3 retries
if [ "$start_step" -le 6 ]; then
    max_retries=3
    attempt=0
    success=false

    while [ $attempt -lt $max_retries ]; do
        echo "Running Python script uploadlamindb.py (Attempt $((attempt+1))/$max_retries)..."
        python $data_pipeline/python/uploadlamindb.py --source_id $source_id
        if [ $? -eq 0 ]; then
            success=true
            break
        else
            echo "Error in step 6: Python script uploadlamindb.py execution failed. Retrying..."
            attempt=$((attempt+1))
        fi
    done

    if ! $success; then
        echo "Error in step 6: Python script uploadlamindb.py failed after $max_retries attempts."
        exit 1
    fi
fi

echo "Pipeline completed successfully."

