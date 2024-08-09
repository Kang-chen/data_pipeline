#!/bin/bash
source_id=$1
# Specify the SharePoint site and folder URL
web_url="https://netorgft12480156.sharepoint.com/sites/cartabio.ai-AIUnitDataPipeline"
folder_url="https://netorgft12480156.sharepoint.com/sites/cartabio.ai-AIUnitDataPipeline/Shared%20Documents/%5BAI%20Unit%5D%20%F0%9F%93%8A%20Data%20Pipeline%20%F0%9F%93%8A"

# List files and extract the download link
download_url=$(m365 file list --webUrl "$web_url" --folderUrl "$folder_url" | jq -r '.[] | select(.name=="Data_collection.xlsx") | .["@microsoft.graph.downloadUrl"]')

# Check if the download link was successfully extracted
if [ -z "$download_url" ]; then
  echo "The specified file was not found."
  exit 1
fi

# Set the output file name
output_file="Data_collection.xlsx"

# Download the file
curl -L -o "$output_file" "$download_url"

# Notify the user of the download completion
echo "The file has been downloaded and saved as $output_file"

echo "----------------start Download script for $source_id------------------"
# Rscript ./convertAnn.R $source_id