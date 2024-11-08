#!/bin/bash

PASSWORD="thisismynewpassword"
REMOTE_USER="rbanderson"
REMOTE_HOST="jakar.utep.edu"
REMOTE_DIR="/home/rbanderson/projectfiles"
LOCAL_DIR="shusolutionoutput"

# Create the local directory if it doesn't exist
mkdir -p $LOCAL_DIR

# List of Ainput values
Ainput_values=(
    "2.1" "2.01" "2.001" "2.0001" "2.00001" 
    "2.000001" "2.0000001" "2.00000001" "2.000000001" 
    "2.0000000001" "2.00000000001" "2.000000000001"
)

# Loop through each Ainput value and download the corresponding .h5 file
for Ainput in "${Ainput_values[@]}"; do
    FILENAME="shusolutionoutput${Ainput}.h5"
    REMOTE_FILE="${REMOTE_DIR}/${FILENAME}"
    LOCAL_FILE="${LOCAL_DIR}/${FILENAME}"
    
    sshpass -p "$PASSWORD" scp ${REMOTE_USER}@${REMOTE_HOST}:${REMOTE_FILE} ${LOCAL_FILE}
done

echo "Download completed."

