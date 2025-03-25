#!/bin/bash

# Define the output file
OUTPUT_FILE="LO_ref.dat"

# Initialize the output file
> $OUTPUT_FILE # Clear the output file or create it if it doesn't exist

# Specify directories explicitly
DIRECTORIES=("run_21" "run_22" "run_23" "run_24" "run_25" "run_26" "run_27") # Add or modify as needed

# Loop through specified directories
for dir in "${DIRECTORIES[@]}"; do
    if [[ -d "$dir" ]]; then
        # Construct the path to result_summary.dat
        SUMMARY_FILE="$dir/summary/result_summary.dat"

        if [[ -f "$SUMMARY_FILE" ]]; then
            # Extract center-of-mass energy
            COM_ENERGY=$(grep "p p --> H" "$SUMMARY_FILE" | grep -oP "\d+(?= TeV)")
            # Convert to GeV
            COM_ENERGY=$((COM_ENERGY * 1000))

            # Extract LO result
            LO_RESULT=$(grep "LO:" "$SUMMARY_FILE" | awk '{print $2}' | sed 's/^LO://' )
            if [[ -z "${LO_RESULT//[[:space:]]/}" || "$LO_RESULT" == "LO:" ]]; then
                LO_RESULT=$(grep "LO:" "$SUMMARY_FILE" | awk '{print $3}')
            fi

            # Append to the output file
            echo -e "${COM_ENERGY}\t ${LO_RESULT}\t $dir" >> "$OUTPUT_FILE"
        else
            echo "File $SUMMARY_FILE does not exist. Skipping..."
        fi
    else
        echo "Directory $dir does not exist. Skipping..."
    fi
done

echo "Data has been extracted and saved to $OUTPUT_FILE."

