#!/bin/bash

# Define the output file
OUTPUT_FILE="LO_ref.dat"
OUTPUT_FILE1="NLO_ref.dat"

# Initialize the output file
> $OUTPUT_FILE1 # Clear the output file or create it if it doesn't exist

# Create multiple temporary files
TEMP_FILE1=$(mktemp)
TEMP_FILE2=$(mktemp)
TEMP_FILE3=$(mktemp)
TEMP_FILE4=$(mktemp)
 
# Specify directories explicitly
#DIRECTORIES=("run_28")
DIRECTORIES=("run_13000" "run_12000" "run_11000" "run_10000" "run_9000" "run_8000" "run_7000" "run_6000" "run_5000" "run_4000" "run_3000" "run_2000" "run_1000" "run_500")


# Loop through specified directories
for dir in "${DIRECTORIES[@]}"; do
    if [[ -d "$dir" ]]; then
        # Construct the path to result_summary.dat
        SUMMARY_FILE="$dir/ihixs_output"

        if [[ -f "$SUMMARY_FILE" ]]; then

            # Extract NLO result (fifth and sixth columns from the 49th line)
	     LO_RESULT=$(sed -n '14p' "$SUMMARY_FILE" | awk '{print $3}')
	    NLO_RESULT=$(sed -n '15p' "$SUMMARY_FILE" | awk '{print $3}')
	    
            COM_ENERGY=$(sed -n '4p' "$SUMMARY_FILE" | awk '{print $3}')
           
            # Write data to the temporary files
            echo -e "${COM_ENERGY}\t ${NLO_RESULT}\t    " >> "$TEMP_FILE1"
            echo -e "${COM_ENERGY}\t ${LO_RESULT}\t    " >> "$TEMP_FILE2"
            
        else
            echo "File $SUMMARY_FILE does not exist. Skipping..."
        fi
    else
        echo "Directory $dir does not exist. Skipping..."
    fi
done
#	cat $TEMP_FILE1
#	cat $TEMP_FILE2
#	cat $TEMP_FILE3
#	cat $TEMP_FILE4
            # Combine temp files into OUTPUT_FILE1
            echo "NLO_RESULT" >> "$OUTPUT_FILE1"
            cat "$TEMP_FILE1" >> "$OUTPUT_FILE1"
            echo "" >> "$OUTPUT_FILE1"
            echo "LO_RESULT" >> "$OUTPUT_FILE1"
            cat "$TEMP_FILE2" >> "$OUTPUT_FILE1"

            # Clean up temporary files if no longer needed
            rm "$TEMP_FILE1" "$TEMP_FILE2"


#echo "LO Data has been extracted and saved to $OUTPUT_FILE."
#echo "NLO Data has been extracted and saved to $OUTPUT_FILE2."

