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
DIRECTORIES=("run_28" "run_29" "run_30" "run_31" "run_32" "run_33" "run_34" "run_35" "run_36" "run_37" "run_38" "run_39" "run_40" "run_41") # Add or modify as needed

# Loop through specified directories
for dir in "${DIRECTORIES[@]}"; do
    if [[ -d "$dir" ]]; then
        # Construct the path to result_summary.dat
        SUMMARY_FILE="$dir/result/result..MATRIX.NLO.result/fixed-scale/complete/NLO.CS.31/overview.NLO.CS.31.dat"
        PARAMETER_FILE="$dir/file_parameter.dat"

        if [[ -f "$SUMMARY_FILE" ]]; then

            # Extract NLO result (fifth and sixth columns from the 49th line)
	    NLO_RESULT=$(sed -n '49p' "$SUMMARY_FILE" | awk '{print $5, $7}')
	    
	    # Extract VA result (fifth and sixth columns from the 49th line)
	    VA_RESULT=$(sed -n '51p' "$SUMMARY_FILE" | awk '{print $5, $7}')
	    
	    # Extract CA result (fifth and sixth columns from the 49th line)
	    CA_RESULT=$(sed -n '52p' "$SUMMARY_FILE" | awk '{print $5, $7}')
	    
	    # Extract RA result (fifth and sixth columns from the 49th line)
	    RA_RESULT=$(sed -n '53p' "$SUMMARY_FILE" | awk '{print $5, $7}')

            COM_ENERGY=$(sed -n '2p' "$PARAMETER_FILE" | awk '{print $3}')
            COM_ENERGY=$(echo "$COM_ENERGY * 2" | bc)

           
            # Write data to the temporary files
            echo -e "${COM_ENERGY}\t ${NLO_RESULT}\t    " >> "$TEMP_FILE1"
            echo -e "${COM_ENERGY}\t ${VA_RESULT}\t    " >> "$TEMP_FILE2"
            echo -e "${COM_ENERGY}\t ${CA_RESULT}\t    " >> "$TEMP_FILE3"
            echo -e "${COM_ENERGY}\t ${RA_RESULT}\t    " >> "$TEMP_FILE4"
            
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
            echo "VA_RESULT" >> "$OUTPUT_FILE1"
            cat "$TEMP_FILE2" >> "$OUTPUT_FILE1"
            echo "" >> "$OUTPUT_FILE1"
            echo "CA_RESULT" >> "$OUTPUT_FILE1"
            cat "$TEMP_FILE3" >> "$OUTPUT_FILE1"
            echo "" >> "$OUTPUT_FILE1"
            echo "RA_RESULT" >> "$OUTPUT_FILE1"
            cat "$TEMP_FILE4" >> "$OUTPUT_FILE1"

            # Clean up temporary files if no longer needed
            rm "$TEMP_FILE1" "$TEMP_FILE2" "$TEMP_FILE3"


#echo "LO Data has been extracted and saved to $OUTPUT_FILE."
#echo "NLO Data has been extracted and saved to $OUTPUT_FILE2."

