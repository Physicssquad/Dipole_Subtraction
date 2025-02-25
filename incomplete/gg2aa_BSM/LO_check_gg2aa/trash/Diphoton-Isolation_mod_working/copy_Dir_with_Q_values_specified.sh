#!/bin/bash

# Define the main directory and the target Q values
main_dir="Diphoton-invmass"
#q_values=(400 700 1000 1300 1600 1900)
q_values=(100 200 300 400 500 600 700 800 900 1000 1100 1200 1300 1400 1500 1600 1700 1800 1900)

# Loop through each Q value, create a new directory, and modify aqstep in run.order.dat
for Q in "${q_values[@]}"; do
    # Create the new directory name with the Q value
    new_dir="${main_dir}_Q_${Q}"
    
    # Copy the main directory to the new directory
    cp -r "$main_dir" "$new_dir"
    
    # Update the aqstep value on line 11 in the copied run.order.dat file
    sed -i "11s/.*/${Q}d0                  ! aqstep; step width of diff. distribution/" "$new_dir/run.order.dat"
    
    echo "Created directory $new_dir with aqstep set to ${Q}d0 in run.order.dat"
done

