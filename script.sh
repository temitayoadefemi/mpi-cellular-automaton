#!/bin/bash
# run_scaling.sh
# This script performs weak and strong scaling tests for Conway's Game of Life automaton.
# It runs each configuration 10 times and captures only the timing data (using the time command).
# The automaton’s own output is discarded.

OUTPUT_FILE="timing_results.txt"
SEED=42  # Fixed seed for all runs

# Remove any existing output file.
rm -f "$OUTPUT_FILE"

# Write a minimal header.
echo "Timing results for Conway's Game of Life scaling tests" > "$OUTPUT_FILE"
echo "Start Time: $(date)" >> "$OUTPUT_FILE"
echo "" >> "$OUTPUT_FILE"

#################################
# Weak Scaling Experiment
#################################
echo "Weak Scaling Test: Fixed 5 processors, varying landscape sizes" >> "$OUTPUT_FILE"
# List of landscape sizes.
landscapes=(500 1000 1500 2000 2500 3000 3500 4000 4500 5000)

for LAND in "${landscapes[@]}"; do
    echo "Landscape Size: $LAND" >> "$OUTPUT_FILE"
    for run in {1..10}; do
        echo "  Run #$run:" >> "$OUTPUT_FILE"
        # Run the program with 5 processors, discard its output, and capture the timing data.
        { time mpirun -n 5 ./automaton "$SEED" -landscape "$LAND" > /dev/null; } 2>> "$OUTPUT_FILE"
        echo "  ------------------------------------" >> "$OUTPUT_FILE"
    done
    echo "" >> "$OUTPUT_FILE"
done

#################################
# Strong Scaling Experiment
#################################
echo "Strong Scaling Test: Fixed landscape size (5000), varying processors" >> "$OUTPUT_FILE"
# List of processor counts.
procs=(2 4 6 8 10 12 14 16)

for NP in "${procs[@]}"; do
    echo "Processor Count: $NP" >> "$OUTPUT_FILE"
    for run in {1..10}; do
        echo "  Run #$run:" >> "$OUTPUT_FILE"
        { time mpirun -n "$NP" ./automaton "$SEED" -landscape 5000 > /dev/null; } 2>> "$OUTPUT_FILE"
        echo "  ------------------------------------" >> "$OUTPUT_FILE"
    done
    echo "" >> "$OUTPUT_FILE"
done

echo "End Time: $(date)" >> "$OUTPUT_FILE"
