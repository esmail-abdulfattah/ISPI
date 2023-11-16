#!/bin/bash

# Directory containing the files
DIR="."

# Iterate over files in the specified directory with the naming pattern Q_num_i_dim_*.bin
for file in $DIR/Q_sympack_num_*_dim_*.bin; do
  
  # Extract the file name from the full path
  fileName=$(basename $file)
  
  # Run the './convert' command with the file name as a parameter
  ./convert $fileName
  
  # Print the command being executed
  echo "Executing: ./convert $fileName"
  echo

done

