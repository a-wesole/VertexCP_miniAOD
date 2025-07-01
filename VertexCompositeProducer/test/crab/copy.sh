#!/bin/bash

# List of values to replace PT-0
#values=("PT-10" "PT-60" "PT-30" "PT-20" "PT-40" "PT-80" "PT-8" "PT-1" "PT-0")
#values=("PT-60" "PT-30" "PT-20" "PT-40" "PT-80" "PT-8")
values=("PT-1")

# Loop over each value in the list
for value in "${values[@]}"; do
  echo "_________________"
  echo "Processing $value"
  echo "_________________"

  # Change directory to the source location
  cd /eos/purdue/store/user/wxie/D0_mc_TTrees

  # Perform the SCP command with the current value
  scp $(find "$value" -type f -name '*.root') $bell:/scratch/bell/awesole/D0MC_Run3_Trees/"$value"

  # Return to the previous directory
  cd -
done


