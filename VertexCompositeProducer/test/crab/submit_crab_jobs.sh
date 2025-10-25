#!/bin/bash

# Ask user for start and end numbers
read -p "Enter start number: " START
read -p "Enter end number: " END

# Validate that START and END are numbers
if ! [[ "$START" =~ ^[0-9]+$ ]] || ! [[ "$END" =~ ^[0-9]+$ ]]; then
    echo "❌ Error: Start and end must be integers."
    exit 1
fi

# Loop from START to END
for i in $(seq "$START" "$END"); do
    FILE="D0_data_submit_parent_Prime_${i}.py"

    # Check if file exists
    if [[ -f "$FILE" ]]; then
        echo "Submitting $FILE..."
        crab submit "$FILE"
    else
        echo "⚠️  File $FILE not found, skipping."
    fi
done

echo "🎉 Done submitting CRAB jobs!"

