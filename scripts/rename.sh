#!/bin/sh

# Go to the directory containing the files
cd .

# The starting index you want
start_index=41

for file in image.*.png; do
    current_index="${file#image.}"     # Remove "image." prefix
    current_index="${current_index%.png}"  # Remove ".png" suffix
    new_index=$(echo "$current_index" | awk '{ printf "%d", $1 }')  # Convert to decimal
    new_index=$((new_index + start_index))
    new_name="image.$(printf "%04d" $new_index).png"
    echo $new_name
    mv "$file" "$new_name"
done
