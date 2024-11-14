#!/bin/bash

# Read the input file composition.txt
input="composition.txt"

# Initialize a variable to store the number of words in a line
i=0

# Loop through each line of the file and count the number of words
while IFS= read -r line
do
    i=$(echo $line | wc -w)
done < "$input"

# Write the number of words to num_Lines.txt
echo $i > num_Lines.txt

# Use awk to extract each word from the line and store it in the holder array
declare -a holder
for ((j=1; j<=$i; j++))
do
    holder[$j]=$(awk '{print $'$j'}' composition.txt)
done

# Write the extracted words back into composition.txt
for ((j=1; j<=$i; j++))
do
    if [ $j == 1 ]; then
        # Overwrite the file with the first word
        echo ${holder[$j]} > composition.txt
    else
        # Append the remaining words to the file
        echo ${holder[$j]} >> composition.txt
    fi
done

# Exit the script
exit 0
