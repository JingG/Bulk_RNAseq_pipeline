#!/bin/bash

#### combine the outcomes of the featureCounts results into one count matrix 

# Check if two arguments are provided
if [ $# -ne 2 ]; then
    echo "Usage: $0 <column_number> <file_containing_list_of_input_files>"
    exit 1
fi

# The first argument is the column number
column_number=$1
# The second argument is the file containing the list of input files
file_list=$2

# Check if the file list exists
if [ ! -f "$file_list" ]; then
    echo "Error: File list $file_list does not exist."
    exit 1
fi

# Create a temporary directory
temp_dir=$(mktemp -d)

outfile="combined_counts.txt"


# Extract the specified column from each file listed in file_list
while IFS= read -r file; do
    if [ -f "$file" ]; then
        # Extract the column and save to a temporary file
        cut -f$column_number "$file" > "$temp_dir/temp_col"
        
        # Create a new file with the filename as the first line, followed by the column data
	echo "$(basename "$file" .count)" > "$temp_dir/$(basename "$file").col"
        tail -n +3 "$temp_dir/temp_col" >> "$temp_dir/$(basename "$file").col"
        
        rm "$temp_dir/temp_col"
    else
        echo "Warning: File $file does not exist. Skipping."
    fi
done < "$file_list"

# Use paste to combine all extracted columns
paste "$temp_dir"/*.col > count.txt

# extract the ID 
echo "ID" > gene_ID.txt
## one sample of the featureCounts result
cut -f 1 4-2.FC_trim/FiGS0078.count |tail -n +3  >>gene_ID.txt

# Use paste to combine all extracted columns
paste gene_ID.txt count.txt > $outfile

# Clean up temporary files
rm -r "$temp_dir" count.txt

echo "Combined output saved to $outfile."



