#!/bin/bash
#SBATCH --partition=batch
#SBATCH --job-name=eukaryomedl
#SBATCH --ntasks=1
#SBATCH --time=2:00:00
#SBATCH --mem=50G
#SBATCH --output=./slurmoutputs/eukaryomedl.out
#SBATCH --error=./slurmoutputs/eukaryomedl.out
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16

# Make directories to hold databases
# -p means the function only runs if the database isn't already there

mkdir -p eukaryome
mkdir -p eukaryome/lsu
mkdir -p eukaryome/ssu
mkdir -p eukaryome/its
mkdir -p eukaryome/longread

#---------------------------- LSU DATABASE ------------------------------------#

# Download the file
wget https://sisu.ut.ee/wp-content/uploads/sites/643/General_EUK_longread_v1.9.3.zip

# Unzip the file
unzip General_EUK_longread_v1.9.3.zip

# Move the downloaded file to the correct place
mv General_EUK_longread_v1.9.3.fasta eukaryome/longread/

rm General_EUK_longread_v1.9.3.zip

# Define input and output files
input_fasta="eukaryome/longread/General_EUK_longread_v1.9.3.fasta"
output_fasta_fungi="eukaryome/longread/General_EUK_longread_v1.9.3_fungi.fasta"
output_fasta_glomeromycota="eukaryome/longread/General_EUK_longread_v1.9.3_glomeromycota.fasta"

# Filter FASTA for entries with k__Fungi
awk '/^>/ {if ($0 ~ /k__Fungi/) print_header=1; else print_header=0}
     print_header {print}
     !/^>/ {if (print_header) print}' $input_fasta > $output_fasta_fungi

# Filter FASTA for entries with p__Glomeromycota
awk '/^>/ {if ($0 ~ /p__Glomeromycota/) print_header=1; else print_header=0}
      print_header {print}
      !/^>/ {if (print_header) print}' $input_fasta > $output_fasta_glomeromycota

# Now we'll count how many entries are in each version of the database

# Output file for the table
output_table="eukaryome/eukaryome_LSU_sequencecounts.txt"

# Header for the table
echo -e "File\tSequence_Count" > $output_table

# Loop through each specified FASTA file
for file in ${input_fasta} ${output_fasta_fungi} ${output_fasta_glomeromycota}; do
    # Count sequences in the file
    count=$(grep -c "^>" "$file")

    # Add the filename and count to the table
    echo -e "$(basename $file)\t$count" >> $output_table
done

# Display the table
cat $output_table

#---------------------------- LSU DATABASE ------------------------------------#

# Download the file
wget https://sisu.ut.ee/wp-content/uploads/sites/643/General_EUK_LSU_v1.9.3.zip

# Unzip the file
unzip General_EUK_LSU_v1.9.3.zip

# Move the downloaded file to the correct place
mv General_EUK_LSU_v1.9.3.fasta eukaryome/lsu/

rm General_EUK_LSU_v1.9.3.zip

# Define input and output files
input_fasta="eukaryome/lsu/General_EUK_LSU_v1.9.3.fasta"
output_fasta_fungi="eukaryome/lsu/General_EUK_LSU_v1.9.3_fungi.fasta"
output_fasta_glomeromycota="eukaryome/lsu/General_EUK_LSU_v1.9.3_glomeromycota.fasta"

# Filter FASTA for entries with k__Fungi
awk '/^>/ {if ($0 ~ /k__Fungi/) print_header=1; else print_header=0}
     print_header {print}
     !/^>/ {if (print_header) print}' $input_fasta > $output_fasta_fungi

# Filter FASTA for entries with p__Glomeromycota
awk '/^>/ {if ($0 ~ /p__Glomeromycota/) print_header=1; else print_header=0}
      print_header {print}
      !/^>/ {if (print_header) print}' $input_fasta > $output_fasta_glomeromycota

# Now we'll count how many entries are in each version of the database

# Output file for the table
output_table="eukaryome/eukaryome_LSU_sequencecounts.txt"

# Header for the table
echo -e "File\tSequence_Count" > $output_table

# Loop through each specified FASTA file
for file in ${input_fasta} ${output_fasta_fungi} ${output_fasta_glomeromycota}; do
    # Count sequences in the file
    count=$(grep -c "^>" "$file")

    # Add the filename and count to the table
    echo -e "$(basename $file)\t$count" >> $output_table
done

# Display the table
cat $output_table

#---------------------------- SSU DATABASE ------------------------------------#

# Download the file
wget https://sisu.ut.ee/wp-content/uploads/sites/643/General_EUK_SSU_v1.9.3.zip

# Unzip the file
unzip General_EUK_SSU_v1.9.3.zip

# Move the downloaded file to the correct place
mv General_EUK_SSU_v1.9.3.fasta eukaryome/ssu/

rm General_EUK_SSU_v1.9.3.zip

# Define input and output files
input_fasta="eukaryome/ssu/General_EUK_SSU_v1.9.3.fasta"
output_fasta_fungi="eukaryome/ssu/General_EUK_SSU_v1.9.3_fungi.fasta"
output_fasta_glomeromycota="eukaryome/ssu/General_EUK_SSU_v1.9.3_glomeromycota.fasta"

# Filter FASTA for entries with k__Fungi
awk '/^>/ {if ($0 ~ /k__Fungi/) print_header=1; else print_header=0}
     print_header {print}
     !/^>/ {if (print_header) print}' $input_fasta > $output_fasta_fungi

# Filter FASTA for entries with p__Glomeromycota
awk '/^>/ {if ($0 ~ /p__Glomeromycota/) print_header=1; else print_header=0}
      print_header {print}
      !/^>/ {if (print_header) print}' $input_fasta > $output_fasta_glomeromycota

# Now we'll count how many entries are in each version of the database

# Output file for the table
output_table="eukaryome/eukaryome_SSU_sequencecounts.txt"

# Header for the table
echo -e "File\tSequence_Count" > $output_table

# Loop through each specified FASTA file
for file in ${input_fasta} ${output_fasta_fungi} ${output_fasta_glomeromycota}; do
    # Count sequences in the file
    count=$(grep -c "^>" "$file")

    # Add the filename and count to the table
    echo -e "$(basename $file)\t$count" >> $output_table
done

# Display the table
cat $output_table

#---------------------------- ITS DATABASE ------------------------------------#

# Download the file
wget https://sisu.ut.ee/wp-content/uploads/sites/643/General_EUK_ITS_v1.9.3.zip

# Unzip the file
unzip General_EUK_ITS_v1.9.3.zip

# Move the downloaded file to the correct place
mv General_EUK_ITS_v1.9.3.fasta eukaryome/its/

rm General_EUK_ITS_v1.9.3.zip

# Define input and output files
input_fasta="eukaryome/its/General_EUK_ITS_v1.9.3.fasta"
output_fasta_fungi="eukaryome/its/General_EUK_ITS_v1.9.3_fungi.fasta"
output_fasta_glomeromycota="eukaryome/its/General_EUK_ITS_v1.9.3_glomeromycota.fasta"

# Filter FASTA for entries with k__Fungi
awk '/^>/ {if ($0 ~ /k__Fungi/) print_header=1; else print_header=0}
     print_header {print}
     !/^>/ {if (print_header) print}' $input_fasta > $output_fasta_fungi

# Filter FASTA for entries with p__Glomeromycota
awk '/^>/ {if ($0 ~ /p__Glomeromycota/) print_header=1; else print_header=0}
      print_header {print}
      !/^>/ {if (print_header) print}' $input_fasta > $output_fasta_glomeromycota

# Now we'll count how many entries are in each version of the database

# Output file for the table
output_table="eukaryome/eukaryome_ITS_sequencecounts.txt"

# Header for the table
echo -e "File\tSequence_Count" > $output_table

# Loop through each specified FASTA file
for file in ${input_fasta} ${output_fasta_fungi} ${output_fasta_glomeromycota}; do
    # Count sequences in the file
    count=$(grep -c "^>" "$file")

    # Add the filename and count to the table
    echo -e "$(basename $file)\t$count" >> $output_table
done

# Display the table
cat $output_table

#------------------------ QIIME2 DATABASES -------------------------------------

#-------- LSU

# Download the file
wget https://sisu.ut.ee/wp-content/uploads/sites/643/QIIME2_EUK_LSU_v1.9.3.zip

# Unzip the file
unzip QIIME2_EUK_LSU_v1.9.3.zip

# Move the downloaded file to the correct place
mv QIIME2_EUK_LSU_v1.9.3* eukaryome/lsu/

#--------- LSU Cleanup
# Input and output file names
input_file="eukaryome/lsu/QIIME2_EUK_LSU_v1.9.3.tsv"
output_file="eukaryome/lsu/QIIME2_EUK_LSU_v1.9.3_cleaned.tsv"

# Process the file
awk 'BEGIN {OFS="\t"} {
    if (NR == 1) {print; next} # Print header as is

    # Extract Feature ID and Taxon columns
    feature_id = $1
    taxon = $2

    # Remove the final semicolon and everything after it
    sub(/;[^;]*$/, "", taxon)

    # Remove only the "." in empty taxonomic levels but keep semicolons
    gsub(/\./, "", taxon)

    # Print cleaned line
    print feature_id, taxon
}' "$input_file" > "$output_file"

echo "Processing complete. Cleaned data saved in $output_file"

#-------- SSU

# Download the file
wget https://sisu.ut.ee/wp-content/uploads/sites/643/QIIME2_EUK_SSU_v1.9.3.zip

# Unzip the file
unzip QIIME2_EUK_SSU_v1.9.3.zip

# Move the downloaded file to the correct place
mv QIIME2_EUK_SSU_v1.9.3* eukaryome/ssu/

#--------- SSU Cleanup
# Input and output file names
input_file="eukaryome/ssu/QIIME2_EUK_SSU_v1.9.3.tsv"
output_file="eukaryome/ssu/QIIME2_EUK_SSU_v1.9.3_cleaned.tsv"

# Process the file
awk 'BEGIN {OFS="\t"} {
    if (NR == 1) {print; next} # Print header as is

    # Extract Feature ID and Taxon columns
    feature_id = $1
    taxon = $2

    # Remove the final semicolon and everything after it
    sub(/;[^;]*$/, "", taxon)

    # Remove only the "." in empty taxonomic levels but keep semicolons
    gsub(/\./, "", taxon)

    # Print cleaned line
    print feature_id, taxon
}' "$input_file" > "$output_file"

echo "Processing complete. Cleaned data saved in $output_file"

#-------- ITS

# Download the file
wget https://sisu.ut.ee/wp-content/uploads/sites/643/QIIME2_EUK_ITS_v1.9.3.zip

# Unzip the file
unzip QIIME2_EUK_ITS_v1.9.3.zip

# Move the downloaded file to the correct place
mv QIIME2_EUK_ITS_v1.9.3* eukaryome/its/

#--------- ITS Cleanup
# Input and output file names
input_file="eukaryome/its/QIIME2_EUK_ITS_v1.9.3.tsv"
output_file="eukaryome/its/QIIME2_EUK_ITS_v1.9.3_cleaned.tsv"

# Process the file
awk 'BEGIN {OFS="\t"} {
    if (NR == 1) {print; next} # Print header as is

    # Extract Feature ID and Taxon columns
    feature_id = $1
    taxon = $2

    # Remove the final semicolon and everything after it
    sub(/;[^;]*$/, "", taxon)

    # Remove only the "." in empty taxonomic levels but keep semicolons
    gsub(/\./, "", taxon)

    # Print cleaned line
    print feature_id, taxon
}' "$input_file" > "$output_file"

echo "Processing complete. Cleaned data saved in $output_file"

# Now for the sake of working with a smaller database we want to subset only the fungal entries

# Define input files
FASTA_FILE="eukaryome/its/QIIME2_EUK_ITS_v1.9.3.fasta"
TAXONOMY_FILE="eukaryome/its/QIIME2_EUK_ITS_v1.9.3_cleaned.tsv"
OUTPUT_TAX="eukaryome/its/QIIME2_EUK_ITS_v1.9.3_cleaned_fungi.tsv"
OUTPUT_FASTA="eukaryome/its/QIIME2_EUK_ITS_v1.9.3_fungi.fasta"


# Step 1: Filter the taxonomy file to keep only fungal entries, preserving the header
awk 'BEGIN {OFS="\t"}
     NR==1 {print; next} # Print header
     $2 ~ /^Fungi/ {print}' "$TAXONOMY_FILE" > "$OUTPUT_TAX"

# Step 2: Extract unique identifiers (Feature IDs) from the filtered taxonomy file
cut -f1 "$OUTPUT_TAX" > ids.txt

# Step 3: Filter the fasta file based on extracted IDs
awk 'BEGIN {while ((getline < "ids.txt") > 0) ids[$1]=1}
     /^>/ {header=$0; gsub(/^>/, "", header); keep=ids[header]}
     keep {print}' "$FASTA_FILE" > "$OUTPUT_FASTA"

# Clean up temporary file
rm ids.txt

echo "Filtering complete. Results saved in $OUTPUT_TAX and $OUTPUT_FASTA."

#-------- Longread

# Download the file
wget https://sisu.ut.ee/wp-content/uploads/sites/643/QIIME2_EUK_longread_v1.9.3.zip

# Unzip the file
unzip QIIME2_EUK_longread_v1.9.3.zip

# Move the downloaded file to the correct place
mv QIIME2_EUK_longread_v1.9.3* eukaryome/longread/

#--------- Longread Cleanup
# Convert lowercase sequences to uppercase
# Input and output file names
input_file="eukaryome/longread/QIIME2_EUK_longread_v1.9.3.fasta"
output_file="eukaryome/longread/QIIME2_EUK_longread_v1.9.32.fasta"

# Read the input file line by line
while IFS= read -r line; do
    # If the line starts with '>', it's a header line, so just print it as is
    if [[ "$line" =~ ^\> ]]; then
        echo "$line" >> "$output_file"
    else
        # For sequence lines, convert lowercase to uppercase and print
        echo "$line" | tr 'a-z' 'A-Z' >> "$output_file"
    fi
done < "$input_file"

#mv eukaryome/longread/QIIME2_EUK_longread_v1.9.32.fasta eukaryome/longread/QIIME2_EUK_longread_v1.9.3.fasta

# Use awk to exclude the sequence with header >CP048731
# This entry has "-" in it which messes with the qiime2 taxonomy... I don't know how this made its wau into the file
# Process the input file and remove sequences containing '-'
awk '/^>/ {if (seq && !found_dash) print header "\n" seq; header = $0; seq = ""; found_dash = 0; next}
     !/^>/ {seq = seq $0; if ($0 ~ /-/) found_dash = 1}
     END {if (seq && !found_dash) print header "\n" seq}' "$output_file" > "$input_file"


# Input and output file names
input_file="eukaryome/longread/QIIME2_EUK_longread_v1.9.3.tsv"
output_file="eukaryome/longread/QIIME2_EUK_longread_v1.9.3_cleaned.tsv"

# Process the file
awk 'BEGIN {OFS="\t"} {
    if (NR == 1) {print; next} # Print header as is

    # Extract Feature ID and Taxon columns
    feature_id = $1
    taxon = $2

    # Remove the final semicolon and everything after it
    sub(/;[^;]*$/, "", taxon)

    # Remove only the "." in empty taxonomic levels but keep semicolons
    gsub(/\./, "", taxon)

    # Print cleaned line
    print feature_id, taxon
}' "$input_file" > "$output_file"

echo "Processing complete. Cleaned data saved in $output_file"

# Now for the sake of working with a smaller database we want to subset only the fungal entries

# Define input files
FASTA_FILE="eukaryome/longread/QIIME2_EUK_longread_v1.9.3.fasta"
TAXONOMY_FILE="eukaryome/longread/QIIME2_EUK_longread_v1.9.3_cleaned.tsv"
OUTPUT_TAX="eukaryome/longread/QIIME2_EUK_longread_v1.9.3_cleaned_fungi.tsv"
OUTPUT_FASTA="eukaryome/longread/QIIME2_EUK_longread_v1.9.3_fungi.fasta"


# Step 1: Filter the taxonomy file to keep only fungal entries, preserving the header
awk 'BEGIN {OFS="\t"}
     NR==1 {print; next} # Print header
     $2 ~ /^Fungi/ {print}' "$TAXONOMY_FILE" > "$OUTPUT_TAX"

# Step 2: Extract unique identifiers (Feature IDs) from the filtered taxonomy file
cut -f1 "$OUTPUT_TAX" > ids.txt

# Step 3: Filter the fasta file based on extracted IDs
awk 'BEGIN {while ((getline < "ids.txt") > 0) ids[$1]=1}
     /^>/ {header=$0; gsub(/^>/, "", header); keep=ids[header]}
     keep {print}' "$FASTA_FILE" > "$OUTPUT_FASTA"

# Clean up temporary file
rm ids.txt

echo "Filtering complete. Results saved in $OUTPUT_TAX and $OUTPUT_FASTA."
