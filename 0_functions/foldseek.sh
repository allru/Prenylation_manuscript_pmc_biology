#!/bin/sh

# Set the search directory to the current directory
SEARCH_DIR="data/Foldseek"

# Find all directories and execute a command for each one
find "$SEARCH_DIR" -type d | while read -r dir; do
    echo "Processing directory: $dir"
    
    # Change to the directory being processed
    cd "$dir" || continue

    # Check if c_positions.csv exists in the directory
    if [ -f c_positions.csv ]; then
        # Create necessary subdirectories
        mkdir -p pdbfiles database cluster_subdatabases fastafiles plots

        # Make a list of Uniprot IDs from the c_positions.csv file and save it as a separate text file
        awk -F "\"*,\"*" '{print $1}' c_positions.csv | tail -n +2 | sort -u > uniprot_ids.txt

        # Download the PDB files from Alphafold and save them in the pdbfiles directory
        for i in $(cat uniprot_ids.txt); do
            wget -P pdbfiles/ "https://alphafold.ebi.ac.uk/files/AF-${i}-F1-model_v4.pdb"
        done

        # Create a 3Di Sequence Database using Foldseek and Convert it to FASTA format
        foldseek createdb pdbfiles/*.pdb database/db > database/createdb.log

        # Get the 3Di sequences in FASTA format
        foldseek lndb database/db_h database/db_ss_h
        foldseek convert2fasta database/db_ss fastafiles/unmasked_3di_seq.fasta

        # Remove the database and pdb files that are not necessary anymore
        rm -rf pdbfiles database cluster_subdatabases
    else
        echo "c_positions.csv not found in $dir"
    fi

    # Change back to the original directory
    cd - > /dev/null || exit
done
