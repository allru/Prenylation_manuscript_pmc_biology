import pandas as pd
from Bio import SeqIO # NOTE: Biopython needs to be installed!
import re


def is_protein_id_valid(protein_id, c_positions):
    
    is_valid = False
    
    if protein_id in c_positions.index:
        is_valid = True
    else:
        print(f"ERROR: {protein_id} is not in c position list. Removing protein from analysis.")
    
    return(is_valid)



def is_sequence_valid(sequence_id, sequence_length, c_positions):
    
    is_valid = False
    reference_sequence_length = c_positions.loc[sequence_id]["len"].unique()
    
    if sequence_length == reference_sequence_length:
        is_valid = True
    else:
        print("ERROR: Length of Alphafold sequence does not match length of reference sequence.")
        print(f"Masking entire {sequence_id} sequence; protein is ignored in subsequent analysis")
        print("Alphafold sequence length: "+ str(sequence_length))
        print("Reference sequence length: "+ str(reference_sequence_length))
        
    return(is_valid)
      
      

# Reads in a fasta file with 3di or aa sequences from Foldseek and masks all 3di states/aa residues outside range
# At the moment, range is defined by a certain distance (in aa) from a C of interest
def mask_all_sequences_in_fastafile(input_fasta_filepath, 
                                    c_positions, 
                                    local_range,
                                    output_state_filepath = None):
    
    # OUTPUT: csv file with 3di sequences around C of interest       
    if not output_state_filepath is None:
        with open(output_state_filepath, "w") as out_state_file:
            out_state_file.write("id,c_position,state\n")


    # Parse 3Di FASTA file sequence by sequence
    for seq_record in SeqIO.parse(input_fasta_filepath, "fasta"):
        
        pdb_filename    = seq_record.id           # Get sequence ID from FASTA (= name of PDB file)
        sequence_header = seq_record.description  # Get full header of sequence
        sequence_length = len(seq_record)         # Get length of 3Di sequence
        sequence        = seq_record.seq          # Get 3Di sequence
        
        # Make sure pdb file name matches c_positions ID
        # TODO: make sure this step is correct!! Find more robust way to match names to IDs?
        
        sequence_id = re.findall(r"-([^-]*)-",pdb_filename)[0]
        
        # Check if protein ID from FASTA file matches protein ID name in c_positions file  
        if not is_protein_id_valid(sequence_id, c_positions):
            continue
        
        # Check if length of 3Di sequence matches length of amino acid sequence
        # If length is mismatched, the entire sequence is masked and a warning message is displayed.
        if not is_sequence_valid(sequence_id, sequence_length, c_positions):
            masked_3di_sequence = str(sequence).lower()
            
        else:

            
        # For each protein, process all C sites defined in the input csv
            for c_position in c_positions.loc[sequence_id].index:
        
                array_c_position = c_position - 1  # Minus one because Python starts counting at 0

                start = max(0, array_c_position-local_range)
                end   = min(len(sequence), array_c_position+local_range+1)

                masked_sequence = sequence[start:end]
                
                # Convert 3Di alphabet to greek letters (to avoid confusion with amino acids)
                greek_alphabet = 'αβγδεηθικλμνπρστφχψω'
                latin_alphabet = 'avgdehciklmnprstfqyw'
                greek2latin = str.maketrans(latin_alphabet, greek_alphabet)
                
                greek_masked_seq = str(masked_sequence).lower().translate(greek2latin)

                if not output_state_filepath is None:
                    with open(output_state_filepath, "a", encoding="utf-8") as out_state_file:
                        out_state_file.write(f"{sequence_id},{c_position},{greek_masked_seq}\n")
            
pass
       