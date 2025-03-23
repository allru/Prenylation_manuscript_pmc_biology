import pandas as pd
import tempfile
import os
import structuremap.utils
import numpy as np
from functions import pep_intern
from functions import annotate_pep_internal
from structuremap.processing import download_alphafold_cif, download_alphafold_pae, format_alphafold_data, annotate_accessibility, get_smooth_score

structuremap.utils.set_logger()

# Set output directory to tempdir
output_dir = tempfile.gettempdir()
datafolder = 'data'


def annotate_structuremap(protein_df, bonded_df, group=None, max_dist=30, max_angle=180):

    """
    Annotate protein data using the StructureMap pipeline.

    Parameters:
        protein_df (pd.DataFrame): DataFrame containing protein information.
        group (str): Group label for naming the output .csv.
        bonded_df (pd.DataFrame): DataFrame containing cysteines in disulfide bonds information.
        max_dist (int): Distance value in Angström for calculating cysteine accessibility.
        max_angle (int): Angle for calculating cysteine accessibility.

    Returns:
        None
    """
    test_proteins = (' '.join([s for s in protein_df['ID'].unique()])).split()

    # Download AlphaFold data
    cif_dir = os.path.join(output_dir, 'tutorial_cif')
    pae_dir = os.path.join(output_dir, 'tutorial_pae')
    
    valid_proteins_cif, invalid_proteins_cif, existing_proteins_cif = download_alphafold_cif(proteins=test_proteins, out_folder=cif_dir)
    
    # return invalid_proteins_cif
    valid_proteins_pae, invalid_proteins_pae, existing_proteins_pae = download_alphafold_pae(proteins=test_proteins, out_folder=pae_dir)
    
    # Format AlphaFold data input
    alphafold_annotation = format_alphafold_data(directory=cif_dir, protein_ids=test_proteins)
    
    # Keep only rows with cysteine
    alphafold_annotation = alphafold_annotation[alphafold_annotation['AA'] == 'C']
    
    # Annotate pPSE values (prediction-aware sphere exposure), using full sphere exposure
    full_sphere_exposure = annotate_accessibility(df=alphafold_annotation, max_dist=max_dist, max_angle=max_angle, error_dir=pae_dir)
    alphafold_accessibility = alphafold_annotation.merge(full_sphere_exposure, how='left', 
                                                         on=['protein_id', 'AA', 'position']).reset_index(drop=True)
    
    # Filter cysteines in disulfide bonds
    bonded = pd.read_csv(os.path.join(datafolder, 'UniProt_SPARQL_queries', bonded_df),
                             sep=';')
    bonded = bonded.rename(columns={'ID': 'protein_id', 'bond_Cpos': 'position'})
    discard = bonded.merge(alphafold_accessibility, on=['protein_id', 'position'])

    alphafold_accessibility = pd.concat([alphafold_accessibility, discard],
                                        ignore_index=True).drop_duplicates(keep=False).reset_index(drop=True)
     
    # add cysteine knowledge from previous analysis for each ID
    protein_df = protein_df.rename(columns={'ID': 'protein_id'})
    alphafold_accessibility = pd.merge(alphafold_accessibility, protein_df, how='left', on='protein_id')
    
    # Keep only the rows where the cysteine position from Alphafold matches the cysteine positions determined by us
    alphafold_accessibility = alphafold_accessibility[alphafold_accessibility['position'] == 
                                                      alphafold_accessibility['N_Cpos']].reset_index(drop=True)
    
    # Change seqID to be a better readable name
    alphafold_accessibility['name'] = alphafold_accessibility['seqID'].str.replace('_', '|').str.split('|').str[2]
    # Change back to more commonly used ID
    alphafold_accessibility = alphafold_accessibility.rename(columns={"protein_id": "ID"})
    
    alphafold_accessibility = alphafold_accessibility[['ID', 'name', 'quality', 'structure_group',
                                                       'Ccount', 'Count_all', f'nAA_{max_dist}_180_pae', 
                                                       'len', 'Cpos', 'N_Cpos', 'pep', 'pepCCC', 'pepCC', 'pepCXC', 'pepC',]]
        
    # save df
    alphafold_accessibility.to_csv(os.path.join(datafolder, 'AlphaFold_annotated', f'{group}_structure.csv'), sep=',',
                                   index=False)
    
    return alphafold_accessibility
