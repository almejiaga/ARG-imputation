#!/usr/bin/env python
# coding: utf-8

import argparse
import pandas as pd
import tskit

def extract_numbers(sample):
    """Extracts numerical sample IDs from a sample list."""
    if pd.isna(sample):
        return []
    return [int(x.split('_')[1]) for x in sample.split(',') if x.startswith('Sample_')]

#importing all the libraries needed for the analysis
import numpy as np
import tskit
import time
import math
import itertools

def find_min_subset_tmrca_with_mrcas(ids_list, tree, tmrca_threshold=None):
    min_tmrca = float('inf')
    min_subset = None
    
    # Generate all possible subsets of non-consecutive IDs with size equal to half of the original list size
    subset_size = len(ids_list) // 2
    all_subsets = itertools.combinations(ids_list, subset_size)
    
    # Iterate over all subsets
    for subset in all_subsets:
        # Check if the IDs are consecutive, skip if they are
        if any(abs(subset[i] - subset[i+1]) == 1 for i in range(len(subset)-1)):
            continue
        
        # Calculate the total TMRCA for the subset
        total_tmrca = sum(tree.tmrca(int(id1), int(id2)) for id1, id2 in itertools.combinations(subset, 2))
        
        # Update min_tmrca and min_subset if necessary
        if total_tmrca < min_tmrca:
            min_tmrca = total_tmrca
            min_subset = subset
    
    if min_subset:
        # Initialize variables for highest pairwise TMRCA and corresponding pair
        highest_pair_tmrca = -float('inf')
        best_pair = None
        
        # Calculate pairwise TMRCA within min_subset
        for pair in itertools.combinations(min_subset, 2):
            tmrca = tree.tmrca(pair[0], pair[1])
            if tmrca > highest_pair_tmrca and (tmrca_threshold is None or tmrca < tmrca_threshold):
                highest_pair_tmrca = tmrca
                best_pair = pair
        
        # Find MRCA of the pair with the highest TMRCA
        if best_pair:
            mrca = tree.mrca(best_pair[0], best_pair[1])
            return min_subset, min_tmrca, best_pair, highest_pair_tmrca, mrca
        else:
            return min_subset, min_tmrca, None, None, None
    else:
        return None, None, None, None, None


# Function to extract numbers from 'HetSamples' column in the dataframe (sometimes we wont use it, but it works)
def extract_numbers(sample):
    if pd.isna(sample):
        return []
    else:
        numbers = [int(x.split('_')[1]) for x in sample.split(',') if x.startswith('Sample_')]
        return numbers

# Function to find node with time closest before a target time t, which in this case is 100
def find_node_with_time_closest_before(tree, node_id, target_time):
    current_node_id = node_id
    prev_node_id = None
    
    while True:
        current_time = tree.time(int(current_node_id))  # Cast current_node_id to int
        
        if current_time >= target_time:
            return prev_node_id, tree.num_samples(int(prev_node_id))
        
        prev_node_id = current_node_id
        parent_node_id = tree.parent(int(current_node_id))  # Cast current_node_id to int
        
        if parent_node_id is None:
            return prev_node_id
        
        current_node_id = parent_node_id
# Function to find node with 75% clade size and also the time of that node (the founder score)
def find_node_with_75_percent_clade_size(tree, node_id, clade_size_threshold):
    if clade_size_threshold < 4:
        return math.nan, math.nan  # Output nan, nan if clade size threshold is less than 4
    
    current_node_id = node_id
    parent_node_id = tree.parent(int(current_node_id))  # Cast current_node_id to int
    
    while parent_node_id is not None:
        parent_clade_size = tree.num_samples(int(parent_node_id))
        if parent_clade_size >= clade_size_threshold * 0.749:
            return parent_node_id, tree.time(parent_node_id)
        current_node_id = parent_node_id
        parent_node_id = tree.parent(int(current_node_id))  # Cast current_node_id to int


    # If we reach here, it means we reached the root without finding a node with the desired clade size
    return math.nan, math.nan  # Output nan, nan if no node found

def clean_node_ids(node_ids):
    """Remove NaNs and convert all elements to integers."""
    return [int(node_id) for node_id in node_ids if not pd.isna(node_id)]

def find_min_tmrca(tree, node_of_interest, all_ids):
    """Find the node with the smallest TMRCA with the node of interest."""
    min_tmrca = float('inf')
    min_tmrca_node = None
    for node_id in all_ids:
        tmrca_value = tree.tmrca(node_of_interest, node_id)
        if tmrca_value < min_tmrca:
            min_tmrca = tmrca_value
            min_tmrca_node = node_id
    return min_tmrca_node, min_tmrca

def calculate_pairwise_tmrca(tree, node_of_interest, nodes):
    """Calculate pairwise TMRCA with node of interest."""
    return {node: tree.tmrca(node_of_interest, node) for node in nodes}

def assign_probabilities(min_tmrca, tmrca_results, tree, node_of_interest):
    """Assign probability of being a carrier."""
    node_time = tree.time(node_of_interest)
    probabilities = {}
    for node, tmrca_result in tmrca_results.items():
        probability = (min_tmrca - tmrca_result) / (min_tmrca - node_time)
        probabilities[node] = probability
    return probabilities

def main_pipeline(finaldf, tree, node_of_interest):
    # Step 1: Clean Data
    row_index = finaldf.index[0]
    all_ids = finaldf.loc[row_index, 'all_hom_ids']
    all_ids = clean_node_ids(all_ids)
    imputed_carriers = list(tree.samples(node_of_interest))
    set_all_ids = set(all_ids)
    set_imputed_carriers = set(imputed_carriers)

    # Find intersection
    intersection = set_all_ids.intersection(set_imputed_carriers)
    intersection_list = list(intersection)
    all_ids = list(set(all_ids) - set(intersection_list)) 
    # Step 2: Find Node with Smallest TMRCA
    min_tmrca_node, min_tmrca = find_min_tmrca(tree, node_of_interest, all_ids)
    print(f"The node with the smallest TMRCA is {min_tmrca_node} with a TMRCA of {min_tmrca}")
    
    # Step 3: Find Node with Time Closest Before
    result = find_node_with_time_closest_before(tree, node_of_interest, min_tmrca)
    
    # Step 4: Generate Sample Lists
    all_possible = list(tree.samples(result[0]))
    carriers = list(tree.samples(node_of_interest))
    
    # Step 5: Calculate Pairwise TMRCA
    non_carriers = set(all_possible) - set(carriers)
    tmrca_results = calculate_pairwise_tmrca(tree, node_of_interest, non_carriers)
    
    # Step 6: Assign Probability
    probabilities = assign_probabilities(min_tmrca, tmrca_results, tree, node_of_interest)
    data = []
    
    # Add non-carriers with their calculated probabilities
    for node_id, probability in probabilities.items():
        data.append({'node id': node_id, 'posterior probability': probability})
    
    # Add carriers with a probability of 1
    for carrier in carriers:
        data.append({'node id': carrier, 'posterior probability': 1})
    
    # Create DataFrame
    df = pd.DataFrame(data)
    return df
    

def main(chromosome, variant):
    base_path = "/lustre06/project/6068353/almejiaga/projects/Founder_score/results/"
    cartagene_path = "/lustre06/project/6068353/almejiaga/cartagenegenotypedata/VCFs_2/ARG_files/"

    # Load chromosome-specific data
    hets_file = f"{base_path}data/processed/hetsids/homhetimputeddata_chr{chromosome}.txt"
    df = pd.read_csv(hets_file, comment="#", sep="\t", header=0)

    # Select relevant columns
    finaldf2 = df[["ID", "POS", "nHet", "HetSamples", "nHomAlt", "HomSamplesAlt", "HomSamplesRef"]]

    # Filter by the specified variant
    finaldf = finaldf2[finaldf2['ID'] == variant]

    # Process het and hom carriers
    finaldf['all_het_carriers'] = finaldf['HetSamples'].apply(extract_numbers)
    finaldf['all_hom_samples'] = finaldf['HomSamplesRef'].apply(extract_numbers)
    finaldf['all_carriers'] = finaldf['HetSamples'].apply(extract_numbers)
    # Read the file into a DataFrame
    dictionary = pd.read_csv('/lustre06/project/6068353/almejiaga/projects/Founder_score/results/data/raw/dictionary_args.txt', delim_whitespace=True)

    #checking the clusters
    #reading the csv metadata from alex
    metadata = pd.read_csv('~/projects/ctb-sgravel/cartagene/research/quebec_structure_936028/Alex_CaG-UMAP-clusters/cartagene_clusters_eth_cob.csv', comment="#", sep=",", header=0)
    # Load ARG tree
    ts_file = f"{cartagene_path}chr{chromosome}_final.trees"
    ts = tskit.load(ts_file)

    # Load position correction map
    map_file = f"{cartagene_path}chr{chromosome}_corrected.map"
    MAP = pd.read_csv(map_file, delim_whitespace=True, header=None)
    correction_value = MAP.iloc[0, 3]

    # Correct positions
    finaldf["Corrected_POS"] = finaldf["POS"] - correction_value
    finaldf["Corrected_POS"] = finaldf["Corrected_POS"].clip(lower=1)
    print(finaldf)
    # Create empty lists to store the firstID and secondID values
    first_ids = []
    second_ids = []

    # Iterate over the rows of the nonsingletons DataFrame
    for index, row in finaldf.iterrows():
        # Extract the list of IDs from the all_het_carriers column
        id_list = row['all_carriers']

        # Initialize variables to store the firstID and secondID values
        first_id = []
        second_id = []

        # Look up each ID in the dictionary DataFrame and store the corresponding firstID and secondID values
        for id_ in id_list:
            # Check if the ID exists in the dictionary DataFrame
            if id_ in dictionary['ID_2'].values:
                # Get the corresponding firstID and secondID values
                first_id.append(dictionary.loc[dictionary['ID_2'] == id_, 'firstID'].iloc[0])
                second_id.append(dictionary.loc[dictionary['ID_2'] == id_, 'secondID'].iloc[0])
            else:
                # If the ID does not exist in the dictionary, append None values
                first_id.append(None)
                second_id.append(None)

        # Append the lists of firstID and secondID values to the respective lists
        first_ids.append(first_id)
        second_ids.append(second_id)

    # Create a new list to store the combined IDs
    all_ids = []

    # Iterate over the rows of the nonsingletons DataFrame
    for first, second in zip(first_ids, second_ids):
        # Combine the first and second lists and remove any None values
        combined_ids = [x for x in first + second if x is not None]

        # Append the combined list to the all_ids list
        all_ids.append(combined_ids)

    # Add the all_ids list as a new column in the nonsingletons DataFrame
    finaldf['all_carriers_ids'] = all_ids
    # Create empty lists to store the firstID and secondID values
    first_ids = []
    second_ids = []

    # Iterate over the rows of the nonsingletons DataFrame
    for index, row in finaldf.iterrows():
        # Extract the list of IDs from the all_het_carriers column
        id_list = row['all_hom_samples']

        # Initialize variables to store the firstID and secondID values
        first_id = []
        second_id = []

        # Look up each ID in the dictionary DataFrame and store the corresponding firstID and secondID values
        for id_ in id_list:
            # Check if the ID exists in the dictionary DataFrame
            if id_ in dictionary['ID_2'].values:
                # Get the corresponding firstID and secondID values
                first_id.append(dictionary.loc[dictionary['ID_2'] == id_, 'firstID'].iloc[0])
                second_id.append(dictionary.loc[dictionary['ID_2'] == id_, 'secondID'].iloc[0])
            else:
                # If the ID does not exist in the dictionary, append None values
                first_id.append(None)
                second_id.append(None)

        # Append the lists of firstID and secondID values to the respective lists
        first_ids.append(first_id)
        second_ids.append(second_id)

    # Create a new list to store the combined IDs
    all_ids = []

    # Iterate over the rows of the nonsingletons DataFrame
    for first, second in zip(first_ids, second_ids):
        # Combine the first and second lists and remove any None values
        combined_ids = [x for x in first + second if x is not None]

        # Append the combined list to the all_ids list
        all_ids.append(combined_ids)

    # Add the all_ids list as a new column in the nonsingletons DataFrame
    finaldf['all_hom_ids'] = all_ids
    row_index = finaldf.index[0]
    print(row_index)
    final_list_mut =finaldf.loc[row_index,"all_carriers_ids"]
    final_list_mut
    final_list_mut.sort()
    final_list_mut
    #Converting to a list of integers
    int_list = [int(value) for value in final_list_mut]
    print(int_list)
    len(int_list)
    tree3= ts.at(finaldf.loc[row_index, "Corrected_POS"])
        # Example usage with a threshold
    ids_list = int_list  # Replace with your list of IDs
    tmrca_threshold = 100  # Set your desired threshold here
    min_subset, min_tmrca, best_pair, highest_pair_tmrca, mrca = find_min_subset_tmrca_with_mrcas(ids_list, tree3, tmrca_threshold)

    # Print results
    print("Smallest TMRCA subset:", min_subset)
    print("TMRCA of the smallest subset:", min_tmrca)
    print("Pair with highest TMRCA below", tmrca_threshold, "in the subset:", best_pair)
    print("Highest TMRCA within the subset below", tmrca_threshold, ":", highest_pair_tmrca)
    print("MRCA of the pair with highest TMRCA:", mrca)
    

    # Assuming your DataFrame is named 'finaldf' and the 'all_ids' column contains the lists
    all_ids = finaldf.loc[row_index, 'all_hom_ids']
    # Remove NaNs and convert all elements in the list to integers (remove .0)
    cleaned_all_ids = all_ids

    # Assuming tree3 is your tree object and it has a method `tmrca`
    node_of_interest = mrca

    # Remove NaNs and convert all elements to integers
    all_ids = [int(node_id) for node_id in all_ids if not pd.isna(node_id)]
    # Assuming all_ids is your list
    # Initialize variables to store the minimum TMRCA and corresponding node
    min_tmrca = float('inf')
    min_tmrca_node = None

    # Iterate through each node_id in the list of IDs
    for node_id in all_ids:
        tmrca_value = tree3.tmrca(node_of_interest, node_id)
        if tmrca_value < min_tmrca:
            min_tmrca = tmrca_value
            min_tmrca_node = node_id

    print(f"The node with the smallest TMRCA is {min_tmrca_node} with a TMRCA of {min_tmrca}")
    # Convert both lists to sets to find the intersection
    imputed_carriers = list(tree3.samples(mrca))
    set_all_ids = set(all_ids)
    set_imputed_carriers = set(imputed_carriers)

    # Find intersection
    intersection = set_all_ids.intersection(set_imputed_carriers)

    # Alternatively, using the & operator
    # intersection = set_all_ids & set_imputed_carriers

    # Convert back to a list if needed
    intersection_list = list(intersection)

    # Check the result
    print("Intersection:", intersection_list)

    node_of_interest = mrca

    # Run the pipeline
    results = main_pipeline(finaldf, tree3, node_of_interest)

# Create the DataFrame
#probability_df = create_probability_dataframe(non_carriers, probabilities, carriers)

# Print the DataFrame
    results
    exclude_ids = intersection_list

    # Exclude rows where 'node id' is in exclude_ids
    results = results[~results['node id'].isin(exclude_ids)]

    # Display the filtered DataFrame
    #print(results)
    # Merge based on 'firstID'
    merged_df1 = pd.merge(dictionary, results, left_on='firstID', right_on='node id', how='left')
    merged_df1.rename(columns={'node id': 'ID', 'posterior probability': 'posterior'}, inplace=True)

    # Merge based on 'secondID'
    merged_df2 = pd.merge(dictionary, results, left_on='secondID', right_on='node id', how='left')
    merged_df2.rename(columns={'node id': 'ID', 'posterior probability': 'posterior'}, inplace=True)

    # Combine the results, ensuring no duplicates
    merged_dfinal = pd.concat([merged_df1, merged_df2]).drop_duplicates().reset_index(drop=True)

    # Fill NaN values in 'ID' and 'posterior' with appropriate defaults
    merged_dfinal['ID'] = merged_dfinal['ID'].fillna(0)  # Fill with 0 or any default value as needed
    merged_dfinal['posterior'] = merged_dfinal['posterior'].fillna(0)  # Fill with 0 or any default value as needed

    # Display the final merged dataframe
    print(merged_dfinal)


    # Read the gzipped CSV file into a DataFrame
    df2 = pd.read_csv('~/projects/ctb-sgravel/cartagene/research/quebec_structure_936028/data/phenotypes/Gravel936028_5/Gravel_936028_5.mesures_biochimiques.csv.gz', compression='gzip')

    # Now you can work with the DataFrame 'df'
    #carriers = pd.read_csv('/lustre06/project/6068353/almejiaga/projects/Founder_score/results/scripts/NR1H3carriers.txt', header=0)
    # Example data
    # df = pd.DataFrame({'file111': ['id1', 'id2', 'id3', 'id4', 'id5']})
    # carriers = pd.DataFrame({'carrier_ids': ['id1', 'id3']})

    # Replace the above example with your actual data
    # df = pd.read_csv('your_dataframe_file.csv')
    # carriers = pd.read_csv('your_carriers_file.csv')

    # List of carrier IDs
    carriers = merged_dfinal["ID_2"].tolist()
    carrier_ids = [int(s) for s in carriers]

    #checking the clusters
    #reading the csv metadata from alex
    
    metadata = pd.read_csv('~/projects/ctb-sgravel/cartagene/research/quebec_structure_936028/Alex_CaG-UMAP-clusters/cartagene_clusters_eth_cob.csv', comment="#", sep=",", header=0)
    metadata
    filtered_metadata3 = metadata[metadata['IID'].isin(carrier_ids)]
    carriers_list= filtered_metadata3["PROJECT_CODE"].to_list()
    len(carriers_list)
    merged_df3 = pd.merge(merged_dfinal, metadata, left_on='ID_2', right_on='IID', how='left')
    merged_df3
    frequency_table = merged_df3['cluster'].value_counts()
    print(frequency_table)
    # Assuming merged_df3 is your dataframe and it contains columns 'Cluster' and 'posterior'
    # Group by 'Cluster' and compute mean posterior value
    cluster_posterior_stats = merged_df3['posterior'].sum()

    # Optionally, you can compute other statistics such as median, standard deviation, etc.
    # cluster_posterior_median = merged_df3.groupby('Cluster')['posterior'].median()
    # cluster_posterior_std = merged_df3.groupby('Cluster')['posterior'].std()

    # Display the results
    print("Mean Posterior Value per Cluster:")
    print(cluster_posterior_stats)
    # Assuming merged_df3 is your DataFrame
    subset_non_zero_posterior = merged_df3[merged_df3['posterior'] != 0]

    # Display the subset or use it as needed
    subset_non_zero_posterior.groupby('cluster')['posterior'].sum()
    subset_non_zero_posterior_final = subset_non_zero_posterior[['ID_2', 'posterior']]

    # Display the subset or use it as needed
    subset_non_zero_posterior_final
    subset_non_zero_posterior_final.to_csv(f'posterior_{variant}.csv', index=False)
    frequency_table = merged_df3['cluster'].value_counts()
    print(frequency_table)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Pipeline for analyzing SACs with ARG integration.")
    parser.add_argument("--chr", required=True, type=int, help="Chromosome number (1-22).")
    parser.add_argument("--variant", required=True, help="Variant in format chrX_POS_REF_ALT (e.g., chr8_19954279_C_T).")

    args = parser.parse_args()
    main(args.chr, args.variant)
