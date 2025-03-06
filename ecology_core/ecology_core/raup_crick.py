"""
Raup-Crick dissimilarity implementation.

This module provides a Python implementation of the Raup-Crick dissimilarity metric,
which is a probabilistic measure of community similarity that accounts for differences
in richness and the selection of species from a common species pool.
"""

import numpy as np
import pandas as pd
from scipy.stats import hypergeom
from tqdm.auto import tqdm
from .data_structures import EcologicalData


def raup_crick_abundance(data, n_sim=999, null_model="proportional", rarefy_depth=None, seed=None):
    """
    Calculate Raup-Crick dissimilarity accounting for species abundances.
    
    Parameters:
        data (EcologicalData): Ecological data object
        n_sim (int): Number of simulations for null model
        null_model (str): Null model type:
            - "proportional": Species occurrence probability proportional to frequency
            - "equiprobable": All species equally likely to occur
        rarefy_depth (int): Depth for rarefaction before calculating dissimilarity
        seed (int): Random seed for reproducibility
    
    Returns:
        pd.DataFrame: Raup-Crick dissimilarity matrix
    """
    # Ensure data is an EcologicalData object
    if not isinstance(data, EcologicalData):
        raise TypeError("data must be an EcologicalData object")
    
    # Set random seed
    np.random.seed(seed)
    
    # Rarefy if specified
    if rarefy_depth is not None:
        data = data.rarefy(depth=rarefy_depth)
    
    # Extract abundance data and convert to presence/absence
    abundance = data.abundance
    presence_absence = (abundance > 0).astype(int)
    
    # Number of samples and features
    n_samples = presence_absence.shape[1]
    n_features = presence_absence.shape[0]
    
    # Calculate feature frequencies (proportion of samples where each feature is present)
    feature_freq = presence_absence.mean(axis=1)
    
    # Calculate sample richness (number of features in each sample)
    sample_richness = presence_absence.sum(axis=0)
    
    # Calculate observed shared features between sample pairs
    observed_shared = np.zeros((n_samples, n_samples))
    for i in range(n_samples):
        for j in range(i+1, n_samples):
            observed_shared[i, j] = np.sum(presence_absence.iloc[:, i] & presence_absence.iloc[:, j])
            observed_shared[j, i] = observed_shared[i, j]  # Symmetric matrix
    
    # Calculate expected shared features and variance under null model
    # using randomization simulations
    simulated_shared = np.zeros((n_samples, n_samples, n_sim))
    
    # Progress bar for simulations
    for sim in tqdm(range(n_sim), desc="Raup-Crick simulations"):
        # Create simulated communities
        sim_communities = np.zeros((n_features, n_samples), dtype=int)
        
        # For each sample, draw features based on its richness
        for j in range(n_samples):
            richness = sample_richness.iloc[j]
            
            # Set probability based on null model
            if null_model == "proportional":
                probs = feature_freq.values
            else:  # equiprobable
                probs = np.ones(n_features) / n_features
            
            # Draw features without replacement
            features = np.random.choice(
                np.arange(n_features), 
                size=richness, 
                replace=False, 
                p=probs / sum(probs)
            )
            sim_communities[features, j] = 1
        
        # Calculate shared features for this simulation
        for i in range(n_samples):
            for j in range(i+1, n_samples):
                shared = np.sum(sim_communities[:, i] & sim_communities[:, j])
                simulated_shared[i, j, sim] = shared
                simulated_shared[j, i, sim] = shared  # Symmetric
    
    # Calculate Raup-Crick index
    rc_matrix = np.zeros((n_samples, n_samples))
    
    for i in range(n_samples):
        for j in range(i+1, n_samples):
            # Count number of simulations with fewer shared species than observed
            num_less = np.sum(simulated_shared[i, j, :] < observed_shared[i, j])
            
            # Count number of simulations with same number of shared species as observed
            num_equal = np.sum(simulated_shared[i, j, :] == observed_shared[i, j])
            
            # Calculate RC index following Chase et al. 2011
            # RC = (number less than observed + half of number equal to observed) / total number of simulations
            # Scale from 0-1 to -1 to 1
            rc = ((num_less + (0.5 * num_equal)) / n_sim) * 2 - 1
            
            rc_matrix[i, j] = rc
            rc_matrix[j, i] = rc  # Symmetric
    
    # Convert to DataFrame
    rc_df = pd.DataFrame(
        rc_matrix,
        index=abundance.columns,
        columns=abundance.columns
    )
    
    return rc_df


def raup_crick_beta_diversity(rc_matrix):
    """
    Interpret Raup-Crick values in terms of beta diversity components.
    
    Parameters:
        rc_matrix (pd.DataFrame): Raup-Crick dissimilarity matrix
    
    Returns:
        pd.DataFrame: Beta diversity interpretation with stochastic, deterministic values
    """
    # Extract unique sample pairs (upper triangle of the matrix)
    n_samples = rc_matrix.shape[0]
    pairs = []
    
    for i in range(n_samples):
        for j in range(i+1, n_samples):
            pairs.append({
                'sample1': rc_matrix.index[i],
                'sample2': rc_matrix.index[j],
                'rc_value': rc_matrix.iloc[i, j]
            })
    
    results = pd.DataFrame(pairs)
    
    # Add interpretation columns
    # RC values close to 0 (-0.1 to 0.1) indicate stochastic assembly
    # RC values < -0.1 indicate communities more similar than expected (underdispersion)
    # RC values > 0.1 indicate communities less similar than expected (overdispersion)
    results['assembly'] = 'stochastic'
    results.loc[results['rc_value'] < -0.1, 'assembly'] = 'underdispersed'
    results.loc[results['rc_value'] > 0.1, 'assembly'] = 'overdispersed'
    
    # Calculate summary statistics
    summary = {
        'underdispersed': sum(results['assembly'] == 'underdispersed') / len(results),
        'stochastic': sum(results['assembly'] == 'stochastic') / len(results),
        'overdispersed': sum(results['assembly'] == 'overdispersed') / len(results),
        'mean_rc': results['rc_value'].mean(),
        'median_rc': results['rc_value'].median()
    }
    
    return results, summary