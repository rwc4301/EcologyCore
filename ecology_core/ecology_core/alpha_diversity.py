"""
Alpha diversity metrics for ecological data.

This module provides functions to calculate various alpha diversity metrics,
which measure the diversity within a single sample.
"""

import numpy as np
import pandas as pd
import scipy.stats as stats
from skbio.diversity import alpha
from .data_structures import EcologicalData


def calculate_alpha_diversity(data, metrics=None, rarefy_depth=None):
    """
    Calculate alpha diversity metrics for each sample.
    
    Parameters:
        data (EcologicalData): Ecological data object
        metrics (list): List of metrics to calculate, options:
            - "observed": Observed OTUs (richness)
            - "shannon": Shannon diversity index
            - "simpson": Simpson diversity index
            - "invsimpson": Inverse Simpson diversity index
            - "chao1": Chao1 richness estimator
            - "ace": ACE richness estimator
            - "fisher": Fisher's alpha
            - "evenness": Pielou's evenness
            - "all": All available metrics
        rarefy_depth (int): Depth for rarefaction before calculating diversity
            Set to None to use raw counts
    
    Returns:
        pd.DataFrame: Alpha diversity metrics for each sample
    """
    if metrics is None or metrics == "all":
        metrics = ["observed", "shannon", "simpson", "invsimpson", 
                  "chao1", "fisher", "evenness"]
    
    # Ensure data is an EcologicalData object
    if not isinstance(data, EcologicalData):
        raise TypeError("data must be an EcologicalData object")
    
    # Rarefy if specified
    if rarefy_depth is not None:
        data = data.rarefy(depth=rarefy_depth)
    
    # Initialize results DataFrame
    results = pd.DataFrame(index=data.abundance.columns)
    
    # Calculate each metric
    if "observed" in metrics:
        results["observed"] = (data.abundance > 0).sum(axis=0)
    
    if "shannon" in metrics:
        results["shannon"] = data.abundance.apply(lambda x: alpha.shannon(x), axis=0)
    
    if "simpson" in metrics:
        results["simpson"] = data.abundance.apply(lambda x: alpha.simpson(x), axis=0)
    
    if "invsimpson" in metrics:
        results["invsimpson"] = data.abundance.apply(lambda x: 1 / alpha.simpson(x) 
                                                  if alpha.simpson(x) > 0 else np.nan, axis=0)
    
    if "chao1" in metrics:
        results["chao1"] = data.abundance.apply(lambda x: alpha.chao1(x), axis=0)
    
    if "fisher" in metrics:
        # Calculate Fisher's alpha
        def fisher_alpha(counts):
            S = np.sum(counts > 0)  # Number of species
            N = np.sum(counts)      # Total number of individuals
            
            if S <= 1 or N <= 1:
                return np.nan
            
            # Estimate Fisher's alpha using iteration
            alpha = 0.1  # Initial guess
            step = 0.1
            
            # Iteratively improve the estimate
            for _ in range(50):
                S_est = alpha * np.log(1 + N / alpha)
                if abs(S_est - S) < 0.001:
                    break
                
                if S_est < S:
                    alpha += step
                else:
                    alpha -= step
                    step /= 2
                    alpha += step
            
            return alpha
        
        results["fisher"] = data.abundance.apply(fisher_alpha, axis=0)
    
    if "evenness" in metrics:
        # Pielou's evenness: Shannon diversity / log(richness)
        shannon = data.abundance.apply(lambda x: alpha.shannon(x), axis=0)
        richness = (data.abundance > 0).sum(axis=0)
        results["evenness"] = shannon / np.log(richness)
        # Replace inf values with NaN
        results["evenness"].replace([np.inf, -np.inf], np.nan, inplace=True)
    
    # Add sample metadata if available
    if data.metadata is not None:
        results = pd.merge(results, data.metadata, 
                         left_index=True, right_index=True, 
                         how="left")
    
    return results


def alpha_diversity_stats(alpha_div, group_col, metric, test="anova"):
    """
    Calculate statistical tests for alpha diversity metrics between groups.
    
    Parameters:
        alpha_div (pd.DataFrame): Alpha diversity metrics from calculate_alpha_diversity
        group_col (str): Column name in alpha_div for grouping samples
        metric (str): Alpha diversity metric to test
        test (str): Statistical test to perform:
            - "anova": One-way ANOVA
            - "kruskal": Kruskal-Wallis test (non-parametric)
            - "t.test": Two-sample t-test (for exactly two groups)
            - "wilcox": Wilcoxon rank-sum test (for exactly two groups)
    
    Returns:
        dict: Results of the statistical test
    """
    # Check if metric exists in alpha_div
    if metric not in alpha_div.columns:
        raise ValueError(f"Metric '{metric}' not found in alpha diversity DataFrame")
    
    # Check if group_col exists in alpha_div
    if group_col not in alpha_div.columns:
        raise ValueError(f"Group column '{group_col}' not found in alpha diversity DataFrame")
    
    # Get groups
    groups = alpha_div[group_col].unique()
    
    # Perform the appropriate test
    if test == "anova":
        # One-way ANOVA
        group_data = [alpha_div[alpha_div[group_col] == g][metric].dropna() for g in groups]
        result = stats.f_oneway(*group_data)
        return {
            "test": "One-way ANOVA",
            "statistic": result.statistic,
            "p_value": result.pvalue,
            "groups": list(groups)
        }
    
    elif test == "kruskal":
        # Kruskal-Wallis test
        group_data = [alpha_div[alpha_div[group_col] == g][metric].dropna() for g in groups]
        result = stats.kruskal(*group_data)
        return {
            "test": "Kruskal-Wallis test",
            "statistic": result.statistic,
            "p_value": result.pvalue,
            "groups": list(groups)
        }
    
    elif test == "t.test":
        # Two-sample t-test
        if len(groups) != 2:
            raise ValueError("t.test requires exactly two groups")
        
        group1 = alpha_div[alpha_div[group_col] == groups[0]][metric].dropna()
        group2 = alpha_div[alpha_div[group_col] == groups[1]][metric].dropna()
        
        result = stats.ttest_ind(group1, group2)
        return {
            "test": "Two-sample t-test",
            "statistic": result.statistic,
            "p_value": result.pvalue,
            "groups": list(groups)
        }
    
    elif test == "wilcox":
        # Wilcoxon rank-sum test
        if len(groups) != 2:
            raise ValueError("wilcox test requires exactly two groups")
        
        group1 = alpha_div[alpha_div[group_col] == groups[0]][metric].dropna()
        group2 = alpha_div[alpha_div[group_col] == groups[1]][metric].dropna()
        
        result = stats.mannwhitneyu(group1, group2)
        return {
            "test": "Wilcoxon rank-sum test",
            "statistic": result.statistic,
            "p_value": result.pvalue,
            "groups": list(groups)
        }
    
    else:
        raise ValueError(f"Unknown test: {test}")