from math_stats import *
from operator import itemgetter
import random


def select_mapping(mappings):
    """
    Normalises mapping log probabilities and selects one mapping location stochastically
    Note: it modifies the logs in the main list and converts them to normalised probabilities
    """
    # Sort mappings based on likelihood
    mappings.sort(key=itemgetter(-1))
    # Calculating log sum
    log_probs = [mapping[-1] for mapping in mappings]
    log_probs_total = logSum(log_probs)
    for mapping in mappings:
        # Normalizing log likelihood
        mapping[-1] = math.exp(mapping[-1] - log_probs_total)

    rand_num = random.uniform(0, 1)
    for mapping in mappings[:-1]:
        if rand_num < mapping[-1]:
            return mapping
    return mappings[-1]


def best_mapping(mapping_probs):
    """
     Selecting the best mapping among candidate locations for a multi-read after last iteration
    """
    mapping_probs = sorted(mapping_probs, reverse=True)
    selected_prob = mapping_probs[0]
    last_tie_index = 0

    # All max probabilities appear at the beginning of the list; find the index of last one
    for mapping in mapping_probs[1:]:
        # If it has the same probability as highest probability
        if mapping[0] == selected_prob[0]:
            last_tie_index += 1
        else:
            break
    # If there are multiple locations with the same probability
    if last_tie_index > 0:
        selected_prob = mapping_probs[random.randrange(last_tie_index + 1)]

    # For probability normalisation
    log_probs = [mapping[0] for mapping in mapping_probs]
    log_probs_total = logSum(log_probs)

    selected_prob[0] = math.exp(selected_prob[0] - log_probs_total)

    return selected_prob


def select_final_mapping(mappings):
    """
    Returns the most probable mapping location after all runs
    """
    # Sort mappings based on likelihood
    mappings.sort(key=itemgetter(-2), reverse=True)
    selected_mapping = mappings[0]
    last_tie_index = 0

    # All max probabilities appear at the beginning of the list; find the index of last one
    for mapping in mappings[1:]:
        # New idea!
        # If it has the same probability as highest probability
        if mapping[-2] > (selected_mapping[-2] - 0.02):
            last_tie_index += 1
        else:
            break
    # If there are multiple locations with the same probability
    if last_tie_index > 0:
        selected_mapping = mappings[random.randrange(last_tie_index + 1)]
    return selected_mapping


def select_final_mapping_stochastic(mapping_probs):
    """
    Returns a mapping location after all runs stochastically based on normalised probabilities
    """
    mapping_probs_sorted = sorted(mapping_probs)
    rand_num = random.uniform(0, 1)
    for mapping in mapping_probs_sorted[:-1]:
        if rand_num < mapping[0]:
            return mapping
    return mapping_probs_sorted[-1]
