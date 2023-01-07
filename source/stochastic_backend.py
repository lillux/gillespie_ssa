#import torch
import numpy as np
#import matplotlib.pyplot as plt

def calculate_propensity_funct(reag_quant : list, combinatorics : list) -> list:
    propensity_list = []
    for funct in combinatorics:
        prop = funct(*reag_quant)
        propensity_list.append(prop)
    return propensity_list

def calculate_tau(cumulative_propensity: float) -> float:
    r1 = np.random.uniform(0,1)
    tau = (1/cumulative_propensity)*(np.log(1/r1))
    return tau

def calculate_mu(propensity_list: list, cumulative_propensity: float) -> int:
    r2 = np.random.uniform(0,1)
    propensity_progressive_sum = 0
    threshold = r2*cumulative_propensity
    for index, reaction_propensity in enumerate(propensity_list):
        propensity_progressive_sum += reaction_propensity
        if propensity_progressive_sum > threshold:
            return index
        else:
            pass