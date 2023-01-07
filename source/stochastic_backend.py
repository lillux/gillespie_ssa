import numpy as np


def calculate_propensity_funct(reag_quant:list, combinatorics:list) -> list:
    '''
    Calculate the propensity for each reaction to happen

    Parameters
    ----------
    reag_quant : list(int)
        Define the number of molecules for each reagent.
        If we have reagent [a, b, c], each of wich is present in the system with 10 molecules,
        our reag_quant = [10,10,10]
    combinatorics : list
        Define the combinatorial function for each reaction's reagent.
        Use a list of lambda functions to describe the combinatorial rules.
        A description of the combinatorial rules to be used here is found at
        page 11, equations 14a - 14g, of the 1976 paper from Daniel T. Gillespie:
            A General Method for Numerically Simulating 
            the Stochastic Time Evolution of Coupled Chemical Reactions.
        Here are the rule from Gillespie paper, where S is the chemical species,
        i,j,k are species identifiers,
        X is the number of molecule, for the respective identifier.
        If a reaction constants, q, is needed, it can be used to multiply the combinatorial function (hu * q).
            
        *               -> reaction products,                           hu = 1
        Sj              -> reaction products,                           hu = Xj
        Sj + Sk,        -> reaction products, with j != k               hu = Xj*Xi
        2Sj             -> reaction products,                           hu = Xj*(Xj-a)/2
        Si + Sj + Sk    -> reaction products, with i != j != k != i     hu = Xi*Xj*Xk
        Sj + 2Sk        -> reaction products, with j != k               hu = Xj*Xk*(Xk-1)/2
        3Sj             -> reaction products.                           hu = Xj*(Xj-1)*(Xj-2)/6

    Returns
    -------
    list(float)
        Return the propensity for each reaction to happen at the current time point.

    '''
    propensity_list = []
    for funct in combinatorics:
        prop = funct(*reag_quant)
        propensity_list.append(prop)
    return propensity_list

def calculate_tau(cumulative_propensity:float) -> float:
    '''
    Calculate the length of the next timestep.

    Parameters
    ----------
    cumulative_propensity : float
        The sum of the propensity of each reaction to happen.

    Returns
    -------
    float
        A float that define when the next reaction will happen in the system.

    '''
    r1 = np.random.uniform(0,1)
    tau = (1/cumulative_propensity)*(np.log(1/r1))
    return tau

def calculate_mu(propensity_list:list, cumulative_propensity:float) -> int:
    '''
    Calculate wich reaction will happen next

    Parameters
    ----------
    propensity_list : list(float)
        The propensity of each reaction to happen.
    cumulative_propensity : float
        The sum of the propensity of each reaction, that is the sum of the parameter propensity_list.

    Returns
    -------
    int
        The index of the next reaction.

    '''
    r2 = np.random.uniform(0,1)
    propensity_progressive_sum = 0
    threshold = r2*cumulative_propensity
    for index, reaction_propensity in enumerate(propensity_list):
        propensity_progressive_sum += reaction_propensity
        if propensity_progressive_sum > threshold:
            return index
        else:
            pass