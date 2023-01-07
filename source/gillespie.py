"""
Created on Wed Dec 28 19:15:38 2022

@author: lillo
"""

import stochastic_backend

class gillespie_ssa():
    
    def __init__(self,
                 reagent_quantity:list,
                 state_change_vectors:list,
                 combinatorics,
                 time=False,
                 iteration:int = 100,
                 set_fixed_reagents=False):
        '''
        Initialize gillespie simulation class
        
        Parameters
        ----------
        reagent_quantity : list
            The list of number (int) of molecule for each molecular species involved in the reaction
        state_change_vectors : list
            The list of list, list(list(int)), containing the variation of reagent quantity for a certain chemical reaction.
        combinatorics : list(function)
            The list of lambda function, each of which describe the combinatorics for the reagent of a reaction
        iteration : int, optional
            The number of iteration of the algorithm to perform.
            The default is 100.
        set_fixed_reagents : list, optional
            A list(int), each int in the list is the index
            of a corresponding reagent (in reagent_quantity argument),
            that we want to mantain fixed at the initial value through the simulation.
            The default is False.

        Returns
        -------
        None.

        '''
        
        self.actual_reagent_quantity = reagent_quantity
        molecular_species_history = []
        molecular_species_history.append(reagent_quantity)
        self.molecular_species_history = molecular_species_history
        self.state_change_vector = state_change_vectors
        self.reactions_combinatorics = combinatorics
        self.actual_time = 0
        self.timestep_list = []
        self.timestep_list.append(self.actual_time)
        self.max_iteration = iteration
        self.actual_iteration = 0
        
        while self.actual_iteration < self.max_iteration:
            # calculate propensity function for each reaction
            propensity_function_list = stochastic_backend.calculate_propensity_funct(
                reag_quant=self.actual_reagent_quantity, 
                combinatorics=self.reactions_combinatorics)
            # check break points, break if a reagent goes to 0
            if 0 in set(propensity_function_list):
                print('A reagent reached 0')
                break                
            # calculate cumulative propensity
            cumulative_propensity = sum(propensity_function_list)
            # check if there are reaction that can happen
            # this break point almost never comes in, because the one above comes first
            if cumulative_propensity == 0:
                break
            # calculate next timestep
            tau = stochastic_backend.calculate_tau(cumulative_propensity)
            self.actual_time += tau
            self.timestep_list.append(self.actual_time)
            # calculate next reaction
            mu = stochastic_backend.calculate_mu(propensity_list=propensity_function_list, cumulative_propensity=cumulative_propensity)
            self.actual_reagent_quantity = [i+e for i,e in zip(self.actual_reagent_quantity,self.state_change_vector[mu])]
            # check if some reagents have to be fixed
            if set_fixed_reagents:
                for index in set_fixed_reagents:
                    self.actual_reagent_quantity[index] = self.molecular_species_history[0][index]
            # update reagent quantities
            self.molecular_species_history.append(self.actual_reagent_quantity)
            # update iteration counter
            self.actual_iteration += 1
        return