#!/usr/bin/env python3
# -*- coding: utf-8 -*-
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
        
        self.actual_reagent_quantity = reagent_quantity
        molecular_species_history = []
        molecular_species_history.append(reagent_quantity)
        self.molecular_species_history = molecular_species_history
        self.state_change_vector = state_change_vectors
        self.reactions_combinatorics = combinatorics
        if time:
            self.total_time = time
        self.actual_time = 0
        self.timestep_list = []
        self.timestep_list.append(self.actual_time)
        self.max_iteration = iteration
        self.actual_iteration = 0
        
        
        while self.actual_iteration < self.max_iteration:
            
            propensity_function_list = stochastic_backend.calculate_propensity_funct(reag_quant=self.actual_reagent_quantity, combinatorics=self.reactions_combinatorics)
            if 0 in set(propensity_function_list):
                print('A reagent reached 0')
                break                
            cumulative_propensity = sum(propensity_function_list)
            if cumulative_propensity == 0:
                break
            tau = stochastic_backend.calculate_tau(cumulative_propensity)
            self.actual_time += tau
            self.timestep_list.append(self.actual_time)
            mu = stochastic_backend.calculate_mu(propensity_list=propensity_function_list, cumulative_propensity=cumulative_propensity)
            self.actual_reagent_quantity = [i+e for i,e in zip(self.actual_reagent_quantity,self.state_change_vector[mu])]
            if set_fixed_reagents:
                for index in set_fixed_reagents:
                    self.actual_reagent_quantity[index] = self.molecular_species_history[0][index]
            self.molecular_species_history.append(self.actual_reagent_quantity)
            self.actual_iteration += 1
            
            
            
            