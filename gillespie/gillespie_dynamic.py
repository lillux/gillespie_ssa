#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov  3 12:09:21 2024

@author: lillux
"""
import numpy as np
import time
from gillespie.stochastic_backend import calculate_propensity_funct, calculate_mu, calculate_tau
from typing import Dict, List
from collections.abc import Callable
# import warnings

class gillespie_dynamic():
    
    def __init__(self,
                 reagent_quantity: List[int],
                 state_change_vectors: Dict[str:List[List[int]]],
                 combinatorics: Dict[str:List[List[Callable[[List[int]], int]]]],
                 max_time: float = None,
                 max_iteration: int = None,
                 stop_condition: str = 'time',
                 set_fixed_reagents: List = None,
                 rescale: int = None,
                 Ni: int = None,
                 oscillate: bool = False,
                 oscillation_interval: Dict[str:float] = None,
                 start_with: str = None):
        '''
        Initialize gillespie stochastic simulation algorithm class

        Parameters
        ----------
        reagent_quantity : List[int]
            The list of number (int) of molecule for each molecular species involved in the reaction
        state_change_vectors : Dict[str:List[List[int]]]
            Dict containing the state change vectors for each reaction.
        combinatorics : Dict[str:List[List[Callable[[List[int]], int]]]]
            Dict containing combinatorial functions for each reaction.
        max_time : float, optional
            Maximum simulation time. Required if stop_condition='time'.
        max_iteration : int, optional
            Maximum number of iterations. Required if stop_condition='iterations'.
        stop_condition : str, optional
            Specify 'time' or 'iterations' as the stopping condition. Default is 'time'.
        set_fixed_reagents : List, optional
            List of indices of reagents to keep fixed throughout the simulation.
            The default is False.
        rescale : int, optional
            Threshold for population rescaling.ù
        Ni : int, optional
            Population size after rescaling.
        oscillate : bool, optional
            Whether to oscillate between states.
        oscillation_interval : Dict[str:float], optional
            Time intervals for state oscillations.
        start_with : str, optional
            Initial state to start the simulation with.
        Returns
        -------
        None.

        '''
        
        self.actual_reagent_quantity = reagent_quantity.copy()
        self.molecular_species_history = [reagent_quantity.copy()]
        self.state_change_vector = state_change_vectors
        self.reactions_combinatorics = combinatorics
        self.set_fixed_reagents = set_fixed_reagents
        self.rescale = rescale
        self.Ni = Ni
        self.oscillate = oscillate
        self.oscillation_interval = oscillation_interval
        self.start_with = start_with
        self.stop_condition = stop_condition
        self.max_time = max_time
        self.actual_time = 0
        self.timestep_list = [self.actual_time]
        self.max_iteration = max_iteration
        self.actual_iteration = 0

        # Validate stop_condition and required parameters
        if self.stop_condition == 'time':
            if self.max_time is None:
                raise ValueError("max_time must be specified when stop_condition is 'time'")
        elif self.stop_condition == 'iterations':
            if self.max_iteration is None:
                raise ValueError("max_iteration must be specified when stop_condition is 'iterations'")
        else:
            raise ValueError("Invalid stop_condition. Choose 'time' or 'iterations'.")

        # Initialize state tracking
        state_time = 0
        running_state = self.start_with
        self.time_tracker = {i:{'start':[], 'end':[]} for i in self.reactions_combinatorics.keys()}
        self.time_tracker[running_state]['start'].append(self.actual_time)
        if self.oscillation_interval is None:
            self.oscillation_interval = {}
        
        # Set loop condition based on stopping criterion
        if self.stop_condition == 'time':
            loop_condition = lambda: self.actual_time < self.max_time
        elif self.stop_condition == 'iterations':
            loop_condition = lambda: self.actual_iteration < self.max_iteration

        # while self.actual_iteration < self.max_iteration:
        while loop_condition():
            
            # check if there are still reagents
            if np.sum(self.actual_reagent_quantity) <= 0:
                self.time_tracker[running_state]['end'].append(self.actual_time)
                print('No reagent left, stopping simulation.')
                break
            
            # switch state if time interval has passed
            if self.oscillation_interval and state_time > self.oscillation_interval[running_state]:
                self.time_tracker[running_state]['end'].append(self.actual_time)
                print(state_time)
                state_time = 0
                # code below works only if you have 2 states
                running_state = list(set(self.reactions_combinatorics.keys()) - set([running_state]))[0]
                self.time_tracker[running_state]['start'].append(self.actual_time)
                
            # Print progress each million of iteration
            if (self.actual_iteration % 1000000) == 0:
                print(f'Actual state is: {running_state}, actual iteration is: {self.actual_iteration}, simulation time is {self.actual_time}, real time is: {time.ctime()}.')

                        
            # calculate propensity function for each reaction
            propensity_function_list = calculate_propensity_funct(
                reag_quant=self.actual_reagent_quantity, 
                combinatorics=self.reactions_combinatorics[running_state])
              
            # calculate cumulative propensity
            cumulative_propensity = sum(propensity_function_list)
            # check if there are reaction that can happen
            # this break point almost never comes in, because the one above comes first
            if cumulative_propensity <= 0:
                self.time_tracker[running_state]['end'].append(self.actual_time)
                print('Cumulative propensity is zero, stopping simulation.')
                break
            
            # calculate next timestep
            tau = calculate_tau(cumulative_propensity)
            state_time += tau
            self.actual_time += tau
            self.timestep_list.append(self.actual_time)
            
            # calculate next reaction
            mu = calculate_mu(propensity_list=propensity_function_list, cumulative_propensity=cumulative_propensity)
            # print(f'Combinatoric is {running_state}, mu is {mu}, cumcumulative_propensity is {cumulative_propensity}')
            self.actual_reagent_quantity = [i+e for i,e in zip(self.actual_reagent_quantity,self.state_change_vector[running_state][mu])]
            
            # Ensure reagent quantities are non-negative
            # assert np.all(np.array(self.actual_reagent_quantity) > 0), f'Reagent with index {np.where(np.array(self.actual_reagent_quantity) <= 0)} are zero or negative.'
            # self.actual_reagent_quantity = [max(0, qty) for qty in self.actual_reagent_quantity]

            # check if some reagents have to be fixed
            if self.set_fixed_reagents:
                for index in self.set_fixed_reagents:
                    self.actual_reagent_quantity[index] = self.molecular_species_history[0][index]
            if self.rescale and np.sum(self.actual_reagent_quantity) > self.rescale:
                scale_factor = self.Ni / self.rescale
                self.actual_reagent_quantity = [np.random.poisson(reag * scale_factor) for reag in self.actual_reagent_quantity]
        
            # update reagent quantities
            self.molecular_species_history.append(self.actual_reagent_quantity)
            
            # update iteration counter
            self.actual_iteration += 1
        self.time_tracker[running_state]['end'].append(self.actual_time)
        return