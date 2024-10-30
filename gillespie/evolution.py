#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 30 15:40:13 2024

@author: lillux
"""

import numpy as np
from joblib import Parallel, delayed

def update_gen(wt: int, sv: int, r: float, s: float, u: float) -> tuple[int, int]:
    '''
    Perform one step of the evolutionary simulation.

    Parameters
    ----------
    wt : int
        The number of Wild Type individuals.
    sv : int
        The number of standing variation individuals (mutants).
    r : float
        The impact of an environmental stressor that reduce the fitness. Must be 0 < r < 1.
    s : float
        The impact of the mutation, that improve the fitness. Must be s > r
    u : float
        Probability tha a mutation appear in the next generation.

    Returns
    -------
    tuple[int, int]
        Return the updated values respectively for the WT and the SV populations.

    '''
    # Update fitness values
    fit_wt = 1 - r
    fit_sv = fit_wt * (1 + s)

    # Calculate new Standing Variation population
    if sv == 0:
        sv_new = 0
    else:
        # Total mean for Poisson distribution
        sv_mean = sv * fit_sv
        sv_new = np.random.poisson(lam=sv_mean)

    # Calculate new Wild Type population
    if wt == 0:
        wt_new = 0
    else:
        wt_mean = wt * fit_wt
        wt_new = np.random.poisson(lam=wt_mean)

        if wt_new > 0:
            # Mutations in Wild Type population
            mut_new = np.random.binomial(n=wt_new, p=u)
            wt_new -= mut_new
            sv_new += mut_new

    return wt_new, sv_new
    

class EvoSim:
    

    def __init__(self, wt: int, sv: int, r: float, s: float, u: float, max_gen: int):
        self.wt_start = wt
        self.sv_start = sv
        self.r = r
        self.s = s
        self.u = u
        self.max_gen = max_gen

        # instantiate array to hold generations data
        self.tot_array = np.zeros(max_gen, dtype=int)
        self.wt_array = np.zeros(max_gen, dtype=int)
        self.sv_array = np.zeros(max_gen, dtype=int)
        self.tot_array[0] = self.wt_start + self.sv_start
        self.wt_array[0] = self.wt_start
        self.sv_array[0] = self.sv_start

        # Instantiate outcome variables
        self.exinct = False
        self.saved = False

        self.wt = self.wt_start
        self.sv = self.sv_start
        # run population evolution
        for gen in range(1, self.max_gen):

            self.wt, self.sv = update_gen(self.wt, self.sv, self.r, self.s, self.u)
            tot = self.wt + self.sv
            self.tot_array[gen] = tot
            self.wt_array[gen] = self.wt
            self.sv_array[gen] = self.sv
            
            if tot == 0:
                self.exinct = True
                break
            if self.sv > self.wt_start:
                self.saved = True
                break
        # trim results arrays to the number of generations
        self.tot_array = self.tot_array[:gen+1]
        self.wt_array = self.wt_array[:gen+1]
        self.sv_array = self.sv_array[:gen+1]
        self.generations = gen+1

        return
    
    
class MultiSim:

    def __init__(self, epochs:int, params:dict[str:float], ncore=-1, replicable=False):

        # Prepare the list of arguments
        self.epochs = epochs
        self.params = params
        simulation_args = [(epoch, self.params, replicable) for epoch in range(self.epochs)]

        
        def run_simulation(epoch, params, replicable=False):
            
            if replicable:
                np.random.seed(epoch)  # Optional: Set seed for reproducibility
            
            evo = EvoSim(
                wt=params['wt'],
                sv=params['sv'],
                r=params['r'],
                s=params['s'],
                u=params['u'],
                max_gen=params['max_gen']
            )
            return (epoch, evo)

                # Run simulations
        self.results = Parallel(n_jobs=ncore)(
            delayed(run_simulation)(epoch, params, replicable) for epoch, params, replicable in simulation_args
        )

        # self.rescued = sum([i[1].saved for i in self.results])
        self.rescued = {i[0]:i[1] for i in self.results if i[1].saved}
        self.rescued_n = len(self.rescued)
        return