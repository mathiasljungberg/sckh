from __future__ import print_function
import numpy as np
from ase import Atoms
from ase.optimize.optimize import Dynamics
from copy import copy
import ase.units as units
import weakref

class Mc(Dynamics):
    """Simple Metropolis Monte Carlo class. All property calculations should be done with observers, 
    preferrably inherited from Mc_observer"""

    def __init__(self, atoms, max_step = 0.5, T= 298.0, accept_ratio=None):
        self.max_step = max_step
        self.beta = 1.0 / (T * units.kB)
        self.accept_ratio = accept_ratio

        Dynamics.__init__(self, atoms,logfile=None, trajectory=None)
        
    def run(self, steps=50):
        """Do Metropolis Monte Carlo run""" 

        self.E = self.atoms.get_potential_energy()
        self.E_new = 0.0
        
        for step in range(steps):
            self.atoms_old = copy(self.atoms)
            self.mc_step()
            self.nsteps += 1
            self.call_observers()
            
    def mc_step(self):
        # do a random move for a random coordinate
        a = int(np.random.random() * len(self.atoms))
        coord = int(np.random.random() * 3)
        pos= self.atoms.get_positions()
        pos_old = copy(pos)
        h = 2.0*(np.random.random()-0.5) * self.max_step
        pos[a][coord] += h
        self.atoms.set_positions(pos)

        # get new energy
        self.E_new = self.atoms.get_potential_energy()

        # test for acceptance
        E_diff = self.E_new - self.E

        #print(E_diff)
        
        if(np.random.random() < np.exp(-self.beta * E_diff)):
            self.accept =True
        else:
            self.accept= False
            
        if self.accept:
            self.E = self.E_new
        else:
            self.atoms.set_positions(pos_old)

# base observer class, only counts the number of accepted moves 
class Mc_observer:
    def __init__(self, dyn, atoms):
        self.dyn = weakref.proxy(dyn)
        self.atoms = atoms
        self.natoms = atoms.get_number_of_atoms()  

        self.n_steps = 0.0
        self.n_acc = 0.0
        
    def __call__(self):
        self.n_steps += 1.0
        if self.dyn.accept:
            self.n_acc += 1.0
