from __future__ import print_function

"""This module defines an ASE interface to a potential for a single water molecule:

C     Potential PJT2 due Polyansky, Jensen and Tennyson,                                                                                                                                                           C     J. Chem. Phys., 105, 6490-6497 (1996)                                                                                                                                                                        C     Update of Polyansky, Jensen and Tennyson, J Chem Phys 101, 7651 (1994)) 

"""

import numpy as np
from ase.units import Bohr, Ha
import ase.data
from ase.calculators.calculator import Calculator
import h2o_pjt2

class H2o_pjt2_calculator(Calculator):
    """ An ASE interface to the PJT2 potential for a single water molecule"""
    
    implemented_properties = [
        'energy',
        'forces']

    def calculate(self, atoms=None, properties=['energy'],
                  system_changes=[]):

        Calculator.calculate(self, atoms, properties, system_changes)

        # only one property allowed
        assert(len(properties) == 1)
        
        # make sure we have a water molecule, and identify O, H1, H2
        chemical_formula = self.atoms.get_chemical_formula(mode='all')
        
        if chemical_formula == 'OHH':
            o = 0
            h1= 1
            h2 =2
        elif chemical_formula == 'HOH':
            o = 1
            h1= 0
            h2 =2
        elif chemical_formula == 'HHO':
            o = 2
            h1= 0
            h2 =1
        else:
            raise RuntimeError("chemical formula '{}' is not a water molecule!".format(chemical_formula))
        
        # compute internal coordinates, both stretch coordinates q1, q2 in [bohr] and theta in radians  
        pos=atoms.get_positions()
        q1_tmp = (pos[h1] -  pos[o]) / Bohr
        q2_tmp = (pos[h2] -  pos[o]) / Bohr

        q1 = np.linalg.norm(q1_tmp)
        q2 = np.linalg.norm(q2_tmp)
        
        q1_norm = q1_tmp / q1 
        q2_norm = q2_tmp / q2
        q_norm_sum = q1_norm + q2_norm  
        q_norm_sum_norm = q_norm_sum / np.linalg.norm(q_norm_sum)

        theta = 2.0 * abs(np.arcsin(np.linalg.norm(np.cross(q1_norm,q_norm_sum_norm))))

        if properties[0] == 'energy':

            V =0.0
            V= h2o_pjt2.pots(q1,q2,theta)
            
            self.results['energy'] = V * Ha

        if properties[0] == 'forces':
            self.results['forces'] = self.calculate_numerical_forces(atoms, d=0.00000001)

        
