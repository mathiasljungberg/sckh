from __future__ import print_function
import numpy as np
from ase import Atoms
from h2o_pjt2_calculator import H2o_pjt2_calculator
from ase.optimize import BFGS
from ase.optimize.bfgslinesearch import BFGSLineSearch
from ase.vibrations import Vibrations

#r_OH1 = 0.95792059
#r_OH2 = r_OH1
#theta = 104.4996469200000 * np.pi / 180.0
r_OH1 = 0.99792059
r_OH2 = 1.90792059
theta = 114.4996469200000 * np.pi / 180.0



print(theta)
atoms = Atoms('OHH', positions=[(0.0, 0.0, 0.0),
                              (np.sin(theta/2.0) * r_OH1, np.cos(theta/2.0) * r_OH1, 0.0 ),
                              (-np.sin(theta/2.0) * r_OH2, np.cos(theta/2.0) * r_OH2, 0.0 )])


calc  = H2o_pjt2_calculator()
atoms.set_calculator(calc)

E = atoms.get_potential_energy()

print(E)

# optimize geometry
dyn = BFGSLineSearch(atoms, logfile='h2o_opt.log')
dyn.run(fmax=0.00001)

pos = atoms.get_positions()

q1_tmp = (pos[1] -  pos[0]) 
q2_tmp = (pos[2] -  pos[0]) 

q1 = np.linalg.norm(q1_tmp)
q2 = np.linalg.norm(q2_tmp)
        
q1_norm = q1_tmp / q1 
q2_norm = q2_tmp / q2
q_norm_sum = q1_norm + q2_norm  
q_norm_sum_norm = q_norm_sum / np.linalg.norm(q_norm_sum)

theta = 2.0 * abs((180.0 / np.pi)* np.arcsin(np.linalg.norm(np.cross(q1_norm,q_norm_sum_norm))))
        
coords = [q1, q2, theta]

ref=np.array([0.957920759532, 0.957920773357, 104.499646202])

error = np.sqrt(np.sum((coords - ref)**2))
print('coords')
print(coords)
print('diff from reference:')
print(error)


# compute frequencies
vib = Vibrations(atoms, delta=0.00001)
vib.run()

vib.summary(method='frederiksen')

