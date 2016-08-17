from __future__ import print_function
import numpy as np
from ase import Atoms
from h2o_pjt2_calculator import H2o_pjt2_calculator, get_internal_coords_h2o, get_cart_coords_h2o
from ase.optimize import BFGS
from ase.optimize.bfgslinesearch import BFGSLineSearch
from ase.vibrations import Vibrations
from ase.units import Bohr, Ha

#r_OH1 = 0.95792059
#r_OH2 = r_OH1
#theta = 104.4996469200000 * np.pi / 180.0
r_OH1 = 0.99792059
r_OH2 = 1.90792059
theta = 114.4996469200000 * np.pi / 180.0

r_o, r_h1, r_h2 = get_cart_coords_h2o(r_OH1, r_OH2, theta, [0.0,0.0, 0.0])

atoms = Atoms('OHH', positions=[r_o, r_h1, r_h2])

calc  = H2o_pjt2_calculator()
atoms.set_calculator(calc)

E = atoms.get_potential_energy()

print(E)

# optimize geometry
dyn = BFGSLineSearch(atoms, logfile='h2o_opt.log')
dyn.run(fmax=0.00001)

pos = atoms.get_positions()

q1,q2,theta = get_internal_coords_h2o(pos[0], pos[1], pos[2])
coords = [q1, q2, theta * 180/np.pi]

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

